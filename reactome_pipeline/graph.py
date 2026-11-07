"""Disk-backed exact-ID joins, with one evidence edge per entity and role."""

import json
import os
import sqlite3
from collections import defaultdict
from pathlib import Path

from .common import digest, dumps, file_sha, read_json, read_jsonl, write_json, write_jsonl
from .currency import load_currency


def build_graph(root, currency_path=None):
    root = Path(root)
    status = read_json(root / "normalized" / "status.json")
    config, config_sha, chebi_ids, explicit_ids = load_currency(currency_path)
    out = root / "graph"
    out.mkdir(parents=True, exist_ok=True)
    temp = out / "graph.building.sqlite"
    if temp.exists():
        temp.unlink()
    db = sqlite3.connect(temp)
    db.executescript("""
        CREATE TABLE reactions (id TEXT PRIMARY KEY, data TEXT NOT NULL);
        CREATE TABLE participants (reaction TEXT, role TEXT, kind TEXT, entity TEXT, pointer TEXT,
                                   PRIMARY KEY(reaction,role,kind,entity,pointer));
        CREATE TABLE entities (id TEXT PRIMARY KEY, name TEXT, chebi_ids TEXT, excluded INTEGER NOT NULL);
        CREATE TABLE edges (source TEXT,target TEXT,kind TEXT,entity TEXT,
                            PRIMARY KEY(source,target,kind,entity));
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT);
    """)
    entity_data, seen_reactions = {}, set()
    for reaction in read_jsonl(root / "normalized" / "reactions.jsonl"):
        rid = reaction["id"]
        if rid in seen_reactions:
            raise ValueError(f"Repeated reaction ID: {rid}")
        seen_reactions.add(rid)
        db.execute("INSERT INTO reactions VALUES (?,?)", (rid, dumps(reaction)))
        for role, kind in (("outputs", "output"), ("inputs", "material_flow"), ("catalysts", "catalyst_dependency"), ("regulators", None)):
            for entity in reaction[role]:
                eid = entity["entity_id"]
                if not eid:
                    continue
                info = entity_data.setdefault(eid, {"names": set(), "chebi_ids": set()})
                if entity["display_name"]:
                    info["names"].add(entity["display_name"])
                info["chebi_ids"].update(entity["chebi_ids"])
                db.execute("INSERT OR IGNORE INTO participants VALUES (?,?,?,?,?)", (rid, role, kind or entity["kind"], eid, entity["source_pointer"]))
    for eid, info in sorted(entity_data.items()):
        excluded = eid in explicit_ids or bool(info["chebi_ids"] & chebi_ids)
        db.execute("INSERT INTO entities VALUES (?,?,?,?)", (eid, sorted(info["names"])[0] if info["names"] else eid, dumps(sorted(info["chebi_ids"])), int(excluded)))
    db.executescript("""
        CREATE INDEX participants_entity ON participants(entity,role);
        CREATE INDEX participants_reaction ON participants(reaction);
        INSERT OR IGNORE INTO edges
        SELECT a.reaction,b.reaction,b.kind,a.entity
        FROM participants a JOIN participants b ON a.entity=b.entity
        JOIN entities e ON a.entity=e.id
        WHERE a.role='outputs' AND b.role!='outputs' AND e.excluded=0;
        CREATE INDEX edges_target ON edges(target);
    """)
    # Include self-loops in the graph; non-repetition is applied only during path enumeration.
    fingerprint = digest({"normalization": file_sha(root / "normalized" / "reactions.jsonl"), "currency": config_sha, "edge_rules": "exact-stid-v1"})
    metadata = {"format_version": 1, "graph_sha256": fingerprint, "currency_sha256": config_sha,
                "currency_version": config["version"], "input_complete": status["complete"],
                "reactions": len(seen_reactions), "evidence_edges": db.execute("SELECT count(*) FROM edges").fetchone()[0],
                "reaction_pairs": db.execute("SELECT count(*) FROM (SELECT DISTINCT source,target FROM edges)").fetchone()[0]}
    db.executemany("INSERT INTO metadata VALUES (?,?)", ((k, dumps(v)) for k, v in sorted(metadata.items())))
    db.commit()
    write_json(out / "currency.frozen.json", config)
    write_json(out / "status.json", metadata)
    write_jsonl(out / "edges.jsonl", iter_edges(db))
    def exclusions():
        for eid, name, chebi in db.execute("SELECT id,name,chebi_ids FROM entities WHERE excluded=1 ORDER BY id"):
            producers = db.execute("SELECT count(DISTINCT reaction) FROM participants WHERE entity=? AND role='outputs'", (eid,)).fetchone()[0]
            consumers = db.execute("SELECT count(DISTINCT reaction) FROM participants WHERE entity=? AND role!='outputs'", (eid,)).fetchone()[0]
            targets = db.execute("SELECT count(*) FROM (SELECT DISTINCT reaction,kind FROM participants WHERE entity=? AND role!='outputs')", (eid,)).fetchone()[0]
            candidates = producers * targets
            yield {"entity_id": eid, "display_name": name, "chebi_ids": json.loads(chebi), "producers": producers, "consumers": consumers, "excluded_evidence_edges": candidates}
    write_jsonl(out / "excluded_entities.jsonl", exclusions())
    db.close()
    os.replace(temp, out / "graph.sqlite")
    return metadata


def iter_edges(db):
    for source, target, kind, entity, name in db.execute("SELECT source,target,kind,entity,name FROM edges JOIN entities ON entity=id ORDER BY source,target,kind,entity"):
        evidence = {"source": source, "target": target, "kind": kind, "entity_id": entity}
        evidence["id"] = digest(evidence)
        evidence["entity_name"] = name
        # Pointers plus source records locate both roles in their original nested JSON.
        evidence["source_pointers"] = [x[0] for x in db.execute("SELECT DISTINCT pointer FROM participants WHERE reaction=? AND entity=? AND role='outputs' ORDER BY pointer", (source, entity))]
        evidence["target_pointers"] = [x[0] for x in db.execute("SELECT DISTINCT pointer FROM participants WHERE reaction=? AND entity=? AND kind=? AND role!='outputs' ORDER BY pointer", (target, entity, kind))]
        yield evidence


def load_graph(root):
    db = sqlite3.connect(Path(root) / "graph" / "graph.sqlite")
    reactions = {rid: json.loads(data) for rid, data in db.execute("SELECT id,data FROM reactions ORDER BY id")}
    adjacency = {rid: [] for rid in reactions}
    reverse = {rid: [] for rid in reactions}
    for source, target in db.execute("SELECT DISTINCT source,target FROM edges ORDER BY source,target"):
        adjacency[source].append(target)
        reverse[target].append(source)
    metadata = {k: json.loads(v) for k, v in db.execute("SELECT key,value FROM metadata")}
    return db, reactions, adjacency, reverse, metadata


def components(adjacency, reverse):
    """Iterative Kosaraju, avoiding Python recursion limits on long pathways."""
    seen, order = set(), []
    for start in sorted(adjacency):
        if start in seen:
            continue
        seen.add(start)
        stack = [(start, 0)]
        while stack:
            node, index = stack[-1]
            if index < len(adjacency[node]):
                child = adjacency[node][index]
                stack[-1] = node, index + 1
                if child not in seen:
                    seen.add(child)
                    stack.append((child, 0))
            else:
                order.append(node)
                stack.pop()
    seen, result = set(), []
    for start in reversed(order):
        if start in seen:
            continue
        seen.add(start)
        group, stack = [], [start]
        while stack:
            node = stack.pop()
            group.append(node)
            for parent in reverse[node]:
                if parent not in seen:
                    seen.add(parent)
                    stack.append(parent)
        result.append(sorted(group))
    return sorted(result, key=lambda group: group[0])


def weak_components(adjacency, reverse):
    """Independent graph components, ordered by size then stable identifier."""
    unseen = set(adjacency)
    groups = []
    for node in sorted(adjacency):
        if node not in unseen:
            continue
        unseen.remove(node)
        group, stack = [], [node]
        while stack:
            current = stack.pop()
            group.append(current)
            for neighbor in adjacency[current] + reverse[current]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        groups.append(sorted(group))
    return sorted(groups, key=lambda group: (len(group), group[0]))
