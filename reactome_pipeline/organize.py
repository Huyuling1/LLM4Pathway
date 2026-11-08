"""Index global chains by pathway; store references and read full chains on demand."""

import argparse
import json
import os
import re
import sqlite3
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path

from .common import digest, dumps, file_sha, read_json, read_jsonl, write_json, write_jsonl
from .graph import EDGE_RULES
from .fetch import load_snapshot
from .locking import SnapshotLock


def slug(name):
    value = re.sub(r"[^\w.-]+", "_", unicodedata.normalize("NFKC", name or "unrecorded"))
    return value.strip("._")[:80] or "unrecorded"


def view_locations(entry):
    species = entry.get("species") or [{"tax_id": None, "name": None}]
    result = []
    for organism in species:
        taxid = organism.get("tax_id")
        suffix = "taxon-" + slug(str(taxid)) if taxid else "unrecorded-" + digest(organism)[:12]
        species_folder = slug(organism.get("name")) + "__" + suffix
        pathway_folder = ("unassigned" if entry.get("unassigned") else slug(entry.get("display_name"))) + "__" + slug(entry["id"])
        result.append("by_species/" + species_folder + "/" + pathway_folder)
    return sorted(set(result))


def index_dataset(dataset):
    dataset = Path(dataset)
    root = dataset / "global"
    catalog_status = read_json(dataset / "catalog/status.json")
    if not catalog_status["complete"]:
        raise ValueError("Pathway catalog is unfinished")
    hierarchy = read_json(dataset / "pathway_hierarchy.json")
    snapshot = load_snapshot(root)
    if hierarchy["reaction_inventory_sha256"] != snapshot["inventory_sha256"]:
        raise ValueError("Hierarchy and reaction inventory differ")
    if (file_sha(dataset / "pathway_hierarchy.json") != catalog_status["hierarchy_sha256"] or
            file_sha(dataset / "catalog/pathway_membership.jsonl") != catalog_status["membership_sha256"]):
        raise ValueError("Frozen catalog checksum differs")
    graph = read_json(root / "graph/status.json")
    normalized_sha = file_sha(root / "normalized/reactions.jsonl")
    if digest({"normalization": normalized_sha, "currency": graph["currency_sha256"], "edge_rules": EDGE_RULES}) != graph["graph_sha256"]:
        raise ValueError("Normalized reactions and graph differ")
    pathways = hierarchy["pathways"]
    fingerprint = digest({"graph": graph["graph_sha256"], "hierarchy": catalog_status["hierarchy_sha256"],
                          "membership": catalog_status["membership_sha256"], "rules": "any-member-full-chain-v1"})
    target = dataset / "index.sqlite"
    if not target.exists():
        temp = dataset / "index.building.sqlite"
        if temp.exists():
            temp.unlink()
        db = sqlite3.connect(temp)
        db.executescript("""
            CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
            CREATE TABLE pathways(id TEXT PRIMARY KEY,data TEXT);
            CREATE TABLE members(pathway TEXT,reaction TEXT,PRIMARY KEY(pathway,reaction));
            CREATE INDEX members_reaction ON members(reaction,pathway);
            CREATE TABLE chain_refs(pathway TEXT,chain_id TEXT,ordinal INTEGER,PRIMARY KEY(pathway,chain_id));
            CREATE INDEX chain_refs_order ON chain_refs(pathway,ordinal);
        """)
        for row in read_jsonl(dataset / "catalog/pathway_membership.jsonl"):
            db.executemany("INSERT OR IGNORE INTO members VALUES (?,?)", ((row["pathway_id"], rid) for rid in row["reaction_ids"]))
        for reaction in read_jsonl(root / "normalized/reactions.jsonl"):
            if db.execute("SELECT 1 FROM members WHERE reaction=? LIMIT 1", (reaction["id"],)).fetchone():
                continue
            for organism in reaction.get("species") or [{"tax_id": None, "name": reaction.get("species_name")}]:
                taxid = organism.get("tax_id")
                pid = "unassigned-" + (str(taxid) if taxid else digest(organism)[:12])
                pathways.setdefault(pid, {"id": pid, "display_name": "Unassigned reactions", "unassigned": True,
                                         "species": [organism], "children": [], "membership_sources": []})
                db.execute("INSERT OR IGNORE INTO members VALUES (?,?)", (pid, reaction["id"]))
        db.executemany("INSERT INTO pathways VALUES (?,?)", ((pid, dumps(entry)) for pid, entry in sorted(pathways.items())))
        db.execute("INSERT INTO metadata VALUES ('fingerprint',?)", (fingerprint,))
        db.execute("INSERT INTO metadata VALUES ('last_ordinal','0')")
        db.commit()
        db.close()
        os.replace(temp, target)
    db = sqlite3.connect(target)
    try:
        if db.execute("SELECT value FROM metadata WHERE key='fingerprint'").fetchone()[0] != fingerprint:
            raise ValueError("Dataset index belongs to another graph/catalog; use a fresh export dataset")
        path_status_file = root / "paths/status.json"
        path_status = read_json(path_status_file) if path_status_file.exists() else {"enumeration_complete": False}
        if path_status.get("graph_sha256", graph["graph_sha256"]) != graph["graph_sha256"]:
            raise ValueError("Chains and graph differ")
        last = int(db.execute("SELECT value FROM metadata WHERE key='last_ordinal'").fetchone()[0])
        if path_status_file.exists():
            version = path_status["record_version"] + ":" + path_status.get("start_order", "stable-id")
            saved = db.execute("SELECT value FROM metadata WHERE key='chain_version'").fetchone()
            if saved and saved[0] != version:
                raise ValueError("Chain version/order changed; rebuild the view index")
            db.execute("INSERT OR REPLACE INTO metadata VALUES ('chain_version',?)", (version,))
            mapping = defaultdict(set)
            for pid, rid in db.execute("SELECT pathway,reaction FROM members"):
                mapping[rid].add(pid)
            with sqlite3.connect((root / "paths/paths.sqlite").resolve().as_uri() + "?mode=ro", uri=True) as chains:
                for ordinal, cid, data in chains.execute("SELECT ordinal,id,data FROM chains WHERE ordinal>? ORDER BY ordinal", (last,)):
                    row = json.loads(data)
                    memberships = set().union(*(mapping[rid] for rid in row["reaction_ids"]))
                    db.executemany("INSERT OR IGNORE INTO chain_refs VALUES (?,?,?)", ((pid, cid, ordinal) for pid in sorted(memberships)))
                    last = ordinal
                    db.execute("UPDATE metadata SET value=? WHERE key='last_ordinal'", (str(last),))
                    if ordinal % 1000 == 0:
                        db.commit()
        db.commit()
        def entries():
            for pid, data in db.execute("SELECT id,data FROM pathways ORDER BY id"):
                entry = json.loads(data)
                yield {"pathway_id": pid, "display_name": entry["display_name"], "locations": view_locations(entry),
                       "member_reactions": db.execute("SELECT count(*) FROM members WHERE pathway=?", (pid,)).fetchone()[0],
                       "available_chains": db.execute("SELECT count(*) FROM chain_refs WHERE pathway=?", (pid,)).fetchone()[0]}
        write_jsonl(dataset / "pathway_index.jsonl", entries())
        status = {"format_version": 1, "reactome_release": snapshot["reactome_release"],
                  "graph_sha256": graph["graph_sha256"],
                  "global_directory": "global", "pathway_directory": "by_species", "index_sha256": fingerprint,
                  "catalog_complete": True, "input_complete": graph["input_complete"],
                  "enumeration_complete": path_status["enumeration_complete"], "indexed_through_ordinal": last,
                  "chain_record_version": path_status.get("record_version"),
                  "pathways_and_unassigned_groups": db.execute("SELECT count(*) FROM pathways").fetchone()[0],
                  "membership_rule": "any member reaction; full global chain retained; descendant pathway membership included",
                  "materialization": "on demand; pathway_index.jsonl lists all available views"}
        write_json(dataset / "manifest.json", status)
        return status
    finally:
        db.close()


def pathway_chains(dataset, pathway):
    """Resolve indexed IDs to full global records without creating a copy."""
    dataset = Path(dataset)
    with sqlite3.connect((dataset / "index.sqlite").resolve().as_uri() + "?mode=ro", uri=True) as index:
        if not index.execute("SELECT 1 FROM pathways WHERE id=?", (pathway,)).fetchone():
            raise ValueError("Pathway not in index: " + pathway)
        members = {r[0] for r in index.execute("SELECT reaction FROM members WHERE pathway=?", (pathway,))}
        refs = index.execute("SELECT chain_id FROM chain_refs WHERE pathway=? ORDER BY ordinal", (pathway,))
        path_db = dataset / "global/paths/paths.sqlite"
        if not path_db.exists():
            if refs.fetchone():
                raise ValueError("Indexed chains exist but the global chain database is missing")
            return
        with sqlite3.connect(path_db.resolve().as_uri() + "?mode=ro", uri=True) as chains:
            for (cid,) in refs:
                record = chains.execute("SELECT data FROM chains WHERE id=?", (cid,)).fetchone()
                if record is None:
                    raise ValueError("Indexed chain missing from global records: " + cid)
                row = json.loads(record[0])
                yield row, [rid for rid in row["reaction_ids"] if rid in members]


def show_pathway(dataset, pathway, stream):
    count = 0
    for row, _ in pathway_chains(dataset, pathway):
        if count:
            stream.write("\n")
        stream.write(row["text"] + "\n")
        count += 1
    return count


def export_full_chains(dataset, pathway, out):
    """Explicit standalone copy for sharing; never part of normal pathway export."""
    out = Path(out)
    if out.resolve().is_relative_to(Path(dataset).resolve()):
        raise ValueError("Standalone full-chain export must be outside the dataset directory")
    out.mkdir(parents=True, exist_ok=True)
    with (out / "chains.txt.tmp").open("w", encoding="utf-8", newline="\n") as text:
        count = 0
        def rows():
            nonlocal count
            for row, members in pathway_chains(dataset, pathway):
                if count:
                    text.write("\n")
                text.write(row["text"] + "\n")
                count += 1
                yield {**row, "pathway_view": {"pathway_id": pathway, "member_reaction_ids": members}}
        write_jsonl(out / "chains.jsonl", rows())
    os.replace(out / "chains.txt.tmp", out / "chains.txt")
    return {"pathway_id": pathway, "chains": count, "output": str(out)}


def export_pathway(dataset, pathway):
    dataset = Path(dataset)
    root = dataset / "global"
    status = read_json(dataset / "manifest.json")
    db = sqlite3.connect((dataset / "index.sqlite").resolve().as_uri() + "?mode=ro", uri=True)
    try:
        row = db.execute("SELECT data FROM pathways WHERE id=?", (pathway,)).fetchone()
        if not row:
            raise ValueError("Pathway not in index: " + pathway)
        entry = json.loads(row[0])
        members = {r[0] for r in db.execute("SELECT reaction FROM members WHERE pathway=?", (pathway,))}
        graph = sqlite3.connect((root / "graph/graph.sqlite").resolve().as_uri() + "?mode=ro", uri=True)
        try:
            graph_sha = json.loads(graph.execute("SELECT value FROM metadata WHERE key='graph_sha256'").fetchone()[0])
            if graph_sha != status["graph_sha256"] or graph_sha != read_json(root / "graph/status.json")["graph_sha256"]:
                raise ValueError("Graph metadata mismatch")
            for location in view_locations(entry):
                out = dataset / location
                out.mkdir(parents=True, exist_ok=True)
                # Remove only unchanged full-text copies previously generated by us.
                old_manifest = read_json(out / "manifest.json") if (out / "manifest.json").exists() else {}
                obsolete = [out / name for name in ("chains.jsonl", "chains.txt") if (out / name).exists()]
                for old in obsolete:
                    if old_manifest.get("files", {}).get(old.name) != file_sha(old):
                        raise ValueError("Existing chain copy was modified or has no receipt; retain it separately: " + str(old))
                write_json(out / "manifest.json", {"pathway_id": pathway, "export_complete": False,
                                                   "files": old_manifest.get("files", {})})
                present = 0
                def reactions():
                    nonlocal present
                    for rid in sorted(members):
                        result = graph.execute("SELECT data FROM reactions WHERE id=?", (rid,)).fetchone()
                        if result:
                            present += 1
                            yield json.loads(result[0])
                write_jsonl(out / "reactions.jsonl", reactions())
                def edges():
                    for edge in read_jsonl(root / "graph/edges.jsonl"):
                        a, b = edge["source"] in members, edge["target"] in members
                        if a or b:
                            yield {**edge, "source_in_pathway": a, "target_in_pathway": b,
                                   "view_scope": "internal" if a and b else "boundary"}
                write_jsonl(out / "edges.jsonl", edges())
                count = 0
                def references():
                    nonlocal count
                    for row, matched in pathway_chains(dataset, pathway):
                        count += 1
                        yield {"chain_id": row["id"], "member_reaction_ids": matched}
                write_jsonl(out / "chain_refs.jsonl", references())
                for old in obsolete:
                    old.unlink()
                manifest = {**entry, "reactome_release": status["reactome_release"], "graph_sha256": graph_sha,
                            "global_directory": os.path.relpath(root, out),
                            "expected_member_reactions": len(members), "present_member_reactions": present,
                            "membership_complete": present == len(members), "chain_count": count,
                            "input_complete": status["input_complete"], "enumeration_complete": status["enumeration_complete"],
                            "indexed_through_ordinal": status["indexed_through_ordinal"], "export_complete": True,
                            "chain_storage": "references_only",
                            "chain_scope": "references to full global chains touching this pathway; external reactions retained on lookup",
                            "edge_scope": "internal and boundary edges; full graph in global_directory",
                            "files": {name: file_sha(out / name) for name in ("reactions.jsonl", "edges.jsonl", "chain_refs.jsonl")}}
                write_json(out / "manifest.json", manifest)
            return {"pathway_id": pathway, "locations": view_locations(entry), "chains": count}
        finally:
            graph.close()
    finally:
        db.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("index", "export", "show", "export-full"))
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--pathway", action="append", default=[])
    parser.add_argument("--out", help="Standalone full-chain export directory; required for export-full")
    args = parser.parse_args()
    if args.stage != "index" and not args.pathway:
        parser.error("This stage requires --pathway")
    if args.stage in ("show", "export-full") and len(args.pathway) != 1:
        parser.error("show and export-full require exactly one pathway")
    if args.stage == "export-full" and not args.out:
        parser.error("export-full requires --out")
    with SnapshotLock(Path(args.dataset) / "global"), SnapshotLock(Path(args.dataset) / "views"):
        print(index_dataset(args.dataset), file=sys.stderr if args.stage == "show" else sys.stdout)
        if args.stage == "export":
            if not args.pathway:
                parser.error("export requires --pathway")
            for pid in args.pathway:
                print(export_pathway(args.dataset, pid))
        elif args.stage == "show":
            show_pathway(args.dataset, args.pathway[0], sys.stdout)
        elif args.stage == "export-full":
            print(export_full_chains(args.dataset, args.pathway[0], args.out))
