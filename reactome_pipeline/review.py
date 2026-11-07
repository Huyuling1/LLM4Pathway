"""Freeze a complete pathway case and export source-linked human review materials."""

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import sqlite3
from pathlib import Path

from .common import digest, file_sha, read_json, read_jsonl, write_bytes, write_json, write_jsonl
from .fetch import Client, QUERY, download, load_snapshot, raw_path
from .locking import SnapshotLock


def prepare(source, root, pathway, client):
    source, root = Path(source), Path(root)
    parent = load_snapshot(source)
    manifest_path = root / "snapshot.json"
    if manifest_path.exists():
        manifest = load_snapshot(root)
        if manifest.get("pathway_id") != pathway or manifest.get("parent_inventory_sha256") != parent["inventory_sha256"]:
            raise ValueError("Existing case belongs to another pathway or parent inventory")
        for receipt in manifest["membership_sources"]:
            if file_sha(root / receipt["file"]) != receipt["sha256"]:
                raise ValueError("Frozen pathway membership source was modified")
    else:
        if client.base != parent["base_url"] or client.version() != parent["reactome_release"]:
            raise ValueError("Pathway API and parent snapshot release differ")
        membership_sources = []
        objects = {}
        for name, endpoint in (("pathway", "/data/query/" + pathway),
                               ("contained", "/data/pathway/" + pathway + "/containedEvents")):
            body, obj = client.get(endpoint)
            filename = "raw/membership/" + name + ".json"
            write_bytes(root / filename, body)
            membership_sources.append({"file": filename, "url": client.url(endpoint),
                                       "sha256": hashlib.sha256(body).hexdigest()})
            objects[name] = obj
        if client.version() != parent["reactome_release"]:
            raise ValueError("API release changed while freezing pathway membership")
        pathway_obj, contained = objects["pathway"], objects["contained"]
        if pathway_obj.get("stId") != pathway or pathway_obj.get("schemaClass") != "Pathway" or not isinstance(contained, list):
            raise ValueError("Unexpected pathway response")
        # The top-level array includes numeric dbId references as well as objects.
        ids = {int(x["dbId"] if isinstance(x, dict) else x) for x in contained}
        all_reactions = {x["dbId"]: x for x in parent["inventory"]}
        rows = sorted((all_reactions[x] for x in ids & all_reactions.keys()), key=lambda x: (x.get("stId", ""), x["dbId"]))
        if not rows:
            raise ValueError("Pathway has no reactions in the complete parent inventory")
        manifest = {"format_version": 1, "reactome_release": parent["reactome_release"],
                    "base_url": parent["base_url"], "query": QUERY,
                    "scope": "complete contained ReactionLikeEvents of one pathway; internal audit graph",
                    "pathway_id": pathway, "pathway_name": pathway_obj["displayName"],
                    "expected_reactions": len(rows), "inventory": rows, "inventory_sha256": digest(rows),
                    "parent_inventory_sha256": parent["inventory_sha256"],
                    "contained_event_db_ids": sorted(ids), "non_reaction_contained_db_ids": sorted(ids - all_reactions.keys()),
                    "membership_sources": membership_sources}
        write_json(manifest_path, manifest)
    source_receipts = {}
    if (source / "downloads.sqlite").exists():
        with sqlite3.connect((source / "downloads.sqlite").resolve().as_uri() + "?mode=ro", uri=True) as db:
            source_receipts = dict(db.execute("SELECT db_id,sha256 FROM receipts"))
    for row in manifest["inventory"]:
        src, dst = raw_path(source, row["dbId"]), raw_path(root, row["dbId"])
        if src.exists() and not dst.exists():
            body = gzip.decompress(src.read_bytes())
            obj = json.loads(body)
            if obj.get("dbId") != row["dbId"] or obj.get("stId") != row.get("stId"):
                raise ValueError("Cached reaction identity mismatch")
            if row["dbId"] in source_receipts and hashlib.sha256(body).hexdigest() != source_receipts[row["dbId"]]:
                raise ValueError("Cached reaction differs from its receipt")
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
    return download(root, client, workers=2)


def pointer_value(raw, pointer):
    obj = raw
    for part in pointer.split("/")[1:]:
        key = part.replace("~1", "/").replace("~0", "~")
        obj = obj[int(key)] if isinstance(obj, list) else obj[key]
    return obj


def prepare_classes(source, root, client):
    """Freeze one actual reaction per inventory class; prefer verified cached records."""
    source, root = Path(source), Path(root)
    parent = load_snapshot(source)
    if (root / "snapshot.json").exists():
        manifest = load_snapshot(root)
        if manifest.get("selection_algorithm") != "one-per-class-cached-first-stid-v1" or manifest.get("parent_inventory_sha256") != parent["inventory_sha256"]:
            raise ValueError("Existing output uses another class selection or parent inventory")
    else:
        candidates = sorted(parent["inventory"], key=lambda x: (x.get("schemaClass", ""),
                            not raw_path(source, x["dbId"]).exists(), x.get("stId", ""), x["dbId"]))
        selected = {}
        for row in candidates:
            selected.setdefault(row["schemaClass"], row)
        rows = sorted(selected.values(), key=lambda x: (x.get("stId", ""), x["dbId"]))
        manifest = {"format_version": 1, "reactome_release": parent["reactome_release"], "base_url": parent["base_url"],
                    "query": QUERY, "scope": "one real reaction per ReactionLikeEvent subclass in the parent inventory",
                    "selection_algorithm": "one-per-class-cached-first-stid-v1", "inventory": rows,
                    "inventory_sha256": digest(rows), "expected_reactions": len(rows),
                    "parent_inventory_sha256": parent["inventory_sha256"],
                    "parent_schema_classes": sorted(selected)}
        write_json(root / "snapshot.json", manifest)
    receipts = {}
    if (source / "downloads.sqlite").exists():
        with sqlite3.connect((source / "downloads.sqlite").resolve().as_uri() + "?mode=ro", uri=True) as db:
            receipts = dict(db.execute("SELECT db_id,sha256 FROM receipts"))
    for row in manifest["inventory"]:
        src, dst = raw_path(source, row["dbId"]), raw_path(root, row["dbId"])
        if src.exists() and not dst.exists():
            body = gzip.decompress(src.read_bytes())
            if row["dbId"] in receipts and hashlib.sha256(body).hexdigest() != receipts[row["dbId"]]:
                raise ValueError("Cached class example differs from receipt")
            write_bytes(dst, src.read_bytes())
    return download(root, client, workers=2)


def freeze(root, destination):
    root, destination = Path(root), Path(destination)
    manifest = load_snapshot(root)
    files = ["snapshot.json", "download_status.json"]
    files += [x["file"] for x in manifest.get("membership_sources", [])]
    files += [raw_path(root, x["dbId"]).relative_to(root).as_posix() for x in manifest["inventory"]]
    hashes = {}
    for filename in sorted(files):
        src, dst = root / filename, destination / filename
        if dst.exists() and file_sha(dst) != file_sha(src):
            raise ValueError("Frozen export already contains different data: " + filename)
        write_bytes(dst, src.read_bytes())
        hashes[filename] = file_sha(dst)
    write_json(destination / "checksums.json", hashes)
    return {"complete": True, "files": len(files), "scope": manifest["scope"], "pathway_id": manifest.get("pathway_id")}


def packet(root, context=None):
    root = Path(root)
    manifest = load_snapshot(root)
    reactions = {r["id"]: r for r in read_jsonl(root / "normalized" / "reactions.jsonl")}
    checks, rows = [], []
    for rid, reaction in sorted(reactions.items()):
        raw = read_json(root / reaction["source"]["file"])
        reference_index, stack = {}, [raw]
        while stack:
            obj = stack.pop()
            if isinstance(obj, dict):
                if "@id" in obj:
                    reference_index[str(obj["@id"])] = obj
                stack.extend(obj.values())
            elif isinstance(obj, list):
                stack.extend(obj)
        for role in ("inputs", "outputs", "catalysts", "regulators"):
            for item in reaction[role]:
                value = pointer_value(raw, item["source_pointer"])
                if isinstance(value, dict) and "@ref" in value:
                    value = reference_index[str(value["@ref"])]
                expected = item["entity_id"]
                actual = value.get("stId") if isinstance(value, dict) else None
                checks.append({"reaction_id": rid, "role": role, "entity_id": expected,
                               "source_st_id": actual, "source_pointer": item["source_pointer"],
                               "status": "pass" if expected is not None and actual == expected else "needs_review"})
        rows.append({"reaction_id": rid, "display_name": reaction["display_name"], "text": reaction["text"],
                     "source": reaction["source"], "inputs": reaction["inputs"], "outputs": reaction["outputs"],
                     "catalysts": reaction["catalysts"], "regulators": reaction["regulators"],
                     "manual_status": "pending"})
    out = root / "review"
    out.mkdir(parents=True, exist_ok=True)
    write_jsonl(out / "reactions.jsonl", rows)
    write_jsonl(out / "pointer_checks.jsonl", checks)
    edges = list(read_jsonl(root / "graph" / "edges.jsonl"))
    write_jsonl(out / "edges.jsonl", edges)
    boundary = []
    if context:
        for edge in read_jsonl(Path(context) / "graph" / "edges.jsonl"):
            if (edge["source"] in reactions) != (edge["target"] in reactions):
                boundary.append(edge)
    write_jsonl(out / "context_boundary_edges.jsonl", boundary)
    # User edits are kept separately and are never overwritten when regenerating.
    review_path = out / "manual_review.csv"
    if not review_path.exists():
        with review_path.open("w", encoding="utf-8", newline="") as stream:
            fields = ["item_type", "item_id", "status", "reviewer", "notes"]
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for kind, ids in (("reaction", sorted(reactions)), ("edge", [x["id"] for x in edges])):
                writer.writerows({"item_type": kind, "item_id": rid, "status": "pending"} for rid in ids)
    with review_path.open(encoding="utf-8", newline="") as stream:
        reviews = list(csv.DictReader(stream))
    expected_items = {("reaction", rid) for rid in reactions} | {("edge", edge["id"]) for edge in edges}
    statuses = {key: "pending" for key in expected_items}
    for row in reviews:
        key = (row["item_type"], row["item_id"])
        if key in statuses:
            statuses[key] = row["status"]
    counts = {status: sum(v == status for v in statuses.values()) for status in sorted(set(statuses.values()))}
    summary = {"pathway_id": manifest.get("pathway_id"), "pathway_name": manifest.get("pathway_name"),
               "expected_reactions": manifest["expected_reactions"], "present_reactions": len(reactions),
               "membership_complete": len(reactions) == manifest["expected_reactions"],
               "source_pointer_checks": len(checks), "source_pointer_passes": sum(c["status"] == "pass" for c in checks),
               "internal_evidence_edges": len(edges), "context_boundary_evidence_edges": len(boundary),
               "context_graph_sha256": read_json(Path(context) / "graph" / "status.json")["graph_sha256"] if context else None,
               "manual_review_status_counts": counts}
    write_json(out / "summary.json", summary)
    lines = ["# Pathway review packet", "", f"{manifest.get('pathway_id')}: {manifest.get('pathway_name')}", "",
             "Internal chains use the complete pathway case. Context boundary edges use the supplied cached graph.", "",
             "Review inputs, outputs, stoichiometry, catalysts and regulation against each source JSON. Then review each edge's shared stable ID and role. Record decisions in manual_review.csv.", ""]
    for row in rows:
        lines.extend(["## " + row["reaction_id"], "", row["display_name"], "", row["text"], "",
                      "Source: `" + row["source"]["file"] + "`; SHA-256: `" + row["source"]["sha256"] + "`.", ""])
        for role in ("inputs", "outputs", "catalysts", "regulators"):
            lines += ["- " + role + ": " + "; ".join(str(x["entity_id"]) + " @ `" + x["source_pointer"] + "`" for x in row[role])]
        lines.append("")
    lines.extend(["## Edge evidence", ""])
    for edge in edges:
        lines.append(f"- {edge['source']} → {edge['target']}: {edge['kind']}; {edge['entity_id']} ({edge['entity_name']}); source {edge['source_pointers']}; target {edge['target_pointers']}")
    write_bytes(out / "review.md", ("\n".join(lines) + "\n").encode())
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "prepare-classes", "packet", "freeze"))
    parser.add_argument("--out", required=True)
    parser.add_argument("--source", help="Parent full-inventory snapshot, also used as cache")
    parser.add_argument("--pathway")
    parser.add_argument("--context", help="Existing graph for external-edge context")
    parser.add_argument("--export", help="Frozen example destination")
    args = parser.parse_args()
    with SnapshotLock(args.out):
        if args.stage == "prepare":
            if not args.source or not args.pathway:
                parser.error("prepare requires --source and --pathway")
            parent = load_snapshot(args.source)
            result = prepare(args.source, args.out, args.pathway, Client(parent["base_url"], timeout=45, retries=2, interval=0.5))
        elif args.stage == "prepare-classes":
            if not args.source:
                parser.error("prepare-classes requires --source")
            parent = load_snapshot(args.source)
            result = prepare_classes(args.source, args.out, Client(parent["base_url"], timeout=45, retries=2, interval=0.5))
        elif args.stage == "packet":
            result = packet(args.out, args.context)
        else:
            if not args.export:
                parser.error("freeze requires --export")
            result = freeze(args.out, args.export)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0 if result.get("complete", result.get("membership_complete")) else 2


if __name__ == "__main__":
    raise SystemExit(main())
