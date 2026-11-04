"""Create a local synthetic snapshot for inspecting every pipeline stage."""

import argparse
import gzip
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from reactome_pipeline.common import digest, dumps, read_json, write_bytes, write_json
from reactome_pipeline.fetch import QUERY, raw_path


def create_example(out):
    out = Path(out)
    rows = read_json(Path(__file__).parent / "input" / "reactions.json")
    inventory = sorted([{k: row[k] for k in ("dbId", "stId", "schemaClass", "displayName")} for row in rows], key=lambda row: (row["stId"], row["dbId"]))
    manifest = {"format_version": 1, "reactome_release": "synthetic-example", "base_url": "https://example.invalid/ContentService",
                "scope": "six synthetic reactions illustrating edge types, a cycle and currency exclusion",
                "query": QUERY, "expected_reactions": len(rows), "inventory": inventory, "inventory_sha256": digest(inventory)}
    if (out / "snapshot.json").exists() and read_json(out / "snapshot.json") != manifest:
        raise ValueError("Output directory contains another snapshot")
    write_json(out / "snapshot.json", manifest)
    for row in rows:
        write_bytes(raw_path(out, row["dbId"]), gzip.compress(dumps(row).encode(), mtime=0))
    return out


def create_real_example(out):
    source = Path(__file__).parent / "real_snapshot"
    out = Path(out)
    if (out / "snapshot.json").exists() and read_json(out / "snapshot.json") != read_json(source / "snapshot.json"):
        raise ValueError("Output directory contains another snapshot")
    shutil.copytree(source, out, dirs_exist_ok=True, ignore=shutil.ignore_patterns("*.sqlite", "normalized", "graph", "paths", "audit"))
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="runs/example")
    parser.add_argument("--real", action="store_true", help="Copy the frozen two-reaction Reactome example")
    args = parser.parse_args()
    (create_real_example if args.real else create_example)(args.out)
