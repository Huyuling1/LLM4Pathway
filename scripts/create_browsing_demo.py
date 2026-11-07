"""Build the species/pathway export offline from the frozen seven-reaction case."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from examples.create_example import create_pathway_example
from reactome_pipeline.bulk import process
from reactome_pipeline.common import file_sha, read_json, write_json, write_jsonl


def create(dataset):
    dataset = Path(dataset)
    root = dataset / "global"
    create_pathway_example(root)
    snapshot = read_json(root / "snapshot.json")
    source = read_json(root / "raw/membership/pathway.json")
    pid = snapshot["pathway_id"]
    species = source.get("species", [])
    if isinstance(species, dict):
        species = [species]
    entry = {"id": pid, "display_name": source["displayName"],
             "species": [{"tax_id": str(s["taxId"]) if s.get("taxId") else None,
                          "name": s.get("displayName")} for s in species],
             "membership_sources": snapshot["membership_sources"]}
    write_json(dataset / "pathway_hierarchy.json", {
        "scope": "offline single-pathway browsing demonstration; not the full Reactome hierarchy",
        "reaction_inventory_sha256": snapshot["inventory_sha256"], "pathways": {pid: entry}})
    write_jsonl(dataset / "catalog/pathway_membership.jsonl", [
        {"pathway_id": pid, "reaction_ids": [r["stId"] for r in snapshot["inventory"]]}])
    write_json(dataset / "catalog/status.json", {"complete": True, "scope": "single frozen pathway",
        "hierarchy_sha256": file_sha(dataset / "pathway_hierarchy.json"),
        "membership_sha256": file_sha(dataset / "catalog/pathway_membership.jsonl")})
    return process(dataset, [pid], min_free_gb=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="runs/browsing-demo")
    args = parser.parse_args()
    print(create(args.dataset))
