"""Freeze the complete pathway inventory and source-backed membership hierarchy."""

import argparse
import math
from collections import defaultdict, deque
from pathlib import Path

from .common import digest, file_sha, read_json, write_bytes, write_json, write_jsonl
from .fetch import Client, load_snapshot
from .locking import SnapshotLock
from .normalize import Resolver


def contained_reactions(pathways, children, reaction_ids):
    """Leaf-to-root closure; shared subpathways retain every parent."""
    members = {pid: set(children.get(pid, set())) & reaction_ids for pid in pathways}
    parents = defaultdict(set)
    remaining = {}
    for pid in pathways:
        subpaths = set(children.get(pid, set())) & pathways.keys()
        remaining[pid] = len(subpaths)
        for child in subpaths:
            parents[child].add(pid)
    queue = deque(sorted(pid for pid, n in remaining.items() if n == 0))
    visited = 0
    while queue:
        child = queue.popleft()
        visited += 1
        for parent in sorted(parents[child]):
            members[parent].update(members[child])
            remaining[parent] -= 1
            if remaining[parent] == 0:
                queue.append(parent)
    if visited != len(pathways):
        raise ValueError("Pathway containment contains a cycle; inspect frozen hierarchy sources")
    return members


def freeze_catalog(dataset, client):
    dataset = Path(dataset)
    snapshot = load_snapshot(dataset / "global")
    release = snapshot["reactome_release"]
    if client.base != snapshot["base_url"] or client.version() != release:
        raise ValueError("Catalog endpoint/release differs from frozen reactions")
    out = dataset / "catalog"
    out.mkdir(parents=True, exist_ok=True)
    receipts_path = out / "receipts.json"
    receipt_info = read_json(receipts_path) if receipts_path.exists() else {
        "reactome_release": release, "base_url": client.base, "sources": {}}
    if receipt_info["reactome_release"] != release or receipt_info["base_url"] != client.base:
        raise ValueError("Existing catalog belongs to another source")
    receipts = receipt_info["sources"]

    def get(name, endpoint, params=None):
        path = out / "raw" / (name + ".json")
        url = client.url(endpoint, params)
        if name in receipts:
            receipt = receipts[name]
            if receipt["url"] != url or file_sha(path) != receipt["sha256"]:
                raise ValueError("Catalog source checksum/URL mismatch: " + name)
            return read_json(path)
        body, obj = client.get(endpoint, params)
        if client.version() != release:
            write_json(out / "release_mismatch.json", {"expected": release, "observed_during": name})
            raise ValueError("Reactome release changed while downloading catalog")
        write_bytes(path, body)
        receipts[name] = {"file": path.relative_to(dataset).as_posix(), "url": url, "sha256": file_sha(path)}
        write_json(receipts_path, receipt_info)
        print("catalog cached " + name, flush=True)
        return obj

    if (out / "release_mismatch.json").exists():
        raise ValueError("Catalog has a recorded release mismatch")
    write_json(out / "status.json", {"complete": False, "stage": "fetching", "reactome_release": release})
    count = int(get("pathway-count", "/data/schema/Pathway/count"))
    inventory = {}
    for size in (20000, 19969, 19873):
        for page in range(1, math.ceil(count / size) + 1):
            rows = get(f"inventory-{size}-{page}", "/data/schema/Pathway/min", {"offset": size, "page": page})
            for row in rows:
                previous = inventory.setdefault(row["stId"], row)
                if previous != row:
                    raise ValueError("Pathway inventory changed: " + row["stId"])
        if len(inventory) == count:
            break
    if len(inventory) != count:
        raise ValueError(f"Pathway inventory incomplete: {len(inventory)}/{count}")
    write_json(out / "pathway_inventory.json", [inventory[k] for k in sorted(inventory)])
    species = get("main-species", "/data/species/main")
    pathways = {pid: {"id": pid, "db_id": row["dbId"], "display_name": row.get("displayName"),
                      "species": [], "membership_sources": []} for pid, row in sorted(inventory.items())}
    children, roots, seen = defaultdict(set), defaultdict(set), set()
    reactions = {r["stId"] for r in snapshot["inventory"] if r.get("stId")}
    event_by_db = {r["dbId"]: r.get("stId") for r in snapshot["inventory"]}
    event_by_db.update({r["dbId"]: pid for pid, r in inventory.items()})
    unknown = set()
    for organism in sorted(species, key=lambda s: str(s["taxId"])):
        taxid = str(organism["taxId"])
        name = "hierarchy-" + taxid
        hierarchy = get(name, "/data/eventsHierarchy/" + taxid)
        stack = [(node, None) for node in hierarchy]
        while stack:
            node, parent = stack.pop()
            eid = node["stId"]
            if parent:
                children[parent].add(eid)
            if eid in pathways:
                seen.add(eid)
                entry = pathways[eid]
                sp = {"tax_id": taxid, "name": organism["displayName"], "source": "browser_species_context"}
                if sp not in entry["species"]:
                    entry["species"].append(sp)
                if name not in entry["membership_sources"]:
                    entry["membership_sources"].append(name)
                if parent is None:
                    roots[taxid].add(eid)
                stack.extend((child, eid) for child in node.get("children", []))
            elif eid not in reactions:
                unknown.add(eid)
    # Main-species browser trees are not assumed to cover the complete inventory.
    for pid in sorted(pathways.keys() - seen):
        name = "pathway-" + pid
        obj = get(name, "/data/query/" + pid)
        if obj.get("stId") != pid or obj.get("schemaClass") not in ("Pathway", "TopLevelPathway"):
            raise ValueError("Pathway query identity mismatch: " + pid)
        resolver = Resolver(obj)
        for event, _, _ in resolver.objects(obj, "hasEvent"):
            eid = event.get("stId")
            if not eid:
                dbid = event.get("dbId")
                eid = event_by_db.get(dbid)
            if eid not in pathways and eid not in reactions:
                unknown.add(eid or "db:" + str(event.get("dbId")))
            else:
                children[pid].add(eid)
        pathways[pid]["membership_sources"].append(name)
        pathways[pid]["species"] = [{"tax_id": str(s["taxId"]) if s.get("taxId") else None,
                                       "name": s.get("displayName"), "source": "pathway_species"}
                                      for s, _, _ in resolver.objects(obj, "species")]
    if unknown:
        write_json(out / "unresolved_events.json", sorted(unknown))
        raise ValueError("Hierarchy contains events absent from frozen inventories; see unresolved_events.json")
    members = contained_reactions(pathways, children, reactions)
    for pid, entry in pathways.items():
        entry["children"] = sorted(children[pid])
        entry["species"].sort(key=digest)
        entry["membership_sources"].sort()
    hierarchy = {"format_version": 1, "reactome_release": release,
                 "reaction_inventory_sha256": snapshot["inventory_sha256"],
                 "pathway_inventory_sha256": digest([inventory[k] for k in sorted(inventory)]),
                 "roots_by_taxon": {k: sorted(v) for k, v in sorted(roots.items())},
                 "pathways": pathways, "sources": receipts}
    if client.version() != release:
        raise ValueError("Reactome release changed before catalog completion")
    write_json(dataset / "pathway_hierarchy.json", hierarchy)
    write_jsonl(out / "pathway_membership.jsonl", ({"pathway_id": pid, "reaction_ids": sorted(members[pid])}
                                                   for pid in sorted(pathways)))
    assigned = set().union(*members.values()) if members else set()
    write_json(out / "unassigned_reaction_ids.json", sorted(reactions - assigned))
    status = {"complete": True, "reactome_release": release, "pathways": count,
              "browser_species": len(species), "pathways_queried_separately": count - len(seen),
              "unassigned_reactions": len(reactions - assigned),
              "hierarchy_sha256": file_sha(dataset / "pathway_hierarchy.json"),
              "membership_sha256": file_sha(out / "pathway_membership.jsonl")}
    write_json(out / "status.json", status)
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    args = parser.parse_args()
    root = Path(args.dataset)
    with SnapshotLock(root / "catalog"):
        print(freeze_catalog(root, Client(load_snapshot(root / "global")["base_url"])))
