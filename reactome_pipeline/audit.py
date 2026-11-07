"""Coverage, topology, currency exclusions and curated-order comparisons."""

from collections import Counter
from pathlib import Path

from .common import file_sha, read_json, read_jsonl, write_json, write_jsonl
from .graph import components, load_graph, weak_components
from .audit_details import chain_metrics, currency_metrics


def audit(root):
    root = Path(root)
    db, reactions, adjacency, reverse, metadata = load_graph(root)
    out = root / "audit"
    out.mkdir(parents=True, exist_ok=True)
    groups = components(adjacency, reverse)
    independent = weak_components(adjacency, reverse)
    write_jsonl(out / "weak_components.jsonl", ({"size": len(group), "reaction_ids": group} for group in independent))
    cyclic = [group for group in groups if len(group) > 1 or group[0] in adjacency[group[0]]]
    write_jsonl(out / "cyclic_components.jsonl", ({"size": len(group), "reaction_ids": group} for group in cyclic))
    pairs = {(source, target) for source in adjacency for target in adjacency[source]}
    preceding, rejected, outside = set(), [], []
    for rid, reaction in reactions.items():
        for event in reaction["preceding_events"]:
            if event["id"] in reactions:
                preceding.add((event["id"], rid))
            else:
                outside.append({"target": rid, **event})
        for event in reaction["negative_preceding_events"]:
            if (event["id"], rid) in pairs:
                rejected.append({"source": event["id"], "target": rid, "reasons": event["reasons"]})
    write_jsonl(out / "preceding_only.jsonl", ({"source": a, "target": b} for a, b in sorted(preceding - pairs)))
    write_jsonl(out / "matched_rejected_connections.jsonl", sorted(rejected, key=lambda x: (x["source"] or "", x["target"])))
    write_jsonl(out / "preceding_outside_inventory.jsonl", sorted(outside, key=lambda x: (x["target"], x.get("id") or "")))
    sources = sorted(node for node in adjacency if not reverse[node])
    sinks = sorted(node for node in adjacency if not adjacency[node])
    isolated = sorted(node for node in adjacency if not adjacency[node] and not reverse[node])
    def species_set(reaction):
        values = {str(item["tax_id"]) for item in reaction["species"] if item.get("tax_id")}
        return values or ({reaction["species_name"]} if reaction.get("species_name") else set())
    cross_species = sum(bool(species_set(reactions[a]) and species_set(reactions[b]) and not (species_set(reactions[a]) & species_set(reactions[b]))) for a, b in pairs)
    cross_pathway = 0
    for a, b in pairs:
        pa = {x["id"] for x in reactions[a]["pathways"] if x["id"]}
        pb = {x["id"] for x in reactions[b]["pathways"] if x["id"]}
        cross_pathway += bool(pa and pb and not pa & pb)
    report = {"scope": read_json(root / "snapshot.json").get("scope"),
              "graph": metadata, "normalization": read_json(root / "normalized" / "status.json"),
              "species": dict(sorted(Counter(r.get("species_name") or "[not recorded]" for r in reactions.values()).items())),
              "schema_classes": dict(sorted(Counter(r.get("schema_class") or "[not recorded]" for r in reactions.values()).items())),
              "inference": dict(sorted(Counter(str(r["is_inferred"]) for r in reactions.values()).items())),
              "edge_types": dict(db.execute("SELECT kind,count(*) FROM edges GROUP BY kind ORDER BY kind")),
              "topology": {"source_nodes": len(sources), "sink_nodes": len(sinks), "isolated_nodes": len(isolated),
                           "branching_nodes": sum(len(x) > 1 for x in adjacency.values()),
                           "merging_nodes": sum(len(x) > 1 for x in reverse.values()),
                           "max_out_degree": max(map(len, adjacency.values()), default=0),
                           "max_in_degree": max(map(len, reverse.values()), default=0),
                           "self_loop_pairs": sum(a == b for a, b in pairs),
                           "strongly_connected_components": len(groups), "cyclic_components": len(cyclic),
                           "weakly_connected_components": len(independent),
                           "weak_component_size_histogram": dict(sorted(Counter(map(len, independent)).items())),
                           "reactions_in_cycles": sum(map(len, cyclic)),
                           "largest_cyclic_component": max(map(len, cyclic), default=0)},
              "preceding_event": {"pairs_within_graph": len(preceding), "overlap": len(preceding & pairs),
                                  "only_preceding": len(preceding - pairs), "only_entity_matching": len(pairs - preceding),
                                  "matched_negative_preceding": len(rejected)},
              "cross_species_pairs_with_disjoint_recorded_species": cross_species,
              "cross_pathway_pairs_with_disjoint_recorded_direct_memberships": cross_pathway,
              "currency": currency_metrics(root, db, reactions, adjacency, reverse)}
    if (root / "paths" / "status.json").exists():
        path_status = read_json(root / "paths" / "status.json")
        report["paths"] = path_status
        report["paths_match_current_graph"] = path_status["graph_sha256"] == metadata["graph_sha256"]
        if report["paths_match_current_graph"]:
            report["chain_metrics"] = chain_metrics(root, reactions, path_status)
    extraction_path = out / "extraction_report.json"
    if extraction_path.exists():
        extraction = read_json(extraction_path)
        report["extraction_audit"] = extraction
        report["extraction_matches_current_normalized"] = extraction["normalized_file_sha256"] == file_sha(root / "normalized/reactions.jsonl")
    write_json(out / "report.json", report)
    db.close()
    return report
