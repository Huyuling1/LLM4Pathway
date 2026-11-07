"""Streaming chain summaries and frozen-currency exclusion attribution."""

import json
import math
import sqlite3
from collections import Counter, defaultdict

from .common import read_json, write_jsonl


def distribution(histogram):
    total = sum(histogram.values())
    def quantile(q):
        rank = max(1, math.ceil(q * total))
        count = 0
        for value, frequency in sorted(histogram.items()):
            count += frequency
            if count >= rank:
                return value
        return None
    return {"count": total, "min": min(histogram, default=None), "max": max(histogram, default=None),
            "mean": sum(k * v for k, v in histogram.items()) / total if total else None,
            "median_nearest_rank": quantile(0.5), "p95_nearest_rank": quantile(0.95),
            "histogram": {str(k): v for k, v in sorted(histogram.items())}}


def chain_metrics(root, reactions, path_status):
    counts, lengths, occurrences = Counter(), defaultdict(Counter), Counter()
    covered = uncovered = missing_id = unannotated = 0
    with sqlite3.connect((root / "paths" / "paths.sqlite").resolve().as_uri() + "?mode=ro", uri=True) as db:
        for (data,) in db.execute("SELECT data FROM chains ORDER BY ordinal"):
            chain = json.loads(data)
            counts[chain["status"]] += 1
            lengths[chain["status"]][len(chain["reaction_ids"])] += 1
            lengths["all"][len(chain["reaction_ids"])] += 1
            occurrences.update(chain["reaction_ids"])
            if "input_coverage" not in chain:
                unannotated += 1
            for step in chain.get("input_coverage", []):
                covered += len(step["covered_inputs"])
                uncovered += len(step["uncovered_inputs"])
                missing_id += sum(x["reason"] == "missing_stable_identifier" for x in step["uncovered_inputs"])
    total = sum(counts.values())
    appearances = sum(occurrences.values())
    covered_reactions = len(occurrences)
    write_jsonl(root / "audit" / "reaction_chain_coverage.jsonl",
                ({"reaction_id": rid, "display_name": reactions[rid].get("display_name"),
                  "chain_occurrences": occurrences[rid]} for rid in sorted(reactions)))
    return {"population": "emitted chains in this graph; DFS prefix if enumeration is unfinished",
            "enumeration_complete": path_status["enumeration_complete"], "input_complete": path_status["input_complete"],
            "total_chains": total, "sink_chains": counts["sink"], "cycle_stopped_chains": counts["cycle_blocked"],
            "sink_fraction": counts["sink"] / total if total else None,
            "cycle_stopped_fraction": counts["cycle_blocked"] / total if total else None,
            "lengths": {key: distribution(value) for key, value in sorted(lengths.items())},
            "reaction_coverage": {"graph_reactions": len(reactions), "covered_reactions": covered_reactions,
                "uncovered_reactions": len(reactions) - covered_reactions, "total_appearances": appearances,
                "mean_appearances_per_covered_reaction": appearances / covered_reactions if covered_reactions else None,
                "reactions_in_multiple_chains": sum(n > 1 for n in occurrences.values()),
                "repeat_appearance_fraction": (appearances - covered_reactions) / appearances if appearances else None,
                "occurrence_distribution": distribution(Counter(occurrences[rid] for rid in reactions)),
                "most_repeated": [{"reaction_id": rid, "count": count} for rid, count in sorted(occurrences.items(), key=lambda x: (-x[1], x[0]))[:10]]},
            "input_coverage": {"covered_input_occurrences": covered, "uncovered_input_occurrences": uncovered,
                "missing_id_occurrences": missing_id, "chains_without_annotations": unannotated,
                "covered_fraction": covered / (covered + uncovered) if covered + uncovered else None}}


def currency_metrics(root, db, reactions, adjacency, reverse):
    config = read_json(root / "graph" / "currency.frozen.json")
    groups = {x["chebi_id"]: {"name": x["name"], "entities": set()} for x in config["molecules"]}
    for eid in config.get("reactome_st_ids", []):
        groups["stId:" + eid] = {"name": eid, "entities": set()}
    sources_by_entity, targets_by_entity, data = {}, {}, {}
    before_in, before_out = set(), set()
    for eid, name, chebi in db.execute("SELECT id,name,chebi_ids FROM entities WHERE excluded=1 ORDER BY id"):
        sources = {r[0] for r in db.execute("SELECT DISTINCT reaction FROM participants WHERE entity=? AND role='outputs'", (eid,))}
        targets = {(r, k) for r, k in db.execute("SELECT DISTINCT reaction,kind FROM participants WHERE entity=? AND role!='outputs'", (eid,))}
        sources_by_entity[eid], targets_by_entity[eid] = sources, targets
        if sources and targets:
            before_out.update(sources)
            before_in.update(r for r, _ in targets)
        data[eid] = {"entity_id": eid, "display_name": name, "chebi_ids": json.loads(chebi),
                     "excluded_evidence_edges": len(sources) * len(targets)}
        for identifier in json.loads(chebi) + ["stId:" + eid]:
            if identifier in groups:
                groups[identifier]["entities"].add(eid)
    lost_in = {r for r in before_in if not reverse[r]}
    lost_out = {r for r in before_out if not adjacency[r]}
    isolated = {r for r in before_in | before_out if not reverse[r] and not adjacency[r]}
    def impact(entities):
        sources, targets, examples = set(), set(), []
        for eid in sorted(entities):
            producers, consumers = sources_by_entity[eid], targets_by_entity[eid]
            if producers and consumers:
                sources.update(producers)
                targets.update(r for r, _ in consumers)
            # Stable first examples, independent of set iteration and network ordering.
            for source in sorted(producers)[:3]:
                for target, kind in sorted(consumers)[:3]:
                    if len(examples) < 3:
                        examples.append({"source": source, "source_name": reactions[source].get("display_name"),
                            "target": target, "target_name": reactions[target].get("display_name"), "kind": kind,
                            "entity_id": eid, "entity_name": data[eid]["display_name"],
                            "pair_retained_by_other_evidence": target in adjacency[source]})
        return {"entity_ids": sorted(entities),
                "excluded_evidence_edges": sum(data[eid]["excluded_evidence_edges"] for eid in entities),
                "affected_reactions_losing_all_upstream": sorted(targets & lost_in),
                "affected_reactions_losing_all_downstream": sorted(sources & lost_out),
                "affected_reactions_becoming_isolated": sorted((sources | targets) & isolated),
                "examples": examples}
    per_entity = [{**data[eid], **impact({eid})} for eid in sorted(data)]
    per_molecule = [{"identifier": identifier, "name": group["name"], **impact(group["entities"])}
                    for identifier, group in sorted(groups.items())]
    write_jsonl(root / "audit" / "currency_entities.jsonl", per_entity)
    write_jsonl(root / "audit" / "currency_molecules.jsonl", per_molecule)
    return {"attribution": "Endpoint losses follow exclusion of the ENTIRE frozen list; per-molecule affected sets overlap. Evidence counts include self-loops and distinguish target roles.",
            "excluded_entities": len(data), "excluded_evidence_edges": sum(x["excluded_evidence_edges"] for x in data.values()),
            "reactions_losing_all_upstream": sorted(lost_in), "reactions_losing_all_downstream": sorted(lost_out),
            "reactions_becoming_isolated": sorted(isolated), "molecules": per_molecule}
