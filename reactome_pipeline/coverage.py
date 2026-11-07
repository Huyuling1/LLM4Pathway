"""Exact input matches to all earlier outputs in a chain; no inventory simulation."""

from collections import defaultdict

COVERAGE_VERSION = "prefix-input-stid-v1"


def input_coverage(reaction_ids, reactions):
    producers = defaultdict(list)
    steps = []
    for step, rid in enumerate(reaction_ids, 1):
        reaction = reactions[rid]
        covered, uncovered = [], []
        for item in reaction.get("inputs", []):
            entity = item.get("entity_id")
            matches = producers.get(entity, []) if entity else []
            entry = dict(item)
            entry["providers"] = [dict(match) for match in matches]
            if matches:
                covered.append(entry)
            else:
                entry["reason"] = "no_earlier_output_match" if entity else "missing_stable_identifier"
                uncovered.append(entry)
        steps.append({"reaction_id": rid, "step": step,
                      "recorded_inputs": len(reaction.get("inputs", [])),
                      "covered_inputs": covered, "uncovered_inputs": uncovered})
        # Insert only after matching: the current reaction cannot supply its own inputs.
        for item in reaction.get("outputs", []):
            if item.get("entity_id"):
                producers[item["entity_id"]].append({
                    "reaction_id": rid, "step": step,
                    "output_source_pointer": item.get("source_pointer"),
                    "stoichiometry": item.get("stoichiometry")})
    return steps
