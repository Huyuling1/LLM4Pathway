"""English templates; names, compartments and stoichiometry come from source fields."""

import json


def reaction_name(row):
    return json.dumps(row.get("display_name") or "[name not recorded]", ensure_ascii=False)


def participant_text(item, coefficient=False):
    text = item.get("display_name") or item.get("entity_id") or ("db:" + str(item.get("db_id")))
    for compartment in item.get("compartments", []):
        if "[" + compartment + "]" not in text:
            text += " [" + compartment + "]"
    count = item.get("stoichiometry")
    if coefficient and count is not None and count != 1:
        text = str(count) + " " + text
    return text


def reaction_text(row):
    def terms(key):
        return " + ".join(participant_text(item, True) for item in row[key]) or "[not recorded]"

    text = reaction_name(row) + ": " + terms("inputs") + " → " + terms("outputs")
    text += "; catalyst: " + (", ".join(sorted({participant_text(x) for x in row["catalysts"]})) or "[not recorded]")
    for kind, label in (("positive_regulation", "positive regulators"), ("negative_regulation", "negative regulators"), ("unknown_regulation", "regulators (sign unspecified)")):
        items = sorted({participant_text(x) for x in row["regulators"] if x["kind"] == kind})
        if items:
            text += "; " + label + ": " + ", ".join(items)
    return text


def chain_text(reaction_ids, reactions, transitions, closing_edges=None):
    lines = []
    labels = {"material_flow": "an input to", "catalyst_dependency": "a catalyst for",
              "positive_regulation": "a positive regulator of", "negative_regulation": "an inhibitor of",
              "unknown_regulation": "a regulator of"}
    for index, reaction_id in enumerate(reaction_ids):
        lines.append(reactions[reaction_id]["text"])
        if index < len(transitions):
            for evidence in transitions[index]["evidence"]:
                lines.append("   Connection: " + reaction_name(reactions[reaction_id])
                             + " produces " + evidence["entity_name"] + ", which is " + labels[evidence["kind"]]
                             + " " + reaction_name(reactions[transitions[index]["target"]])
                             + (" with unspecified sign" if evidence["kind"] == "unknown_regulation" else "") + ".")
    for edge in closing_edges or []:
        lines.append("closed-loop connection: " + reaction_name(reactions[edge["source"]])
                     + " → " + reaction_name(reactions[edge["target"]])
                     + ".")
    return "\n".join(lines)
