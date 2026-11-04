"""English templates; names, compartments and stoichiometry come from source fields."""


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

    text = terms("inputs") + " → " + terms("outputs")
    text += "; catalyst: " + (", ".join(sorted({participant_text(x) for x in row["catalysts"]})) or "[not recorded]")
    for kind, label in (("positive_regulation", "positive regulators"), ("negative_regulation", "negative regulators"), ("unknown_regulation", "regulators (sign unspecified)")):
        items = sorted({participant_text(x) for x in row["regulators"] if x["kind"] == kind})
        if items:
            text += "; " + label + ": " + ", ".join(items)
    return text


def chain_text(reaction_ids, reactions, transitions):
    lines = []
    labels = {"material_flow": "provides an input to", "catalyst_dependency": "provides a catalyst for",
              "positive_regulation": "provides a positive regulator for", "negative_regulation": "provides a negative regulator for",
              "unknown_regulation": "provides a regulator with unspecified sign for"}
    for index, reaction_id in enumerate(reaction_ids):
        lines.append(f"{index + 1}. {reactions[reaction_id]['text']}")
        if index < len(transitions):
            for evidence in transitions[index]["evidence"]:
                lines.append("   Connection: " + labels[evidence["kind"]] + " the next reaction via " + evidence["entity_name"] + ".")
    return "\n".join(lines)
