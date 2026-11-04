"""Resolve JSOG references and extract only the root reaction's roles."""

import gzip
import hashlib
import json
import sqlite3
from collections import Counter
from pathlib import Path

from .common import as_list, digest, identity, read_json, write_json, write_jsonl
from .fetch import load_snapshot, raw_path
from .render import reaction_text


class Resolver:
    def __init__(self, root):
        self.ids, self.db_ids = {}, {}
        self.locations = {}
        self.issues = []
        stack = [(root, "")]
        while stack:
            obj, pointer = stack.pop()
            if isinstance(obj, dict):
                self.locations[id(obj)] = pointer
                if "@id" in obj:
                    self.ids[str(obj["@id"])] = obj
                if "dbId" in obj:
                    key = str(obj["dbId"])
                    if len(obj) > len(self.db_ids.get(key, {})):
                        self.db_ids[key] = obj
                stack.extend((value, pointer + "/" + str(key).replace("~", "~0").replace("/", "~1")) for key, value in obj.items())
            elif isinstance(obj, list):
                stack.extend((value, pointer + "/" + str(index)) for index, value in reversed(list(enumerate(obj))))

    def resolve(self, value):
        if isinstance(value, dict):
            if "@ref" in value:
                key = str(value["@ref"])
                if key not in self.ids:
                    self.issues.append({"type": "unresolved_jsog_reference", "reference": key})
                return self.ids.get(key, {})
            # A stub may precede the richer definition elsewhere in this response.
            rich = self.db_ids.get(str(value.get("dbId")), {})
            return rich if len(rich) > len(value) else value
        if isinstance(value, (int, str)):
            result = self.db_ids.get(str(value))
            if result is None:
                self.issues.append({"type": "unresolved_db_reference", "reference": str(value)})
                return {"dbId": value} if str(value).isdigit() else {"stId": str(value)}
            return result
        return {}

    def field(self, obj, *names):
        for name in names:
            if name in obj:
                return name, as_list(obj[name])
        return names[0], []

    def objects(self, obj, *names):
        field, values = self.field(obj, *names)
        for index, value in enumerate(values):
            pointer = self.locations.get(id(obj), "") + "/" + field
            if isinstance(obj.get(field), list):
                pointer += "/" + str(index)
            wrapper = self.resolve(value)
            if "element" in wrapper:
                yield self.resolve(wrapper["element"]), pointer + "/element", wrapper
            else:
                yield wrapper, pointer, {}

    def entity(self, obj, pointer, wrapper=None):
        obj = self.resolve(obj)
        references = [item for item, _, _ in self.objects(obj, "referenceEntity")]
        chebi_ids = set()
        for reference in references:
            identifier = str(reference.get("identifier", ""))
            database = str(reference.get("databaseName", reference.get("referenceDatabase", ""))).lower()
            if "chebi" in database or str(reference.get("stId", "")).lower().startswith("chebi:"):
                if not identifier:
                    identifier = str(reference.get("stId", ""))
                if identifier:
                    chebi_ids.add("CHEBI:" + identifier.removeprefix("CHEBI:").removeprefix("chebi:"))
        compartments = sorted({x.get("displayName", x.get("name", "")) for x, _, _ in self.objects(obj, "hasCompartment", "compartment") if isinstance(x.get("displayName", x.get("name", "")), str)})
        row = {"entity_id": obj.get("stId"), "db_id": obj.get("dbId"),
               "display_name": obj.get("displayName"), "schema_class": obj.get("schemaClass"),
               "compartments": compartments, "chebi_ids": sorted(chebi_ids), "source_pointer": pointer}
        if wrapper and "stoichiometry" in wrapper:
            row["stoichiometry"] = wrapper["stoichiometry"]
        else:
            row["stoichiometry"] = None
        if not row["entity_id"]:
            self.issues.append({"type": "participant_without_stid", "pointer": pointer, "db_id": row["db_id"]})
        if not row["display_name"]:
            self.issues.append({"type": "participant_without_name", "pointer": pointer, "entity_id": row["entity_id"]})
        return row


def normalize_one(obj, source=None):
    resolver = Resolver(obj)
    row = {"id": identity(obj), "db_id": obj.get("dbId"), "st_id": obj.get("stId"),
           "st_id_version": obj.get("stIdVersion"), "display_name": obj.get("displayName"),
           "schema_class": obj.get("schemaClass"), "is_inferred": obj.get("isInferred"),
           "is_in_disease": obj.get("isInDisease"), "species_name": obj.get("speciesName"),
           "inputs": [], "outputs": [], "catalysts": [], "regulators": [],
           "preceding_events": [], "negative_preceding_events": []}
    for output, fields in (("inputs", ("inputs", "input")), ("outputs", ("outputs", "output"))):
        for entity, pointer, wrapper in resolver.objects(obj, *fields):
            row[output].append(resolver.entity(entity, pointer, wrapper))
    for activity, pointer, _ in resolver.objects(obj, "catalystActivity", "catalystActivities"):
        for entity, subpointer, _ in resolver.objects(activity, "physicalEntity"):
            item = resolver.entity(entity, subpointer)
            item["activity_db_id"] = activity.get("dbId")
            item["activity_source_pointer"] = pointer
            row["catalysts"].append(item)
    for regulation, pointer, _ in resolver.objects(obj, "regulatedBy", "regulations"):
        cls = regulation.get("schemaClass", regulation.get("className", "Regulation"))
        if cls in ("PositiveRegulation", "PositiveGeneExpressionRegulation", "Requirement"):
            kind, sign = "positive_regulation", 1
        elif cls in ("NegativeRegulation", "NegativeGeneExpressionRegulation"):
            kind, sign = "negative_regulation", -1
        else:
            kind, sign = "unknown_regulation", None
        for entity, subpointer, _ in resolver.objects(regulation, "regulator"):
            item = resolver.entity(entity, subpointer)
            item.update({"regulation_db_id": regulation.get("dbId"), "regulation_class": cls, "kind": kind, "sign": sign})
            item["regulation_source_pointer"] = pointer
            row["regulators"].append(item)
    row["compartments"] = sorted({x.get("displayName", "") for x, _, _ in resolver.objects(obj, "hasCompartment", "compartment") if x.get("displayName")})
    row["species"] = sorted([{"db_id": x.get("dbId"), "tax_id": x.get("taxId"), "name": x.get("displayName")} for x, _, _ in resolver.objects(obj, "species")], key=digest)
    for dest, fields in (("pathways", ("eventOfPathways", "eventOf")), ("inferred_from", ("inferredFrom",))):
        row[dest] = sorted([{"id": identity(x), "schema_class": x.get("schemaClass"), "display_name": x.get("displayName")} for x, _, _ in resolver.objects(obj, *fields)], key=digest)
    for event, pointer, _ in resolver.objects(obj, "precedingEvent", "precedingEvents"):
        row["preceding_events"].append({"id": identity(event), "schema_class": event.get("schemaClass"), "source_pointer": pointer})
    for rejected, pointer, _ in resolver.objects(obj, "negativePrecedingEvent"):
        reasons = [reason.get("displayName") for reason, _, _ in resolver.objects(rejected, "reason")]
        for event, subpointer, _ in resolver.objects(rejected, "precedingEvent"):
            row["negative_preceding_events"].append({"id": identity(event), "reasons": sorted(x for x in reasons if x), "source_pointer": subpointer})
    for field in ("inputs", "outputs", "catalysts", "regulators", "preceding_events", "negative_preceding_events"):
        row[field].sort(key=lambda item: (item.get("entity_id") or item.get("id") or "", item.get("source_pointer", "")))
    row["source"] = source or {}
    row["text"] = reaction_text(row)
    return row, resolver.issues


def normalize(root, allow_partial=False):
    root = Path(root)
    manifest = load_snapshot(root)
    if (root / "release_mismatch.json").exists():
        raise ValueError("Snapshot has a recorded release mismatch")
    receipts = {}
    if (root / "downloads.sqlite").exists():
        with sqlite3.connect(root / "downloads.sqlite") as db:
            receipts = dict(db.execute("SELECT db_id,sha256 FROM receipts"))
    missing, issues, rows = [], [], []
    for record in manifest["inventory"]:
        path = raw_path(root, record["dbId"])
        if not path.exists():
            missing.append(record["dbId"])
            continue
        body = gzip.decompress(path.read_bytes())
        sha = hashlib.sha256(body).hexdigest()
        if record["dbId"] in receipts and receipts[record["dbId"]] != sha:
            raise ValueError(f"Raw checksum differs from download receipt: {path}")
        obj = json.loads(body)
        if obj.get("dbId") != record["dbId"] or (record.get("stId") and obj.get("stId") != record["stId"]):
            raise ValueError(f"Raw identity mismatch: {path}")
        source = {"file": path.relative_to(root).as_posix(), "sha256": sha,
                  "url": manifest["base_url"] + "/data/query/enhanced/v2/" + str(record["dbId"]),
                  "query": manifest["query"], "reactome_release": manifest["reactome_release"]}
        row, warnings = normalize_one(obj, source)
        rows.append(row)
        issues.extend({"reaction_id": row["id"], **warning} for warning in warnings)
    if missing and not allow_partial:
        raise ValueError(f"{len(missing)} reactions are not downloaded; resume fetch or explicitly use --allow-partial")
    rows.sort(key=lambda item: item["id"])
    out = root / "normalized"
    write_jsonl(out / "reactions.jsonl", rows)
    write_jsonl(out / "issues.jsonl", sorted(issues, key=digest))
    status = {"expected": manifest["expected_reactions"], "normalized": len(rows), "complete": not missing,
              "missing_db_ids": sorted(missing), "issues": dict(sorted(Counter(x["type"] for x in issues).items())),
              "inventory_sha256": manifest["inventory_sha256"], "reaction_records_sha256": digest(rows)}
    write_json(out / "status.json", status)
    return status
