"""Compare saved normalized roles with raw JSON slots, without calling the normalizer."""

import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from .common import digest, file_sha, read_jsonl, write_json, write_jsonl
from .fetch import load_snapshot, raw_path

AUDIT_VERSION = "raw-role-comparison-v1"
ROLE_FIELDS = {"inputs": ("inputs", "input"), "outputs": ("outputs", "output"),
               "catalysts": ("catalystActivity", "catalystActivities"),
               "regulators": ("regulatedBy", "regulations")}


class SourceDocument:
    """Index reference definitions, then visit only the root reaction's role slots."""

    def __init__(self, raw):
        self.locations, self.jsog, self.db = {}, {}, defaultdict(list)
        self.problems = []
        pending = [(raw, "")]
        while pending:
            value, pointer = pending.pop()
            if isinstance(value, dict):
                self.locations[id(value)] = pointer
                if "@id" in value:
                    self.jsog[str(value["@id"])] = value
                if "dbId" in value:
                    self.db[str(value["dbId"])].append(value)
                pending.extend((v, pointer + "/" + str(k).replace("~", "~0").replace("/", "~1")) for k, v in value.items())
            elif isinstance(value, list):
                pending.extend((v, pointer + "/" + str(i)) for i, v in enumerate(value))

    def resolve(self, value, pointer, seen=None):
        seen = set() if seen is None else seen
        if isinstance(value, dict) and "@ref" in value:
            key = str(value["@ref"])
            if key not in self.jsog or key in seen:
                self.problems.append({"type": "unresolved_source_reference", "pointer": pointer, "reference": key})
                return {}
            return self.resolve(self.jsog[key], pointer, seen | {key})
        if isinstance(value, dict):
            # Inline slots take precedence. Supplement only fields absent from this
            # occurrence; conflicting core identity fields are explicit audit issues.
            candidates = self.db.get(str(value.get("dbId")), [])
            result = dict(value)
            location = self.locations.get(id(value), pointer)
            for field in ("stId", "displayName", "schemaClass"):
                values = {x[field] for x in candidates if isinstance(x.get(field), str)}
                if len(values) > 1:
                    self.problems.append({"type": "conflicting_source_definitions", "pointer": pointer, "field": field,
                                          "values": sorted(values)})
            # Use a real source object's location for role fields, never a synthesized pointer.
            richer = sorted(candidates, key=lambda x: (-len(x), self.locations[id(x)]))
            for candidate in richer:
                for field, content in candidate.items():
                    result.setdefault(field, content)
            self.locations[id(result)] = location
            self.field_locations = getattr(self, "field_locations", {})
            for field in result:
                owner = value if field in value else next(x for x in richer if field in x)
                self.field_locations[(id(result), field)] = self.locations.get(id(owner), pointer)
            return result
        if isinstance(value, (int, str)) and str(value) in self.db:
            return self.resolve(self.db[str(value)][0], pointer, seen)
        self.problems.append({"type": "unresolved_source_reference", "pointer": pointer, "reference": value})
        return {}

    def slots(self, obj, aliases):
        present = [name for name in aliases if name in obj]
        if not present:
            return [], "absent", None
        field = present[0]
        for alternative in present[1:]:
            if obj[alternative] and obj[alternative] != obj[field]:
                self.problems.append({"type": "multiple_role_aliases", "fields": present})
        value = obj[field]
        if value is None:
            return [], "null", field
        if value == []:
            return [], "empty_list", field
        parent = getattr(self, "field_locations", {}).get((id(obj), field), self.locations.get(id(obj), ""))
        pairs = [(item, parent + "/" + field + "/" + str(i)) for i, item in enumerate(value)] if isinstance(value, list) else [(value, parent + "/" + field)]
        result = []
        for item, pointer in pairs:
            resolved = self.resolve(item, pointer)
            coefficient = None
            if "element" in resolved:
                coefficient = resolved.get("stoichiometry")
                pointer += "/element"
                resolved = self.resolve(resolved["element"], pointer)
            result.append((resolved, pointer, coefficient))
        return result, "populated", field

    def participant(self, obj, pointer, coefficient):
        compartments, _, _ = self.slots(obj, ("hasCompartment", "compartment"))
        names = [x.get("displayName", x.get("name")) for x, _, _ in compartments]
        refs, _, _ = self.slots(obj, ("referenceEntity",))
        chebi = set()
        for ref, _, _ in refs:
            database = str(ref.get("databaseName", ref.get("referenceDatabase", ""))).lower()
            if "chebi" in database or str(ref.get("stId", "")).lower().startswith("chebi:"):
                identifier = str(ref.get("identifier", "")) or str(ref.get("stId", ""))
                if identifier:
                    chebi.add("CHEBI:" + identifier.removeprefix("CHEBI:").removeprefix("chebi:"))
        return {"entity_id": obj.get("stId"), "db_id": obj.get("dbId"), "display_name": obj.get("displayName"),
                "schema_class": obj.get("schemaClass"), "stoichiometry": coefficient,
                "compartments": sorted({name for name in names if isinstance(name, str)}),
                "chebi_ids": sorted(chebi), "source_pointer": pointer}

    def role(self, raw, role):
        slots, state, field = self.slots(raw, ROLE_FIELDS[role])
        expected, empty_groups = [], []
        for group, pointer, coefficient in slots:
            if role in ("inputs", "outputs"):
                expected.append(self.participant(group, pointer, coefficient))
                continue
            actors, actor_state, _ = self.slots(group, ("physicalEntity",) if role == "catalysts" else ("regulator",))
            if not actors:
                empty_groups.append({"source_pointer": pointer, "actor_field_state": actor_state})
            for actor, actor_pointer, _ in actors:
                item = self.participant(actor, actor_pointer, None)
                if role == "catalysts":
                    item["activity_db_id"] = group.get("dbId")
                else:
                    cls = group.get("schemaClass", group.get("className", "Regulation"))
                    kind, sign = ("positive_regulation", 1) if cls in {"PositiveRegulation", "PositiveGeneExpressionRegulation", "Requirement"} else ("negative_regulation", -1) if cls in {"NegativeRegulation", "NegativeGeneExpressionRegulation"} else ("unknown_regulation", None)
                    item.update({"regulation_db_id": group.get("dbId"), "regulation_class": cls, "kind": kind, "sign": sign})
                expected.append(item)
        return expected, {"source_field": field, "source_field_state": state, "source_role_slots": len(slots),
                          "groups_without_participants": empty_groups}


def compare_reaction(raw, normalized):
    source = SourceDocument(raw)
    issues, checks, empty_roles = [], [], []
    for field, key in (("dbId", "db_id"), ("stId", "st_id"), ("stIdVersion", "st_id_version"),
                       ("displayName", "display_name"), ("schemaClass", "schema_class"), ("isInferred", "is_inferred")):
        if raw.get(field) != normalized.get(key):
            issues.append({"type": "reaction_field_mismatch", "field": key, "source": raw.get(field), "normalized": normalized.get(key)})
    expected_id = raw.get("stId") or ("db:" + str(raw["dbId"]) if raw.get("dbId") is not None else None)
    if normalized.get("id") != expected_id:
        issues.append({"type": "reaction_field_mismatch", "field": "id", "source": expected_id, "normalized": normalized.get("id")})
    for role in ROLE_FIELDS:
        before = len(source.problems)
        expected, description = source.role(raw, role)
        actual = normalized.get(role, [])
        expected_slots = {item["source_pointer"]: item for item in expected}
        actual_slots = defaultdict(list)
        for item in actual:
            actual_slots[item.get("source_pointer")].append(item)
        role_issues = []
        for pointer, wanted in sorted(expected_slots.items()):
            found = actual_slots.get(pointer, [])
            if not found:
                role_issues.append({"type": "missing_participant", "source_pointer": pointer, "expected": wanted})
            elif len(found) != 1:
                role_issues.append({"type": "duplicate_participant", "source_pointer": pointer, "count": len(found)})
            else:
                for field, value in wanted.items():
                    if found[0].get(field) != value:
                        role_issues.append({"type": "participant_field_mismatch", "source_pointer": pointer,
                                            "field": field, "source": value, "normalized": found[0].get(field)})
        for pointer in actual_slots.keys() - expected_slots.keys():
            role_issues.append({"type": "unexpected_participant", "source_pointer": pointer})
        issues.extend({"role": role, **problem} for problem in role_issues)
        unresolved = any(p["type"] == "unresolved_source_reference" for p in source.problems[before:])
        check = {"role": role, **description, "source_participants": len(expected), "normalized_participants": len(actual),
                 "comparison_status": "mismatch" if role_issues else "source_unresolved" if unresolved else "match"}
        checks.append(check)
        if not actual:
            reason = "unresolved_source_reference" if unresolved else "extraction_omission" if expected else (
                "source_role_present_without_participants" if description["source_role_slots"] else "source_field_" + description["source_field_state"])
            empty_roles.append({**check, "reason": reason})
    issues.extend(source.problems)
    return checks, empty_roles, issues


def audit_extraction(root):
    root = Path(root)
    manifest = load_snapshot(root)
    records, duplicates = {}, []
    for record in read_jsonl(root / "normalized/reactions.jsonl"):
        dbid = record.get("db_id")
        if dbid in records:
            duplicates.append({"type": "duplicate_normalized_reaction", "db_id": dbid})
        records[dbid] = record
    expected_classes = Counter(r.get("schemaClass", "[not recorded]") for r in manifest["inventory"])
    available_classes, checked_classes = Counter(), Counter()
    all_checks, all_empty, all_issues = [], [], list(duplicates)
    unavailable, not_normalized, hashes = [], [], []
    for entry in manifest["inventory"]:
        dbid = entry["dbId"]
        path = raw_path(root, dbid)
        if not path.exists():
            unavailable.append(dbid)
            if dbid in records:
                all_issues.append({"type": "normalized_without_raw", "db_id": dbid})
            continue
        body = gzip.decompress(path.read_bytes())
        raw = json.loads(body)
        sha = hashlib.sha256(body).hexdigest()
        hashes.append({"db_id": dbid, "sha256": sha})
        cls = raw.get("schemaClass", "[not recorded]")
        available_classes[cls] += 1
        context = {"reaction_id": entry.get("stId"), "db_id": dbid, "schema_class": cls,
                   "display_name": raw.get("displayName"), "source_file": path.relative_to(root).as_posix(), "source_sha256": sha}
        if raw.get("dbId") != dbid or raw.get("stId") != entry.get("stId"):
            all_issues.append({**context, "type": "raw_identity_mismatch"})
        if dbid not in records:
            not_normalized.append(dbid)
            all_issues.append({**context, "type": "missing_normalized_reaction"})
            continue
        row = records[dbid]
        if row.get("source", {}).get("file") != path.relative_to(root).as_posix():
            all_issues.append({**context, "type": "source_file_mismatch"})
        if row.get("source", {}).get("sha256") != sha:
            all_issues.append({**context, "type": "source_checksum_mismatch"})
        checks, empty, issues = compare_reaction(raw, row)
        all_checks.extend({**context, **check} for check in checks)
        all_empty.extend({**context, **item} for item in empty)
        all_issues.extend({**context, **item} for item in issues)
        checked_classes[cls] += 1
    known = {row["dbId"] for row in manifest["inventory"]}
    for dbid in records.keys() - known:
        all_issues.append({"type": "normalized_outside_inventory", "db_id": dbid})
    class_rows = [{"schema_class": cls, "inventory": expected_classes[cls], "raw_available": available_classes[cls],
                   "compared": checked_classes[cls]} for cls in sorted(expected_classes.keys() | available_classes.keys())]
    summary = {"audit_version": AUDIT_VERSION, "scope": manifest["scope"], "reactome_release": manifest["reactome_release"],
               "inventory_sha256": manifest["inventory_sha256"], "normalized_file_sha256": file_sha(root / "normalized/reactions.jsonl"),
               "raw_receipts_sha256": digest(hashes), "inventory_reactions": manifest["expected_reactions"],
               "raw_available": sum(available_classes.values()), "reactions_compared": sum(checked_classes.values()),
               "input_complete": not unavailable, "available_records_compared": not not_normalized,
               "source_comparison_passed": bool(sum(checked_classes.values())) and not all_issues and not not_normalized,
               "role_checks": len(all_checks), "matched_role_checks": sum(c["comparison_status"] == "match" for c in all_checks),
               "source_participants": sum(c["source_participants"] for c in all_checks),
               "normalized_participants": sum(c["normalized_participants"] for c in all_checks),
               "issue_counts": dict(sorted(Counter(x["type"] for x in all_issues).items())),
               "empty_roles": {role: dict(sorted(Counter(x["reason"] for x in all_empty if x["role"] == role).items())) for role in ROLE_FIELDS},
               "schema_class_coverage": class_rows,
               "classes_without_raw_examples": sorted(cls for cls in expected_classes if not available_classes[cls]),
               "observed_features": {
                   "regulation_classes": dict(sorted(Counter(p.get("regulation_class", "[not recorded]") for r in records.values() for p in r.get("regulators", [])).items())),
                   "participant_classes": dict(sorted(Counter(p.get("schema_class") or "[not recorded]" for r in records.values() for role in ROLE_FIELDS for p in r.get(role, [])).items())),
                   "reactions_with_multiple_catalysts": sum(len(r.get("catalysts", [])) > 1 for r in records.values()),
                   "reactions_with_multiple_regulators": sum(len(r.get("regulators", [])) > 1 for r in records.values()),
                   "inference": dict(sorted(Counter(str(r.get("is_inferred")) for r in records.values()).items()))},
               "missing_normalized_db_ids": not_normalized,
               "checks": ["role occurrence counts and source pointers", "entity stable IDs, dbIds, names and classes",
                          "stoichiometry, compartments and ChEBI references", "catalyst activity and regulation IDs, classes and signs"],
               "interpretation": "Structural source comparison only; absent/null/empty source fields are reported without assigning a biological cause."}
    out = root / "audit"
    write_jsonl(out / "extraction_checks.jsonl", all_checks)
    write_jsonl(out / "extraction_empty_roles.jsonl", all_empty)
    write_jsonl(out / "extraction_issues.jsonl", sorted(all_issues, key=digest))
    write_json(out / "extraction_report.json", summary)
    return summary
