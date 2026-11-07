import copy
import tempfile
import unittest
from pathlib import Path

from examples.create_example import create_example, create_pathway_example
from reactome_pipeline.common import dumps, read_json, read_jsonl, write_json, write_jsonl
from reactome_pipeline.extraction_audit import audit_extraction, compare_reaction
from reactome_pipeline.normalize import normalize, normalize_one
from reactome_pipeline.review import prepare_classes


class ExtractionAuditTests(unittest.TestCase):
    def test_omissions_extras_duplicates_and_modified_fields_are_detected(self):
        raw = read_json(Path(__file__).parents[1] / "examples/input/reactions.json")[1]
        original, _ = normalize_one(raw)
        _, _, issues = compare_reaction(raw, original)
        self.assertEqual(issues, [])
        modified = copy.deepcopy(original)
        modified["inputs"].pop()
        modified["catalysts"] = []
        modified["outputs"][0]["stoichiometry"] = 999
        modified["outputs"][0]["compartments"] = ["invented"]
        modified["regulators"][0]["regulation_class"] = "PositiveRegulation"
        checks, empty, issues = compare_reaction(raw, modified)
        missing = [x for x in issues if x["type"] == "missing_participant"]
        self.assertEqual({x["role"] for x in missing}, {"inputs", "catalysts"})
        self.assertEqual({x["field"] for x in issues if x["type"] == "participant_field_mismatch"},
                         {"stoichiometry", "compartments", "regulation_class"})
        self.assertEqual(next(x["reason"] for x in empty if x["role"] == "catalysts"), "extraction_omission")
        modified = copy.deepcopy(original)
        modified["inputs"].append(copy.deepcopy(modified["inputs"][0]))
        modified["outputs"].append({"source_pointer": "/invented", "entity_id": "invented"})
        _, _, issues = compare_reaction(raw, modified)
        self.assertEqual({x["type"] for x in issues}, {"duplicate_participant", "unexpected_participant"})

    def test_empty_source_states_and_unresolved_reference_are_distinct(self):
        raw = {"dbId": 1, "stId": "R1", "inputs": None, "outputs": [],
               "catalystActivity": [{"dbId": 2, "schemaClass": "CatalystActivity"}]}
        row, _ = normalize_one(raw)
        _, empty, issues = compare_reaction(raw, row)
        self.assertEqual(issues, [])
        reasons = {x["role"]: x["reason"] for x in empty}
        self.assertEqual(reasons, {"inputs": "source_field_null", "outputs": "source_field_empty_list",
            "catalysts": "source_role_present_without_participants", "regulators": "source_field_absent"})
        raw["inputs"] = [{"element": {"@ref": "missing"}}]
        row["inputs"] = []
        _, empty, issues = compare_reaction(raw, row)
        self.assertIn("unresolved_source_reference", {x["type"] for x in issues})
        self.assertEqual(next(x["reason"] for x in empty if x["role"] == "inputs"), "unresolved_source_reference")

    def test_source_aliases_and_nonroot_roles_cannot_hide_omissions(self):
        raw = {"dbId": 1, "stId": "R1", "input": [{"@ref": "X"}],
               "output": [{"@id": "X", "dbId": 2, "stId": "X", "displayName": "X",
                           "inputFor": [{"dbId": 99, "inputs": [{"stId": "UNRELATED"}]}]}]}
        row, _ = normalize_one(raw)
        checks, _, issues = compare_reaction(raw, row)
        self.assertEqual(issues, [])
        self.assertEqual(sum(x["source_participants"] for x in checks), 2)
        raw["inputs"] = [{"stId": "DIFFERENT", "displayName": "different"}]
        _, _, issues = compare_reaction(raw, row)
        self.assertIn("multiple_role_aliases", {x["type"] for x in issues})

    def test_whole_reaction_omission_and_empty_roles_are_reported(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            create_example(root)
            normalize(root)
            original = audit_extraction(root)
            self.assertTrue(original["source_comparison_passed"])
            self.assertEqual(original["role_checks"], 24)
            self.assertEqual(original["empty_roles"]["inputs"], {"source_field_absent": 1})
            rows = list(read_jsonl(root / "normalized/reactions.jsonl"))
            write_jsonl(root / "normalized/reactions.jsonl", rows[1:])
            changed = audit_extraction(root)
            self.assertFalse(changed["source_comparison_passed"])
            self.assertEqual(changed["missing_normalized_db_ids"], [1])

    def test_class_selection_is_frozen_and_reproducible(self):
        class Client:
            base = "https://example.invalid/ContentService"
            def version(self): return "synthetic-example"
            def url(self, path, params=None): return self.base + path
            def get(self, *args, **kwargs): raise AssertionError("All selected records are cached")
        with tempfile.TemporaryDirectory() as folder:
            source, selected = Path(folder) / "source", Path(folder) / "selected"
            create_example(source)
            manifest = read_json(source / "snapshot.json")
            manifest["base_url"] = Client.base
            write_json(source / "snapshot.json", manifest)
            self.assertTrue(prepare_classes(source, selected, Client())["complete"])
            frozen = (selected / "snapshot.json").read_bytes()
            self.assertTrue(prepare_classes(source, selected, Client())["complete"])
            self.assertEqual(frozen, (selected / "snapshot.json").read_bytes())
            self.assertEqual({x["dbId"] for x in read_json(selected / "snapshot.json")["inventory"]}, {1, 5})

    def test_frozen_real_subclasses_pass_and_audit_repeats_byte_identically(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            create_pathway_example(root, "extraction_snapshot")
            normalize(root)
            result = audit_extraction(root)
            self.assertTrue(result["source_comparison_passed"])
            self.assertTrue(result["input_complete"])
            self.assertEqual(result["source_participants"], 24)
            self.assertEqual({row["schema_class"] for row in result["schema_class_coverage"]},
                             {"Reaction", "BlackBoxEvent", "Polymerisation", "Depolymerisation", "FailedReaction", "CellDevelopmentStep"})
            before = {p.name: p.read_bytes() for p in (root / "audit").glob("extraction_*")}
            audit_extraction(root)
            for name, body in before.items():
                self.assertEqual(body, (root / "audit" / name).read_bytes())


if __name__ == "__main__":
    unittest.main()
