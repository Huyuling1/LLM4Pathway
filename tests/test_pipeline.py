import copy
import gzip
import itertools
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

from examples.create_example import create_example, create_real_example
from reactome_pipeline.audit import audit
from reactome_pipeline.common import digest, dumps, read_json, read_jsonl, write_json, write_jsonl
from reactome_pipeline.currency import load_currency
from reactome_pipeline.fetch import download, inventory
from reactome_pipeline.graph import build_graph, components
from reactome_pipeline.normalize import normalize, normalize_one
from reactome_pipeline.paths import enumerate_chains


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        create_example(self.root)

    def tearDown(self):
        self.temp.cleanup()

    def graph(self):
        normalize(self.root)
        return build_graph(self.root)

    def chains(self):
        return list(read_jsonl(self.root / "paths" / "complete_chains.jsonl")) + list(read_jsonl(self.root / "paths" / "cycle_stopped_chains.jsonl"))

    def test_all_roles_and_currency_exclusion(self):
        status = self.graph()
        edges = list(read_jsonl(self.root / "graph" / "edges.jsonl"))
        actual = {(x["source"], x["target"], x["kind"]) for x in edges}
        self.assertIn(("R-DEMO-5", "R-DEMO-2", "catalyst_dependency"), actual)
        self.assertIn(("R-DEMO-5", "R-DEMO-4", "negative_regulation"), actual)
        self.assertIn(("R-DEMO-1", "R-DEMO-4", "positive_regulation"), actual)
        self.assertIn(("R-DEMO-2", "R-DEMO-3", "material_flow"), actual)
        self.assertFalse(any(x["entity_id"] == "R-DEMO-ATP" for x in edges))
        reactions = list(read_jsonl(self.root / "normalized" / "reactions.jsonl"))
        self.assertIn("ATP [cytosol]", reactions[1]["text"])
        self.assertIn("2 A [cytosol]", reactions[1]["text"])
        self.assertTrue(reactions[1]["is_inferred"])
        self.assertEqual(status["reactions"], 6)
        self.assertTrue(status["input_complete"])

    def test_forward_references_and_root_scope(self):
        raw = {"dbId": 7, "stId": "R7", "schemaClass": "Reaction", "inputs": [{"element": {"@ref": "a"}, "stoichiometry": 3}],
               "outputs": [{"element": {"@id": "a", "dbId": 8, "stId": "A", "displayName": "A", "inputFor": [{"element": {"dbId": 9, "stId": "OTHER", "inputs": [{"element": {"dbId": 10, "stId": "OTHER-INPUT"}}]}}]}}]}
        row, issues = normalize_one(raw)
        self.assertEqual(row["inputs"][0]["entity_id"], "A")
        self.assertEqual(row["inputs"][0]["stoichiometry"], 3)
        self.assertEqual(len(row["inputs"]), 1)
        self.assertFalse(issues)

    def test_flat_fields_and_missing_names(self):
        row, issues = normalize_one({"dbId": 1, "input": [{"stId": "A"}], "output": [], "species": {"dbId": 2, "taxId": "9606", "displayName": "Homo sapiens"}})
        self.assertEqual(row["id"], "db:1")
        self.assertIn("[not recorded]", row["text"])
        self.assertTrue(any(x["type"] == "participant_without_name" for x in issues))
        self.assertEqual(row["species"][0]["tax_id"], "9606")

    def test_complex_components_are_not_matched(self):
        self.graph()
        rows = list(read_jsonl(self.root / "normalized" / "reactions.jsonl"))
        rows[0]["outputs"][0]["entity_id"] = "COMPLEX-CONTAINING-A"
        write_jsonl(self.root / "normalized" / "reactions.jsonl", rows)
        build_graph(self.root)
        edges = list(read_jsonl(self.root / "graph" / "edges.jsonl"))
        self.assertFalse(any(x["source"] == "R-DEMO-1" and x["target"] == "R-DEMO-2" and x["kind"] == "material_flow" for x in edges))

    def test_currency_matches_reference_not_name(self):
        self.graph()
        rows = list(read_jsonl(self.root / "normalized" / "reactions.jsonl"))
        for row in rows:
            for role in ("inputs", "outputs", "regulators"):
                for item in row[role]:
                    if item["entity_id"] == "R-DEMO-ATP":
                        item["chebi_ids"] = []
        write_jsonl(self.root / "normalized" / "reactions.jsonl", rows)
        build_graph(self.root)
        self.assertTrue(any(x["entity_id"] == "R-DEMO-ATP" for x in read_jsonl(self.root / "graph" / "edges.jsonl")))

    def test_resume_is_identical_to_uninterrupted_enumeration(self):
        self.graph()
        status = enumerate_chains(self.root, stop_after=2, checkpoint_every=1)
        self.assertFalse(status["enumeration_complete"])
        self.assertTrue(enumerate_chains(self.root)["complete"])
        first = {name: (self.root / "paths" / name).read_bytes() for name in ("complete_chains.jsonl", "cycle_stopped_chains.jsonl")}
        with tempfile.TemporaryDirectory() as other:
            create_example(other)
            normalize(other)
            build_graph(other)
            enumerate_chains(other)
            for name, content in first.items():
                self.assertEqual(content, (Path(other) / "paths" / name).read_bytes())
        enumerate_chains(self.root)
        for name, content in first.items():
            self.assertEqual(content, (self.root / "paths" / name).read_bytes())

    def test_cycle_records_and_no_repeated_reactions(self):
        self.graph()
        enumerate_chains(self.root)
        chains = self.chains()
        self.assertTrue(any(x["status"] == "cycle_blocked" for x in chains))
        self.assertTrue(any(x["reaction_ids"] == ["R-DEMO-6"] for x in chains))
        for row in chains:
            self.assertEqual(len(row["reaction_ids"]), len(set(row["reaction_ids"])))
            if row["status"] == "cycle_blocked":
                self.assertTrue(row["closing_edges"])
                self.assertTrue(all(edge["target"] in row["reaction_ids"] for edge in row["closing_edges"]))

    def test_curated_order_is_audited_separately(self):
        self.graph()
        report = audit(self.root)
        self.assertEqual(report["preceding_event"]["overlap"], 1)
        self.assertEqual(report["preceding_event"]["matched_negative_preceding"], 1)
        self.assertEqual(report["currency"]["excluded_evidence_edges"], 2)

    def test_partial_coverage_cannot_be_reported_complete(self):
        (self.root / "raw" / "reactions" / "6.json.gz").unlink()
        with self.assertRaises(ValueError):
            normalize(self.root)
        self.assertFalse(normalize(self.root, allow_partial=True)["complete"])
        self.assertFalse(build_graph(self.root)["input_complete"])
        status = enumerate_chains(self.root)
        self.assertTrue(status["enumeration_complete"])
        self.assertFalse(status["complete"])

    def test_graph_change_rejects_old_checkpoint(self):
        self.graph()
        enumerate_chains(self.root, stop_after=1)
        conf, _, _, _ = load_currency()
        conf["reactome_st_ids"] = ["R-DEMO-A"]
        path = self.root / "different-currency.json"
        write_json(path, conf)
        build_graph(self.root, path)
        with self.assertRaises(ValueError):
            enumerate_chains(self.root)

    def test_shifted_inventory_pages_recover_name_ties(self):
        class FakeClient:
            base = "https://example.invalid"
            def version(self): return "98"
            def get(self, path, params=None):
                if path.endswith("/count"): return b"4", 4
                size, page = params["offset"], params["page"]
                ids = ([1, 2] if page == 1 else [2, 4]) if size == 2 else [page]
                rows = [{"dbId": i, "stId": f"R{i}", "displayName": "same"} for i in ids]
                return dumps(rows).encode(), rows
        with tempfile.TemporaryDirectory() as other:
            result = inventory(other, FakeClient(), page_size=2)
            self.assertEqual({x["dbId"] for x in result["inventory"]}, {1, 2, 3, 4})

    def test_download_resume_and_checksum_validation(self):
        fixture = {row["dbId"]: row for row in read_json(Path(__file__).parents[1] / "examples" / "input" / "reactions.json")}
        class FakeClient:
            base = "https://example.invalid/ContentService"
            calls = 0
            def version(self): return "synthetic-example"
            def url(self, path, params=None): return self.base + path
            def get(self, path, params=None):
                self.calls += 1
                obj = fixture[int(path.rsplit("/", 1)[1])]
                return dumps(obj).encode(), obj
        client = FakeClient()
        (self.root / "raw" / "reactions" / "1.json.gz").unlink()
        self.assertTrue(download(self.root, client, workers=1)["complete"])
        self.assertEqual(client.calls, 1)
        self.assertTrue(download(self.root, client, workers=1)["complete"])
        self.assertEqual(client.calls, 1)
        target = self.root / "raw" / "reactions" / "1.json.gz"
        target.unlink()
        result = download(self.root, client, workers=1)
        self.assertTrue(result["complete"])
        self.assertEqual(result["downloaded"], 6)
        self.assertEqual(client.calls, 2)
        changed = copy.deepcopy(fixture[1]); changed["displayName"] = "modified"
        target.write_bytes(gzip.compress(dumps(changed).encode(), mtime=0))
        with self.assertRaises(ValueError):
            download(self.root, client, workers=1)
        with self.assertRaises(ValueError):
            normalize(self.root)

    def test_release_change_blocks_downloads(self):
        class FakeClient:
            base = "https://example.invalid/ContentService"
            def version(self): return "different-release"
        with self.assertRaises(ValueError):
            download(self.root, FakeClient())

    def test_frozen_real_reactions_form_catalyst_connection(self):
        with tempfile.TemporaryDirectory() as other:
            create_real_example(other)
            result = normalize(other)
            self.assertEqual(result["normalized"], 2)
            self.assertEqual(result["issues"], {})
            build_graph(other)
            edges = list(read_jsonl(Path(other) / "graph" / "edges.jsonl"))
            self.assertEqual(len(edges), 1)
            self.assertEqual((edges[0]["source"], edges[0]["target"], edges[0]["kind"], edges[0]["entity_id"]),
                             ("R-HSA-9980233", "R-HSA-6795290", "catalyst_dependency", "R-HSA-9980228"))
            enumerate_chains(other)
            records = list(read_jsonl(Path(other) / "paths" / "complete_chains.jsonl"))
            self.assertEqual(len(records), 1)
            self.assertIn("ATP [cytosol]", records[0]["text"])
            self.assertIn("provides a catalyst", records[0]["text"])

    def test_source_pointers_locate_original_json_values(self):
        rows = read_json(Path(__file__).parents[1] / "examples" / "input" / "reactions.json")
        real = Path(__file__).parents[1] / "examples" / "real_snapshot" / "raw" / "reactions"
        rows.extend(read_json(path) for path in sorted(real.glob("*.json.gz")))
        for raw in rows:
            normalized, _ = normalize_one(raw)
            for role in ("inputs", "outputs", "catalysts", "regulators"):
                for item in normalized[role]:
                    obj = raw
                    for part in item["source_pointer"].split("/")[1:]:
                        key = part.replace("~1", "/").replace("~0", "~")
                        obj = obj[int(key)] if isinstance(obj, list) else obj[key]
                    self.assertIsInstance(obj, dict)
                    if "@ref" not in obj:
                        self.assertEqual(obj.get("stId"), item["entity_id"])

    def test_maximal_paths_match_independent_permutation_oracle(self):
        # All directed graphs on three vertices without self-loops, plus self-loop cases.
        nodes = ["A", "B", "C"]
        possible = list(itertools.permutations(nodes, 2))
        for mask in range(64):
            edges = {edge for i, edge in enumerate(possible) if mask & (1 << i)}
            if mask % 7 == 0:
                edges.add(("A", "A"))
            adjacency = {x: sorted(b for a, b in edges if a == x) for x in nodes}
            reverse = {x: sorted(a for a, b in edges if b == x) for x in nodes}
            expected = set()
            for length in (1, 2, 3):
                for path in itertools.permutations(nodes, length):
                    seen = set(path)
                    if all((a, b) in edges for a, b in zip(path, path[1:])) and set(reverse[path[0]]) <= seen and set(adjacency[path[-1]]) <= seen:
                        expected.add(path)
            with tempfile.TemporaryDirectory() as other:
                out = Path(other) / "graph"
                out.mkdir()
                db = sqlite3.connect(out / "graph.sqlite")
                db.executescript("CREATE TABLE reactions(id,data); CREATE TABLE edges(source,target,kind,entity); CREATE TABLE entities(id,name); CREATE TABLE metadata(key,value);")
                db.executemany("INSERT INTO reactions VALUES (?,?)", ((node, dumps({"text": node})) for node in nodes))
                db.execute("INSERT INTO entities VALUES ('E','E')")
                db.executemany("INSERT INTO edges VALUES (?,?,'material_flow','E')", sorted(edges))
                db.executemany("INSERT INTO metadata VALUES (?,?)", [("graph_sha256", dumps(str(mask))), ("input_complete", "true")])
                db.commit(); db.close()
                enumerate_chains(other)
                actual = {tuple(row["reaction_ids"]) for name in ("complete_chains.jsonl", "cycle_stopped_chains.jsonl") for row in read_jsonl(Path(other) / "paths" / name)}
                self.assertEqual(expected, actual, f"mask={mask}")


if __name__ == "__main__":
    unittest.main()
