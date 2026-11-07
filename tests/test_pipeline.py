import copy
import gzip
import itertools
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

from examples.create_example import create_example, create_real_example, create_pathway_example
from reactome_pipeline.audit import audit
from reactome_pipeline.common import digest, dumps, read_json, read_jsonl, write_json, write_jsonl
from reactome_pipeline.currency import load_currency
from reactome_pipeline.fetch import download, inventory
from reactome_pipeline.graph import build_graph, components
from reactome_pipeline.normalize import normalize, normalize_one
from reactome_pipeline.paths import enumerate_chains
from reactome_pipeline.coverage import input_coverage
from reactome_pipeline.review import prepare, packet


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
                self.assertEqual(row["text"].count("closed-loop connection:"), len(row["closing_edges"]))
                self.assertNotIn("returns to a reaction already in this chain", row["text"])
            else:
                self.assertEqual(row["closing_edges"], [])
                self.assertNotIn("closed-loop connection:", row["text"])
            self.assertIn("input_coverage", row)
            self.assertNotIn("Inputs matched to earlier outputs", row["text"])
            self.assertNotIn("Inputs not matched to earlier outputs", row["text"])
            self.assertNotIn("Input coverage:", row["text"])
            self.assertNotIn("cycle_blocked", row["text"])
            self.assertNotIn("R-DEMO-", row["text"])

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
            self.assertIn(', which is a catalyst for "mTORC2 phosphorylates SGK1".', records[0]["text"])

    def test_prefix_inputs_include_distant_sources_currency_and_missing_ids(self):
        def entity(eid, pointer, count=1):
            return {"entity_id": eid, "source_pointer": pointer, "stoichiometry": count}
        rows = {
            "A": {"inputs": [entity("later", "/inputs/0")],
                  "outputs": [entity("ATP", "/outputs/0"), entity("X", "/outputs/1", 2)]},
            "B": {"inputs": [entity("X", "/inputs/0")], "outputs": [entity("X", "/outputs/0")]},
            "C": {"inputs": [entity("ATP", "/inputs/0"), entity("X", "/inputs/1", 10),
                             entity(None, "/inputs/2"), entity("X-other-compartment", "/inputs/3"),
                             entity("later", "/inputs/4")], "outputs": [entity("later", "/outputs/0")]}}
        steps = input_coverage(["A", "B", "C"], rows)
        self.assertEqual(steps[0]["covered_inputs"], [])
        self.assertEqual([x["entity_id"] for x in steps[2]["covered_inputs"]], ["ATP", "X"])
        providers = steps[2]["covered_inputs"][1]["providers"]
        self.assertEqual([p["reaction_id"] for p in providers], ["A", "B"])
        self.assertEqual(providers[0]["stoichiometry"], 2)
        self.assertEqual(steps[2]["covered_inputs"][1]["stoichiometry"], 10)
        self.assertEqual([x["reason"] for x in steps[2]["uncovered_inputs"]],
                         ["missing_stable_identifier", "no_earlier_output_match", "no_earlier_output_match"])

    def test_parallel_evidence_and_inhibitor_text_survive_enumeration(self):
        self.graph()
        rows = list(read_jsonl(self.root / "normalized" / "reactions.jsonl"))
        inhibitor = copy.deepcopy(rows[3]["regulators"][1])
        inhibitor.update({"entity_id": "R-DEMO-I", "kind": "negative_regulation"})
        rows[3]["inputs"].append(inhibitor)
        write_jsonl(self.root / "normalized" / "reactions.jsonl", rows)
        build_graph(self.root)
        enumerate_chains(self.root)
        chain = next(c for c in self.chains() if c["reaction_ids"] == ["R-DEMO-5", "R-DEMO-4"])
        self.assertEqual({x["kind"] for x in chain["transitions"][0]["evidence"]}, {"negative_regulation", "material_flow"})
        self.assertIn('produces I [cytosol], which is an inhibitor of "B is converted to P".', chain["text"])
        self.assertIn('produces I [cytosol], which is an input to "B is converted to P".', chain["text"])

    def test_chain_and_currency_audit_denominators(self):
        self.graph()
        enumerate_chains(self.root)
        report = audit(self.root)
        chains = self.chains()
        metrics = report["chain_metrics"]
        self.assertEqual(metrics["total_chains"], len(chains))
        self.assertEqual(metrics["lengths"]["all"]["count"], 8)
        self.assertEqual(metrics["sink_fraction"], 6 / 8)
        self.assertEqual(metrics["reaction_coverage"]["total_appearances"], sum(len(c["reaction_ids"]) for c in chains))
        self.assertEqual(metrics["reaction_coverage"]["covered_reactions"], 6)
        currency = report["currency"]
        self.assertEqual(currency["reactions_losing_all_downstream"], ["R-DEMO-6"])
        self.assertEqual(currency["reactions_losing_all_upstream"], [])
        self.assertEqual(currency["reactions_becoming_isolated"], ["R-DEMO-6"])
        atp = next(x for x in currency["molecules"] if x["identifier"] == "CHEBI:30616")
        self.assertEqual(atp["excluded_evidence_edges"], 2)
        self.assertEqual(atp["affected_reactions_becoming_isolated"], ["R-DEMO-6"])

    def test_pathway_membership_numeric_refs_and_review_regeneration(self):
        class FakeClient:
            base = "https://example.invalid/ContentService"
            def version(self): return "synthetic-example"
            def url(self, path, params=None): return self.base + path
            def get(self, path, params=None):
                obj = ([{"dbId": 1}, 2, {"dbId": 999, "schemaClass": "Pathway"}]
                       if path.endswith("containedEvents") else {"stId": "P", "schemaClass": "Pathway", "displayName": "Case"})
                return dumps(obj).encode(), obj
        manifest = read_json(self.root / "snapshot.json")
        manifest["base_url"] = FakeClient.base
        write_json(self.root / "snapshot.json", manifest)
        case = self.root / "case"
        self.assertTrue(prepare(self.root, case, "P", FakeClient())["complete"])
        self.assertEqual(read_json(case / "snapshot.json")["expected_reactions"], 2)
        normalize(case)
        build_graph(case)
        summary = packet(case)
        self.assertTrue(summary["membership_complete"])
        self.assertEqual(summary["source_pointer_checks"], summary["source_pointer_passes"])
        manual = case / "review" / "manual_review.csv"
        content = manual.read_text().replace("pending", "reviewed", 1)
        manual.write_text(content)
        summary = packet(case)
        self.assertEqual(manual.read_text(), content)
        self.assertEqual(summary["manual_review_status_counts"]["reviewed"], 1)

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

    def test_complete_frozen_pathway_feedback_compartments_and_chain_counts(self):
        with tempfile.TemporaryDirectory() as other:
            create_pathway_example(other)
            self.assertTrue(normalize(other)["complete"])
            self.assertEqual(build_graph(other)["evidence_edges"], 11)
            edges = list(read_jsonl(Path(other) / "graph/edges.jsonl"))
            self.assertTrue(any(e["source"] == e["target"] == "R-BTA-174394" and e["kind"] == "negative_regulation" for e in edges))
            self.assertFalse(any("R-BTA-8943279" in (e["source"], e["target"]) for e in edges))
            status = enumerate_chains(other)
            self.assertTrue(status["complete"])
            self.assertEqual((status["sink_chains"], status["cycle_stopped_chains"]), (1, 9))
            report = audit(other)
            self.assertEqual(report["chain_metrics"]["reaction_coverage"]["total_appearances"], 45)
            summary = packet(other)
            self.assertEqual((summary["source_pointer_checks"], summary["source_pointer_passes"]), (40, 40))
            self.assertEqual(summary["manual_review_status_counts"], {"pending": 18})

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
