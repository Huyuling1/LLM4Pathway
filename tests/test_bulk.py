import json
import io
import tempfile
import unittest
from pathlib import Path

from examples.create_example import create_example
from reactome_pipeline.bulk import seed
from reactome_pipeline.catalog import contained_reactions, freeze_catalog
from reactome_pipeline.common import digest, file_sha, read_json, read_jsonl, write_json, write_jsonl
from reactome_pipeline.fetch import download
from reactome_pipeline.graph import build_graph
from reactome_pipeline.normalize import normalize
from reactome_pipeline.organize import export_pathway, index_dataset, pathway_chains, show_pathway, export_full_chains
from reactome_pipeline.paths import enumerate_chains


class BulkTests(unittest.TestCase):
    def test_shared_subpathway_membership_and_containment_cycle(self):
        paths = {p: {} for p in ('P', 'Q', 'S')}
        result = contained_reactions(paths, {'P': {'S'}, 'Q': {'S', 'R2'}, 'S': {'R1'}}, {'R1', 'R2'})
        self.assertEqual(result, {'P': {'R1'}, 'Q': {'R1', 'R2'}, 'S': {'R1'}})
        with self.assertRaisesRegex(ValueError, 'cycle'):
            contained_reactions(paths, {'P': {'Q'}, 'Q': {'P'}}, {'R1'})

    def test_small_component_schedule_work_budget_and_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ['resumed', 'whole']:
                create_example(root / name)
                normalize(root / name)
                build_graph(root / name)
            status = enumerate_chains(root / 'resumed', max_steps=1, small_components_first=True)
            self.assertFalse(status['enumeration_complete'])
            self.assertEqual(status['stop_reason'], 'dfs_step_budget')
            records = list(read_jsonl(root / 'resumed/paths/complete_chains.jsonl'))
            self.assertEqual(records[0]['reaction_ids'], ['R-DEMO-6'])
            self.assertEqual((root / 'resumed/chains.txt').read_text(), records[0]['text'] + '\n')
            with self.assertRaises(ValueError):
                enumerate_chains(root / 'resumed')
            while not enumerate_chains(root / 'resumed', max_steps=3, small_components_first=True)['complete']:
                pass
            enumerate_chains(root / 'whole', small_components_first=True)
            for name in ['complete_chains.jsonl', 'cycle_stopped_chains.jsonl']:
                self.assertEqual((root / 'resumed/paths' / name).read_bytes(), (root / 'whole/paths' / name).read_bytes())
            self.assertEqual((root / 'resumed/chains.txt').read_bytes(), (root / 'whole/chains.txt').read_bytes())
            rows = [r for name in ['complete_chains.jsonl', 'cycle_stopped_chains.jsonl']
                    for r in read_jsonl(root / 'whole/paths' / name)]
            self.assertEqual((root / 'whole/chains.txt').read_text(), '\n\n'.join(r['text'] for r in rows) + '\n')

    def test_view_keeps_external_chain_steps_boundary_edges_and_no_duplicate_indexing(self):
        with tempfile.TemporaryDirectory() as directory:
            dataset = Path(directory)
            root = dataset / 'global'
            create_example(root)
            normalize(root)
            build_graph(root)
            entry = {'id': 'P', 'display_name': '../Same/name', 'species': [{'tax_id': '1', 'name': 'Test species'}]}
            hierarchy = {'reaction_inventory_sha256': read_json(root / 'snapshot.json')['inventory_sha256'],
                         'pathways': {'P': entry}}
            write_json(dataset / 'pathway_hierarchy.json', hierarchy)
            write_jsonl(dataset / 'catalog/pathway_membership.jsonl', [{'pathway_id': 'P', 'reaction_ids': ['R-DEMO-2']}])
            write_json(dataset / 'catalog/status.json', {'complete': True,
                'hierarchy_sha256': file_sha(dataset / 'pathway_hierarchy.json'),
                'membership_sha256': file_sha(dataset / 'catalog/pathway_membership.jsonl')})
            index_dataset(dataset)
            pending = export_pathway(dataset, 'P')
            pending_manifest = read_json(dataset / pending['locations'][0] / 'manifest.json')
            self.assertFalse(pending_manifest['enumeration_complete'])
            self.assertEqual(pending_manifest['chain_count'], 0)
            enumerate_chains(root)
            index_dataset(dataset)
            first = (dataset / 'pathway_index.jsonl').read_bytes()
            index_dataset(dataset)
            self.assertEqual(first, (dataset / 'pathway_index.jsonl').read_bytes())
            result = export_pathway(dataset, 'P')
            out = dataset / result['locations'][0]
            self.assertTrue(out.resolve().is_relative_to(dataset.resolve()))
            self.assertEqual({p.name for p in out.iterdir()}, {'manifest.json', 'reactions.jsonl', 'edges.jsonl', 'chain_refs.jsonl'})
            references = list(read_jsonl(out / 'chain_refs.jsonl'))
            records = list(pathway_chains(dataset, 'P'))
            global_records = {r['id']: r for filename in ('complete_chains.jsonl', 'cycle_stopped_chains.jsonl')
                              for r in read_jsonl(root / 'paths' / filename)}
            self.assertGreater(len(records), 0)
            for reference, (record, members) in zip(references, records):
                self.assertEqual(reference, {'chain_id': record['id'], 'member_reaction_ids': ['R-DEMO-2']})
                self.assertEqual(record['text'], global_records[record['id']]['text'])
                self.assertEqual(record['reaction_ids'], global_records[record['id']]['reaction_ids'])
                self.assertEqual(members, ['R-DEMO-2'])
                self.assertGreater(len(record['reaction_ids']), 1)
            stream = io.StringIO()
            self.assertEqual(show_pathway(dataset, 'P', stream), len(records))
            self.assertEqual(stream.getvalue(), '\n\n'.join(r['text'] for r, _ in records) + '\n')
            with tempfile.TemporaryDirectory() as standalone:
                export_full_chains(dataset, 'P', standalone)
                self.assertEqual((Path(standalone) / 'chains.txt').read_text(), stream.getvalue())
            with self.assertRaises(ValueError):
                export_full_chains(dataset, 'P', out)
            # Upgrade only unchanged copies previously emitted by the old exporter.
            (out / 'chains.txt').write_text(stream.getvalue())
            manifest = read_json(out / 'manifest.json')
            manifest['files']['chains.txt'] = file_sha(out / 'chains.txt')
            write_json(out / 'manifest.json', manifest)
            export_pathway(dataset, 'P')
            self.assertFalse((out / 'chains.txt').exists())
            (out / 'chains.txt').write_text('user notes')
            with self.assertRaisesRegex(ValueError, 'modified or has no receipt'):
                export_pathway(dataset, 'P')
            self.assertEqual((out / 'chains.txt').read_text(), 'user notes')
            edges = list(read_jsonl(out / 'edges.jsonl'))
            self.assertTrue(edges)
            self.assertTrue(all(e['view_scope'] == 'boundary' for e in edges))
            self.assertTrue(read_json(out / 'manifest.json')['enumeration_complete'])
            self.assertTrue(any('unassigned' in r['pathway_id'] for r in read_jsonl(dataset / 'pathway_index.jsonl')))
            write_json(dataset / 'pathway_hierarchy.json', {})
            with self.assertRaises((ValueError, KeyError)):
                index_dataset(dataset)

    def test_catalog_includes_top_level_fallback_and_reuses_source_receipts(self):
        class FakeClient:
            base = 'https://example.invalid/ContentService'
            calls = 0
            def version(self): return 'synthetic-example'
            def url(self, path, params=None): return self.base + path + json.dumps(params, sort_keys=True)
            def get(self, path, params=None):
                self.calls += 1
                if path.endswith('/count'): obj = 2
                elif path.endswith('/min'):
                    obj = [{'stId': 'P', 'dbId': 100, 'displayName': 'P'}, {'stId': 'Q', 'dbId': 101, 'displayName': 'Q'}]
                elif path.endswith('/species/main'): obj = [{'taxId': '1', 'displayName': 'Species'}]
                elif '/eventsHierarchy/' in path:
                    obj = [{'stId': 'P', 'children': [{'stId': 'R-DEMO-1'}]}]
                else:
                    obj = {'stId': 'Q', 'schemaClass': 'TopLevelPathway',
                           'hasEvent': [{'dbId': 2}], 'species': [{'taxId': '2', 'displayName': 'Other'}]}
                return json.dumps(obj).encode(), obj
        with tempfile.TemporaryDirectory() as directory:
            dataset = Path(directory)
            create_example(dataset / 'global')
            client = FakeClient()
            result = freeze_catalog(dataset, client)
            self.assertTrue(result['complete'])
            self.assertEqual(result['pathways_queried_separately'], 1)
            members = {r['pathway_id']: r['reaction_ids'] for r in read_jsonl(dataset / 'catalog/pathway_membership.jsonl')}
            self.assertEqual(members, {'P': ['R-DEMO-1'], 'Q': ['R-DEMO-2']})
            calls = client.calls
            freeze_catalog(dataset, client)
            self.assertEqual(calls, client.calls)

    def test_final_release_network_error_preserves_failure_status(self):
        class FakeClient:
            base = 'https://example.invalid/ContentService'
            versions = 0
            def version(self):
                self.versions += 1
                if self.versions > 1:
                    raise RuntimeError('HTTP 521')
                return 'synthetic-example'
            def url(self, path, params=None): return self.base + path
        with tempfile.TemporaryDirectory() as directory:
            create_example(directory)
            with self.assertRaisesRegex(RuntimeError, '521'):
                download(directory, FakeClient())
            status = read_json(Path(directory) / 'download_status.json')
            self.assertFalse(status['complete'])
            self.assertFalse(status['release_check_complete'])
            self.assertEqual(status['downloaded'], 6)
            self.assertEqual(status['stage'], 'stopped')

    def test_seed_preserves_snapshot_and_does_not_copy_old_processing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            create_example(root / 'source')
            normalize(root / 'source')
            seed(root / 'source', root / 'dataset')
            self.assertEqual(read_json(root / 'source/snapshot.json'), read_json(root / 'dataset/global/snapshot.json'))
            self.assertFalse((root / 'dataset/global/normalized').exists())
            seed(root / 'source', root / 'dataset')
            self.assertEqual(len(list((root / 'dataset/global/raw/reactions').glob('*.json.gz'))), 6)


if __name__ == '__main__':
    unittest.main()
