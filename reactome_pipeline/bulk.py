"""Resumable full-snapshot batch driver; network errors stop with saved progress."""

import argparse
import os
import shutil
import sqlite3
import time
from pathlib import Path

from .common import file_sha, read_json, write_json
from .fetch import Client, download, load_snapshot
from .locking import SnapshotLock


def process(dataset, pathway_ids=(), rounds=None, max_steps=1000000, max_seconds=300,
            stop_after=10000, min_free_gb=10):
    """Process only complete downloads; bounded work slices resume full simple paths."""
    from .normalize import normalize
    from .extraction_audit import audit_extraction
    from .graph import build_graph
    from .paths import enumerate_chains, RECORD_VERSION
    from .audit import audit
    from .organize import index_dataset, export_pathway
    dataset = Path(dataset)
    root = dataset / "global"
    progress_path = dataset / "batch_status.json"
    code = {p.name: file_sha(p) for p in sorted(Path(__file__).parent.glob("*.py"))}
    config = {"chain_record_version": RECORD_VERSION, "start_order": "weak-component-size-v1",
              "max_dfs_steps_per_batch": max_steps, "max_seconds_per_batch": max_seconds,
              "max_new_records_per_batch": stop_after, "min_free_gb": min_free_gb,
              "export_pathway_ids": list(pathway_ids), "code_sha256": code}
    config_path = dataset / "processing_config.json"
    if config_path.exists() and read_json(config_path) != config:
        raise ValueError("Processing code/config changed; retain this run and use a new processing dataset")
    write_json(config_path, config)
    current = {"stage": "checking_inputs", "complete": False, "process_id": os.getpid()}
    def report(stage, **details):
        current.update(stage=stage, **details)
        write_json(progress_path, current)
        print(current, flush=True)
    def disk_check():
        if shutil.disk_usage(dataset).free < min_free_gb * 1024 ** 3:
            raise RuntimeError("Free disk space below configured reserve; saved checkpoints can be resumed")
    try:
        with SnapshotLock(root), SnapshotLock(dataset / "views"):
            if not read_json(root / "download_status.json")["complete"]:
                raise ValueError("Reaction download is incomplete; resume fetch first")
            if not read_json(dataset / "catalog/status.json")["complete"]:
                raise ValueError("Pathway catalog is incomplete; resume catalog first")
            stages_file = dataset / "processed_stages.json"
            stages = read_json(stages_file) if stages_file.exists() else {}
            for name, action in (("normalize", lambda: normalize(root)),
                                 ("extraction_audit", lambda: audit_extraction(root)),
                                 ("graph", lambda: build_graph(root)),
                                 ("graph_audit", lambda: audit(root))):
                if stages.get(name):
                    continue
                disk_check()
                report(name)
                result = action()
                if name == "extraction_audit" and not (result["source_comparison_passed"] and result["input_complete"]):
                    raise ValueError("Raw-versus-normalized audit failed; inspect global/audit/extraction_issues.jsonl")
                stages[name] = True
                write_json(stages_file, stages)
            report("indexing")
            index_dataset(dataset)
            for pid in pathway_ids:
                export_pathway(dataset, pid)
            batch = 0
            while rounds is None or batch < rounds:
                disk_check()
                report("enumerating", batches_this_run=batch)
                result = enumerate_chains(root, stop_after=stop_after, max_steps=max_steps,
                                          max_seconds=max_seconds, small_components_first=True)
                batch += 1
                report("auditing_chains", enumeration_complete=result["enumeration_complete"],
                       chain_count=result["sink_chains"] + result["cycle_stopped_chains"], dfs_steps=result["dfs_steps"])
                audit(root)
                index_dataset(dataset)
                for pid in pathway_ids:
                    export_pathway(dataset, pid)
                if result["complete"]:
                    report("complete", complete=True)
                    return current
                if result["interrupted"]:
                    report("interrupted")
                    return current
            report("checkpointed")
            return current
    except (Exception, KeyboardInterrupt) as exc:
        report("stopped", error=str(exc), complete=False)
        raise


def wait_for_downloads(dataset, jobs_file):
    """Wait for the already-launched workers; never silently restart failed downloads."""
    dataset = Path(dataset)
    jobs = read_json(jobs_file)
    files = {"fetch": dataset / "global/download_status.json", "catalog": dataset / "catalog/status.json"}
    while True:
        pending = []
        for name, path in files.items():
            if path.exists() and read_json(path).get("complete"):
                continue
            pid = jobs[name]["pid"]
            try:
                os.kill(pid, 0)
                # A dead child may remain as a zombie until its launcher reaps it.
                procstat = Path(f"/proc/{pid}/stat")
                if procstat.exists() and procstat.read_text().split(") ", 1)[1].startswith("Z"):
                    raise ProcessLookupError(pid)
            except ProcessLookupError:
                raise RuntimeError(f"{name} worker stopped before completion; inspect {jobs[name]['log']} and resume it manually")
            pending.append(name)
        write_json(dataset / "batch_status.json", {"stage": "waiting_for_downloads" if pending else "inputs_ready",
                   "pending": pending, "complete": False, "process_id": os.getpid()})
        if not pending:
            return
        time.sleep(30)


def status(dataset):
    dataset = Path(dataset)
    root = dataset / "global"
    snapshot = load_snapshot(root)
    result = {"reactome_release": snapshot["reactome_release"], "expected_reactions": snapshot["expected_reactions"]}
    if (root / "downloads.sqlite").exists():
        with sqlite3.connect((root / "downloads.sqlite").resolve().as_uri() + "?mode=ro", uri=True) as db:
            result["downloaded_with_receipts"] = db.execute("SELECT count(*) FROM receipts").fetchone()[0]
    for key, relative in (("download", "global/download_status.json"), ("catalog", "catalog/status.json"),
                          ("processing", "batch_status.json"), ("enumeration", "global/paths/status.json")):
        if (dataset / relative).exists():
            result[key] = read_json(dataset / relative)
    if (dataset / "catalog/receipts.json").exists():
        result["catalog_cached_responses"] = len(read_json(dataset / "catalog/receipts.json")["sources"])
    return result


def seed(source, dataset):
    source, dataset = Path(source), Path(dataset)
    target = dataset / "global"
    manifest = load_snapshot(source)
    if (source / "release_mismatch.json").exists():
        raise ValueError("Source has a recorded release mismatch")
    with SnapshotLock(source), SnapshotLock(target):
        if (target / "snapshot.json").exists():
            if load_snapshot(target) != manifest:
                raise ValueError("Destination belongs to another frozen snapshot")
        else:
            write_json(target / "snapshot.json", manifest)
        for src in sorted((source / "raw").rglob("*")):
            if not src.is_file():
                continue
            dst = target / src.relative_to(source)
            if dst.exists():
                if file_sha(src) != file_sha(dst):
                    raise ValueError("Existing destination raw file differs: " + str(dst))
                continue
            dst.parent.mkdir(parents=True, exist_ok=True)
            # Raw records are immutable. Atomic writes replace, never modify, an inode.
            try:
                os.link(src, dst)
            except OSError:
                shutil.copyfile(src, dst)
        if (source / "downloads.sqlite").exists() and not (target / "downloads.sqlite").exists():
            with sqlite3.connect(source / "downloads.sqlite") as src_db, sqlite3.connect(target / "downloads.sqlite") as dst_db:
                src_db.backup(dst_db)
    info = {"format_version": 1, "reactome_release": manifest["reactome_release"],
            "inventory_sha256": manifest["inventory_sha256"], "expected_reactions": manifest["expected_reactions"],
            "source_snapshot": str(source.resolve()), "global_directory": "global",
            "pathway_directory": "by_species", "status": "seeded"}
    if not (dataset / "manifest.json").exists():
        write_json(dataset / "manifest.json", info)
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("seed", "fetch", "process", "status"))
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--source")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--interval", type=float, default=0.2, help="Minimum seconds between HTTP request starts across all workers")
    parser.add_argument("--after-jobs", help="Wait for these launched download workers, stopping if either fails")
    parser.add_argument("--pathway", action="append", default=[])
    parser.add_argument("--rounds", type=int, help="Limit work batches; rerun to resume")
    parser.add_argument("--max-dfs-steps", type=int, default=1000000)
    parser.add_argument("--max-seconds", type=float, default=300)
    parser.add_argument("--stop-after", type=int, default=10000)
    parser.add_argument("--min-free-gb", type=float, default=10)
    args = parser.parse_args()
    if args.workers < 1 or args.interval < 0 or args.max_dfs_steps < 1 or args.max_seconds <= 0 or args.stop_after < 1 or args.min_free_gb < 0 or (args.rounds is not None and args.rounds < 1):
        parser.error("Invalid worker count or processing budget")
    if args.stage == "status":
        print(status(args.dataset), flush=True)
        return 0
    if args.stage == "seed":
        if not args.source:
            parser.error("seed requires --source")
        result = seed(args.source, args.dataset)
    elif args.stage == "fetch":
        root = Path(args.dataset) / "global"
        manifest = load_snapshot(root)
        with SnapshotLock(root):
            result = download(root, Client(manifest["base_url"], interval=args.interval), workers=args.workers)
    else:
        if args.after_jobs:
            try:
                wait_for_downloads(args.dataset, args.after_jobs)
            except Exception as exc:
                write_json(Path(args.dataset) / "batch_status.json", {"stage": "stopped", "error": str(exc), "complete": False})
                raise
        result = process(args.dataset, args.pathway, args.rounds, args.max_dfs_steps,
                         args.max_seconds, args.stop_after, args.min_free_gb)
    print(result, flush=True)
    return 0 if result.get("complete", True) else 2


if __name__ == "__main__":
    raise SystemExit(main())
