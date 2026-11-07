"""Independent, inspectable processing stages."""

import argparse
import sqlite3
import sys
from pathlib import Path

from .audit import audit
from .common import dumps, read_json
from .fetch import Client, DEFAULT_BASE, download, inventory, load_snapshot
from .graph import build_graph
from .normalize import normalize
from .paths import enumerate_chains
from .locking import SnapshotLock
from .extraction_audit import audit_extraction


def main(argv=None):
    parser = argparse.ArgumentParser(description="Deterministic Reactome reaction graphs and maximal simple chains")
    parser.add_argument("stage", choices=("discover", "fetch", "normalize", "graph", "paths", "audit", "audit-extraction", "status", "run"))
    parser.add_argument("--out", required=True, help="Snapshot directory; reuse it to resume")
    parser.add_argument("--base-url", default=DEFAULT_BASE)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--interval", type=float, default=0.2, help="Minimum seconds between HTTP request starts across all workers")
    parser.add_argument("--timeout", type=float, default=45)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--max-failures", type=int, default=8, help="Stop scheduling requests after consecutive failures")
    parser.add_argument("--page-size", type=int, default=20000)
    parser.add_argument("--currency", help="Alternative versioned currency JSON")
    parser.add_argument("--limit", type=int, help="Fetch at most this many missing records; snapshot remains partial")
    parser.add_argument("--allow-partial", action="store_true", help="Explicitly process downloaded subset; output retains partial status")
    parser.add_argument("--stop-after", type=int, help="Checkpoint after this many NEW chain records; rerun to continue")
    parser.add_argument("--checkpoint-every", type=int, default=1000, help="DFS operations per checkpoint")
    parser.add_argument("--max-dfs-steps", type=int, help="Pause enumeration after this many additional DFS operations")
    parser.add_argument("--max-seconds", type=float, help="Pause enumeration after this many seconds of traversal")
    parser.add_argument("--small-components-first", action="store_true", help="Visit independent small graph components first")
    args = parser.parse_args(argv)
    if args.workers < 1 or args.interval < 0 or args.timeout <= 0 or args.retries < 0 or args.max_failures < 1 or (args.limit is not None and args.limit < 1):
        parser.error("Invalid network or fetch limit argument")
    try:
        if args.stage == "status":
            print(dumps(snapshot_status(args.out)))
            return 0
        with SnapshotLock(args.out):
            return run_stage(args)
    except (ValueError, RuntimeError, OSError, sqlite3.Error) as exc:
        print(str(exc), file=sys.stderr)
        return 1


def snapshot_status(root):
    root = Path(root)
    manifest = load_snapshot(root)
    result = {"reactome_release": manifest["reactome_release"], "expected_reactions": manifest["expected_reactions"]}
    if (root / "downloads.sqlite").exists():
        with sqlite3.connect((root / "downloads.sqlite").resolve().as_uri() + "?mode=ro", uri=True) as db:
            result["downloaded_with_receipts"] = db.execute("SELECT count(*) FROM receipts").fetchone()[0]
    for key, filename in (("download", "download_status.json"), ("normalized", "normalized/status.json"), ("graph", "graph/status.json"), ("paths", "paths/status.json"), ("enumeration_progress", "paths/progress.json"), ("extraction_audit", "audit/extraction_report.json")):
        if (root / filename).exists():
            result[key] = read_json(root / filename)
    return result


def run_stage(args):
    client = Client(args.base_url, args.timeout, args.retries, args.interval)
    try:
        if args.stage in ("discover", "run"):
            result = inventory(args.out, client, args.page_size)
            print(dumps({k: v for k, v in result.items() if k != "inventory"}), flush=True)
        if args.stage in ("fetch", "run"):
            result = download(args.out, client, args.workers, args.max_failures, args.limit)
            print(dumps(result), flush=True)
            if not result["complete"] and not args.allow_partial:
                return 2
        if args.stage in ("normalize", "run"):
            result = normalize(args.out, args.allow_partial)
            print(dumps(result), flush=True)
            if args.stage == "normalize" and not result["complete"]:
                return 2
            if args.stage == "run":
                checks = audit_extraction(args.out)
                print(dumps(checks), flush=True)
                if not checks["source_comparison_passed"]:
                    return 2
        if args.stage in ("graph", "run"):
            print(dumps(build_graph(args.out, args.currency)), flush=True)
        if args.stage == "run":
            print(dumps(audit(args.out)), flush=True)
        if args.stage in ("paths", "run"):
            result = enumerate_chains(args.out, args.stop_after, args.checkpoint_every,
                                      args.max_dfs_steps, args.max_seconds, args.small_components_first)
            print(dumps(result), flush=True)
            if args.stage == "run":
                audit(args.out)
            if not result["complete"]:
                return 2
        if args.stage == "audit":
            print(dumps(audit(args.out)), flush=True)
        if args.stage == "audit-extraction":
            result = audit_extraction(args.out)
            print(dumps(result), flush=True)
            return 0 if result["source_comparison_passed"] and result["input_complete"] else 2
    except (ValueError, RuntimeError, OSError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0
