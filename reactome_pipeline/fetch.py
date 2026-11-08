"""Release-checked inventory and resumable per-reaction nested JSON downloads."""

import concurrent.futures
import gzip
import hashlib
import json
import math
import sqlite3
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from .common import digest, read_json, write_bytes, write_json, write_jsonl

DEFAULT_BASE = "https://reactome.org/ContentService"
RETRYABLE_STATUS = frozenset({403, 408, 429, 500, 502, 503, 504, 520, 521, 522, 523, 524, 525, 526})
QUERY = {"view": "nested", "includeRef": "true", "fetchIncomingRelationships": "true", "summariseReferenceEntity": "true"}


def load_snapshot(root):
    manifest = read_json(Path(root) / "snapshot.json")
    rows = manifest["inventory"]
    if digest(rows) != manifest["inventory_sha256"] or len(rows) != manifest["expected_reactions"]:
        raise ValueError("Frozen inventory checksum or count mismatch")
    if len({row["dbId"] for row in rows}) != len(rows):
        raise ValueError("Repeated dbId in frozen inventory")
    return manifest


class Client:
    def __init__(self, base=DEFAULT_BASE, timeout=45, retries=3, interval=0.2):
        self.base = base.rstrip("/")
        self.timeout, self.retries, self.interval = timeout, retries, interval
        self._lock, self._next = threading.Lock(), 0.0

    def url(self, path, params=None):
        return self.base + path + ("?" + urllib.parse.urlencode(sorted(params.items())) if params else "")

    def get(self, path, params=None, accept="application/json"):
        url = self.url(path, params)
        for attempt in range(self.retries + 1):
            with self._lock:
                now = time.monotonic()
                delay = max(0, self._next - now)
                self._next = max(now, self._next) + self.interval
            if delay:
                time.sleep(delay)
            request = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "LLM4Pathway/0.1"})
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    body = response.read()
                return body, json.loads(body)
            except (urllib.error.URLError, TimeoutError, ConnectionError, json.JSONDecodeError) as exc:
                retryable = not isinstance(exc, urllib.error.HTTPError) or exc.code in RETRYABLE_STATUS
                if not retryable or attempt == self.retries:
                    raise RuntimeError(f"GET {url}: {exc}") from exc
                delay = min(2 ** attempt, 30)
                if isinstance(exc, urllib.error.HTTPError):
                    retry_after = exc.headers.get("Retry-After", "")
                    if retry_after.isdigit():
                        delay = max(delay, min(int(retry_after), 120))
                time.sleep(delay)

    def version(self):
        return str(self.get("/data/database/version", accept="text/plain")[1])


def inventory(root, client, page_size=20000):
    root = Path(root)
    manifest_path = root / "snapshot.json"
    if manifest_path.exists():
        manifest = load_snapshot(root)
        if manifest["base_url"] != client.base:
            raise ValueError("Snapshot base URL differs from this client")
        if digest(manifest["inventory"]) != manifest["inventory_sha256"]:
            raise ValueError("Inventory checksum mismatch")
        return manifest
    if not 1 <= page_size <= 20000:
        raise ValueError("Inventory page size must be between 1 and 20000")
    release = client.version()
    _, count = client.get("/data/schema/ReactionLikeEvent/count")
    count = int(count)
    if count <= 0:
        raise ValueError("Empty ReactionLikeEvent inventory")
    by_db = {}
    # The service orders by displayName; equal names can straddle pages in an
    # unstable order. Shift page boundaries and union dbIds, then check count.
    passes = []
    for size in dict.fromkeys(max(1, page_size - shift) for shift in (0, 31, 127)):
        for page in range(1, math.ceil(count / size) + 1):
            body, rows = client.get("/data/schema/ReactionLikeEvent/min", {"page": page, "offset": size})
            if not isinstance(rows, list):
                raise ValueError("Inventory page is not an array")
            write_bytes(root / "raw" / "inventory" / f"size-{size}-page-{page:05d}.json", body)
            for row in rows:
                key = int(row["dbId"])
                value = {k: row[k] for k in ("dbId", "stId", "schemaClass", "displayName") if k in row}
                if key in by_db and by_db[key] != value:
                    raise ValueError(f"Inventory entry changed for dbId {key}")
                by_db[key] = value
            print(f"inventory page_size={size} page={page} reactions={len(by_db)}/{count}", flush=True)
        passes.append({"page_size": size, "unique_reactions": len(by_db)})
        if len(by_db) >= count:
            break
    if len(by_db) != count:
        raise ValueError(f"Inventory count mismatch: {len(by_db)} != {count}")
    if client.version() != release:
        raise ValueError("Reactome release changed during inventory discovery")
    rows = sorted(by_db.values(), key=lambda row: (row.get("stId", ""), row["dbId"]))
    manifest = {"format_version": 1, "reactome_release": release, "base_url": client.base,
                "scope": "all ReactionLikeEvent subclasses; all species; inferred and disease events included",
                "query": QUERY, "expected_reactions": count, "inventory": rows, "discovery_passes": passes,
                "inventory_sha256": digest(rows)}
    write_json(manifest_path, manifest)
    return manifest


def raw_path(root, db_id):
    return Path(root) / "raw" / "reactions" / f"{db_id}.json.gz"


def download(root, client, workers=4, max_failures=8, limit=None):
    root = Path(root)
    manifest = load_snapshot(root)
    release = manifest["reactome_release"]
    if client.base != manifest["base_url"] or client.version() != release:
        raise ValueError("API endpoint or release differs from the frozen snapshot; use a new output directory")
    if manifest["query"] != QUERY:
        raise ValueError("Query options differ from the frozen snapshot")
    records = manifest["inventory"]
    connection = sqlite3.connect(root / "downloads.sqlite")
    connection.execute("CREATE TABLE IF NOT EXISTS receipts (db_id INTEGER PRIMARY KEY, sha256 TEXT NOT NULL, url TEXT NOT NULL)")
    receipts = {r[0]: r[1:] for r in connection.execute("SELECT db_id,sha256,url FROM receipts")}
    pending = []
    completed = 0
    corrupt = []
    for record in records:
        db_id = record["dbId"]
        path = raw_path(root, db_id)
        if path.exists():
            try:
                body = gzip.decompress(path.read_bytes())
                sha = hashlib.sha256(body).hexdigest()
                obj = json.loads(body)
                if obj.get("dbId") != db_id or (record.get("stId") and obj.get("stId") != record["stId"]):
                    raise ValueError("identity mismatch")
                if db_id in receipts and receipts[db_id][0] != sha:
                    raise ValueError("raw checksum mismatch")
                url = client.url(f"/data/query/enhanced/v2/{db_id}", QUERY)
                connection.execute("INSERT OR IGNORE INTO receipts VALUES (?,?,?)", (db_id, sha, url))
                completed += 1
                continue
            except (ValueError, OSError, EOFError) as exc:
                corrupt.append({"db_id": db_id, "error": str(exc)})
        pending.append(record)
    connection.commit()
    if corrupt:
        write_jsonl(root / "raw_corruption.jsonl", corrupt)
        connection.close()
        raise ValueError("Existing raw files failed validation; see raw_corruption.jsonl")
    if limit is not None:
        pending = pending[:limit]
    errors, consecutive = [], 0
    stopped = False

    def fetch_one(record):
        db_id = record["dbId"]
        body, obj = client.get(f"/data/query/enhanced/v2/{db_id}", QUERY)
        if not isinstance(obj, dict) or obj.get("dbId") != db_id or (record.get("stId") and obj.get("stId") != record["stId"]):
            raise ValueError(f"Downloaded identity differs from inventory for {db_id}")
        sha = hashlib.sha256(body).hexdigest()
        if db_id in receipts and receipts[db_id][0] != sha:
            raise ValueError(f"Redownload differs from original receipt for {db_id}")
        write_bytes(raw_path(root, db_id), gzip.compress(body, mtime=0))
        return db_id, sha, client.url(f"/data/query/enhanced/v2/{db_id}", QUERY)

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
    tasks, position = {}, 0
    try:
        while tasks or (position < len(pending) and not stopped):
            while not stopped and position < len(pending) and len(tasks) < workers:
                item = pending[position]
                tasks[executor.submit(fetch_one, item)] = item
                position += 1
            if not tasks:
                break
            done, _ = concurrent.futures.wait(tasks, return_when=concurrent.futures.FIRST_COMPLETED)
            for task in done:
                item = tasks.pop(task)
                try:
                    result = task.result()
                    connection.execute("INSERT OR REPLACE INTO receipts VALUES (?,?,?)", result)
                    completed += 1
                    consecutive = 0
                except Exception as exc:
                    errors.append({"db_id": item["dbId"], "error": str(exc)})
                    consecutive += 1
                    print(f"download failed db_id={item['dbId']}: {exc}", flush=True)
                if completed % 100 == 0 or errors:
                    connection.commit()
                    write_json(root / "download_status.json", {"expected": len(records), "downloaded": completed,
                               "complete": False, "stage": "downloading", "failed_this_run": len(errors),
                               "reactome_release": release})
                if consecutive >= max_failures:
                    stopped = True
            print(f"download cached={completed}/{len(records)} failures={len(errors)}", flush=True)
    except KeyboardInterrupt:
        stopped = True
        for task in tasks:
            task.cancel()
        print("Download interrupted; completed raw files will be reused", flush=True)
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
        connection.commit()
    # Files written by in-flight workers after interruption are validated on resume.
    status = {"expected": len(records), "downloaded": completed, "failed_this_run": len(errors),
              "complete": False, "stopped_after_failures_or_interrupt": stopped,
              "reactome_release": release, "stage": "checking_release_after_download"}
    write_jsonl(root / "download_errors.jsonl", sorted(errors, key=lambda row: row["db_id"]))
    write_json(root / "download_status.json", status)
    try:
        version_after = client.version()
    except (RuntimeError, OSError, ValueError) as exc:
        status.update(stage="stopped", error=str(exc), release_check_complete=False)
        write_json(root / "download_status.json", status)
        connection.close()
        raise
    if version_after != release:
        write_json(root / "release_mismatch.json", {"expected": release, "observed": version_after})
        status.update(stage="stopped", release_check_complete=False, error="Reactome release changed")
        write_json(root / "download_status.json", status)
        connection.close()
        raise ValueError("Release changed during downloads; this snapshot requires review")
    status = {"expected": len(records), "downloaded": completed, "failed_this_run": len(errors),
              "complete": completed == len(records), "stopped_after_failures_or_interrupt": stopped,
              "reactome_release": release}
    write_jsonl(root / "download_errors.jsonl", sorted(errors, key=lambda row: row["db_id"]))
    write_json(root / "download_status.json", status)
    connection.close()
    return status
