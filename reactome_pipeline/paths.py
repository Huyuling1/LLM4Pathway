"""Checkpointed exhaustive enumeration of directed maximal simple paths."""

import json
import os
import signal
import sqlite3
import threading
import time
from pathlib import Path

from .common import digest, dumps, write_json, write_jsonl
from .graph import components, load_graph, weak_components
from .render import chain_text
from .coverage import COVERAGE_VERSION, input_coverage

ALGORITHM_VERSION = "maximal-simple-dfs-v1"
RECORD_VERSION = "named-reactions-cycle-connections-v4.3"


def start_nodes(adjacency, reverse):
    """A maximal simple path starts at a source or at a vertex in a cycle."""
    starts = {node for node in adjacency if not reverse[node]}
    for group in components(adjacency, reverse):
        if len(group) > 1 or group[0] in adjacency[group[0]]:
            starts.update(group)
    return sorted(starts)


def enumerate_chains(root, stop_after=None, checkpoint_every=1000, max_steps=None,
                     max_seconds=None, small_components_first=False):
    if stop_after is not None and stop_after < 1:
        raise ValueError("stop-after must be positive")
    if checkpoint_every < 1:
        raise ValueError("checkpoint-every must be positive")
    if (max_steps is not None and max_steps < 1) or (max_seconds is not None and max_seconds <= 0):
        raise ValueError("Enumeration work budgets must be positive")
    root = Path(root)
    graph_db, reactions, adjacency, reverse, metadata = load_graph(root)
    out = root / "paths"
    out.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(out / "paths.sqlite")
    db.executescript("""
        CREATE TABLE IF NOT EXISTS progress (id INTEGER PRIMARY KEY CHECK(id=1), state TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS chains (ordinal INTEGER PRIMARY KEY, id TEXT UNIQUE, status TEXT, data TEXT);
    """)
    fingerprint_fields = {"graph": metadata["graph_sha256"], "algorithm": ALGORITHM_VERSION,
                          "record_version": RECORD_VERSION, "coverage_version": COVERAGE_VERSION}
    if small_components_first:
        fingerprint_fields["schedule"] = "weak-component-size-v1"
    fingerprint = digest(fingerprint_fields)
    starts = start_nodes(adjacency, reverse)
    if small_components_first:
        rank = {node: i for i, group in enumerate(weak_components(adjacency, reverse)) for node in group}
        starts.sort(key=lambda node: (rank[node], node))
    saved = db.execute("SELECT state FROM progress WHERE id=1").fetchone()
    state = json.loads(saved[0]) if saved else {"fingerprint": fingerprint, "start_index": 0, "path": [], "frames": [], "steps": 0, "leaves": 0, "complete": False}
    if state["fingerprint"] != fingerprint:
        db.close()
        graph_db.close()
        raise ValueError("Graph, algorithm or chain record version changed; archive the existing paths directory before re-enumerating")
    before = db.execute("SELECT count(*) FROM chains").fetchone()[0]
    last_commit = time.monotonic()
    started_at = last_commit
    initial_steps = state["steps"]
    stop_reason = None
    visited = set(state["path"])
    steps_since_commit = 0
    interrupted = False
    previous_handler = None
    def request_stop(signum, frame):
        nonlocal interrupted
        interrupted = True

    def checkpoint():
        nonlocal last_commit, steps_since_commit
        db.execute("INSERT OR REPLACE INTO progress VALUES (1,?)", (dumps(state),))
        db.commit()
        write_json(out / "progress.json", {"graph_sha256": metadata["graph_sha256"], "enumeration_complete": state["complete"],
                                          "starts": len(starts), "starts_entered": state["start_index"],
                                          "dfs_steps": state["steps"], "candidate_leaves": state["leaves"],
                                          "active_path": state["path"]})
        last_commit = time.monotonic()
        steps_since_commit = 0

    def emit():
        # A forward-maximal path can still be a suffix of another simple path.
        if any(parent not in visited for parent in reverse[state["path"][0]]):
            return False
        path = list(state["path"])
        last = path[-1]
        status = "cycle_blocked" if adjacency[last] else "sink"
        def connection(source, target):
            evidence = []
            for kind, entity, name in graph_db.execute("SELECT kind,entity,name FROM edges JOIN entities ON entity=id WHERE source=? AND target=? ORDER BY kind,entity", (source, target)):
                edge = {"source": source, "target": target, "kind": kind, "entity_id": entity}
                evidence.append({"edge_id": digest(edge), "kind": kind, "entity_id": entity, "entity_name": name})
            return {"source": source, "target": target, "evidence": evidence}
        transitions = [connection(source, target) for source, target in zip(path, path[1:])]
        coverage = input_coverage(path, reactions)
        closing_edges = [connection(last, child) for child in adjacency[last]] if status == "cycle_blocked" else []
        row = {"id": digest({"fingerprint": fingerprint, "reaction_ids": path}), "reaction_ids": path,
               "record_version": RECORD_VERSION, "input_coverage_version": COVERAGE_VERSION,
               "input_coverage": coverage,
               "status": status, "transitions": transitions,
               "closing_edges": closing_edges,
               "text": chain_text(path, reactions, transitions, closing_edges)}
        cursor = db.execute("INSERT OR IGNORE INTO chains(id,status,data) VALUES (?,?,?)", (row["id"], status, dumps(row)))
        return cursor.rowcount > 0

    emitted = 0
    try:
        if threading.current_thread() is threading.main_thread():
            previous_handler = signal.signal(signal.SIGINT, request_stop)
        while not state["complete"] and not interrupted:
            if not state["path"]:
                if state["start_index"] >= len(starts):
                    state["complete"] = True
                    break
                start = starts[state["start_index"]]
                state["start_index"] += 1
                state["path"] = [start]
                state["frames"] = [{"index": 0, "extended": False}]
                visited = {start}
            node = state["path"][-1]
            frame = state["frames"][-1]
            child = None
            while frame["index"] < len(adjacency[node]):
                candidate = adjacency[node][frame["index"]]
                frame["index"] += 1
                if candidate not in visited:
                    child = candidate
                    break
            if child is not None:
                frame["extended"] = True
                state["path"].append(child)
                state["frames"].append({"index": 0, "extended": False})
                visited.add(child)
            else:
                if not frame["extended"]:
                    state["leaves"] += 1
                    emitted += int(emit())
                visited.remove(state["path"].pop())
                state["frames"].pop()
            state["steps"] += 1
            steps_since_commit += 1
            if steps_since_commit >= checkpoint_every or time.monotonic() - last_commit >= 15:
                checkpoint()
            if stop_after is not None and emitted >= stop_after:
                stop_reason = "record_budget"
                break
            if max_steps is not None and state["steps"] - initial_steps >= max_steps:
                stop_reason = "dfs_step_budget"
                break
            if max_seconds is not None and time.monotonic() - started_at >= max_seconds:
                stop_reason = "time_budget"
                break
    except KeyboardInterrupt:
        interrupted = True
    finally:
        checkpoint()
        if previous_handler is not None:
            signal.signal(signal.SIGINT, previous_handler)
    counts = dict(db.execute("SELECT status,count(*) FROM chains GROUP BY status"))
    report = {"algorithm": ALGORITHM_VERSION, "graph_sha256": metadata["graph_sha256"],
              "record_version": RECORD_VERSION, "input_coverage_version": COVERAGE_VERSION,
              "input_complete": metadata["input_complete"], "enumeration_complete": state["complete"],
              "complete": state["complete"] and metadata["input_complete"],
              "starts": len(starts), "starts_entered": state["start_index"], "dfs_steps": state["steps"],
              "candidate_leaves": state["leaves"], "sink_chains": counts.get("sink", 0),
              "cycle_stopped_chains": counts.get("cycle_blocked", 0),
              "new_records": sum(counts.values()) - before, "interrupted": interrupted,
              "stop_reason": stop_reason, "start_order": "weak-component-size-v1" if small_components_first else "stable-id",
              "active_path": state["path"]}
    for status, filename in (("sink", "complete_chains.jsonl"), ("cycle_blocked", "cycle_stopped_chains.jsonl")):
        rows = (json.loads(row[0]) for row in db.execute("SELECT data FROM chains WHERE status=? ORDER BY ordinal", (status,)))
        write_jsonl(out / filename, rows)
    # Stream the readable global view; keep the two structured termination groups.
    text_path = root / "chains.txt"
    temporary = root / "chains.txt.tmp"
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        first = True
        for status in ("sink", "cycle_blocked"):
            for (data,) in db.execute("SELECT data FROM chains WHERE status=? ORDER BY ordinal", (status,)):
                if not first:
                    stream.write("\n")
                stream.write(json.loads(data)["text"] + "\n")
                first = False
    os.replace(temporary, text_path)
    write_json(out / "status.json", report)
    db.close()
    graph_db.close()
    return report
