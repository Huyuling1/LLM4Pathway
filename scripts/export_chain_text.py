"""Export saved chain text verbatim, in sink-then-cycle record order."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from reactome_pipeline.common import read_json, read_jsonl, write_jsonl


def export(root, out):
    root, out = Path(root), Path(out)
    status = read_json(root / "paths/status.json")
    rows = [row for name in ("complete_chains.jsonl", "cycle_stopped_chains.jsonl")
            for row in read_jsonl(root / "paths" / name)]
    out.mkdir(parents=True, exist_ok=True)
    write_jsonl(out / "chains.jsonl", rows)
    (out / "chains.txt").write_text("\n\n".join(row["text"] for row in rows) + ("\n" if rows else ""), encoding="utf-8")
    print({"records": len(rows), "record_version": status["record_version"],
           "input_complete": status["input_complete"], "enumeration_complete": status["enumeration_complete"]})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    export(args.root, args.out)
