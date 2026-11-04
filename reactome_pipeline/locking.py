"""One writer per snapshot; readers can inspect status while it runs."""

import fcntl
from pathlib import Path


class SnapshotLock:
    def __init__(self, root):
        root = Path(root)
        root.mkdir(parents=True, exist_ok=True)
        self.stream = (root / ".pipeline.lock").open("a")

    def __enter__(self):
        try:
            fcntl.flock(self.stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            self.stream.close()
            raise ValueError("Another processing command is writing this snapshot") from exc
        return self

    def __exit__(self, *args):
        fcntl.flock(self.stream, fcntl.LOCK_UN)
        self.stream.close()
