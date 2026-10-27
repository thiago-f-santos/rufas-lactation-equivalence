"""Transactional atomic file writing utilities."""

import json
import os
import tempfile
from pathlib import Path
from typing import Any


def atomic_write_text(target_path: Path, content: str) -> None:
    """
    Writes text to target_path atomically using a temporary file in the same directory,
    os.fsync(), and os.replace().
    """
    target_path = Path(target_path)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        dir=str(target_path.parent),
        prefix=f".{target_path.name}.",
        suffix=".tmp",
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(content)
            f.flush()
            os.fsync(f.fileno())
        os.replace(str(tmp_path), str(target_path))
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)
        raise


def atomic_write_json(target_path: Path, payload: dict[str, Any]) -> None:
    """Serializes payload to JSON in memory first, then writes atomically to target_path."""
    serialized = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    atomic_write_text(target_path, serialized)
