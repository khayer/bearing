from __future__ import annotations

from pathlib import Path
from typing import Iterable


def ensure_paths_exist(paths: Iterable[Path], context: str = "input") -> None:
    """Raise FileNotFoundError for the first missing path."""
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"missing {context}: {path}")
