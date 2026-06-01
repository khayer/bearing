from __future__ import annotations

from pathlib import Path
from typing import Dict

from .validate import ensure_paths_exist


def require_files(paths, context: str) -> None:
    ensure_paths_exist(paths, context=context)


def require_mapping_paths(mapping: Dict[str, str], resolver, context: str) -> None:
    resolved = [resolver(path_value) for path_value in mapping.values()]
    ensure_paths_exist(resolved, context=context)
