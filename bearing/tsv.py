from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .utils import iter_non_comment_lines, missing_columns


def _normalize_header(name: Optional[str]) -> str:
    """Lowercase/strip a header; csv uses ``None`` as the key for overflow cells."""
    return name.strip().lower() if name is not None else ""


def _normalize_tsv_cell(value) -> str:
    """Strip a cell; overflow cells (lists) are comma-joined, ``None`` becomes ''."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ",".join(str(item).strip() for item in value if item is not None)
    return str(value).strip()


def read_tsv_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """Read a TSV and return (normalized_fields, normalized_rows)."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(iter_non_comment_lines(fh), delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("TSV is empty or missing header")

        fields = [_normalize_header(f) for f in reader.fieldnames]
        rows = [
            {_normalize_header(k): _normalize_tsv_cell(v) for k, v in raw.items()}
            for raw in reader
        ]
    return fields, rows


def read_tsv_dict_rows(
    path: Path,
    required_columns: Optional[Iterable[str]] = None,
) -> List[Dict[str, str]]:
    """Read a TSV into normalized lowercase-key dictionaries.

    Comment lines beginning with '#' and blank lines are ignored.
    """
    fields, rows = read_tsv_table(path)
    if required_columns:
        missing = missing_columns(fields, required_columns)
        if missing:
            raise ValueError(f"TSV missing required columns: {', '.join(missing)}")
    return rows
