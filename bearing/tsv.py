from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .utils import iter_non_comment_lines


def _normalize_tsv_cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ",".join(str(item).strip() for item in value if item is not None)
    return str(value).strip()


def read_tsv_table(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """Read a TSV and return (normalized_fields, normalized_rows)."""
    rows: List[Dict[str, str]] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(iter_non_comment_lines(fh), delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("TSV is empty or missing header")

        fields = [f.strip().lower() for f in reader.fieldnames]
        for raw in reader:
            row = {
                (k.strip().lower() if k is not None else ""): _normalize_tsv_cell(v)
                for k, v in raw.items()
            }
            rows.append(row)
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
        missing = [c for c in required_columns if c not in fields]
        if missing:
            raise ValueError(f"TSV missing required columns: {', '.join(missing)}")
    return rows
