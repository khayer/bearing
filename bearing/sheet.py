from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from .tsv import read_tsv_table
from .utils import missing_columns


def load_rows(path: Path) -> List[Dict[str, str]]:
    """Load normalized sheet rows from TSV, preserving all columns."""
    _, rows = read_tsv_table(path)
    return rows


def require_columns(fields: List[str], required: List[str], context: str) -> None:
    """Raise ValueError when one or more required columns are missing."""
    missing = missing_columns(fields, required)
    if missing:
        raise ValueError(f"{context} missing required columns: {', '.join(missing)}")


def require_any_column(fields: List[str], candidates: List[str], context: str) -> None:
    """Raise ValueError when none of the candidate columns are present."""
    if any(c in fields for c in candidates):
        return
    raise ValueError(
        f"{context} must include at least one of: {', '.join(candidates)}"
    )


def sample_name_from_row(row: Dict[str, str]) -> str:
    """Resolve sample name from 'sample' (preferred) or 'name'."""
    return row.get("sample", "") or row.get("name", "")


def extract_bw_values(row: Dict[str, str]) -> List[str]:
    """Extract BigWig path values from 'bw' or bw1,bw2,... style columns."""
    if row.get("bw"):
        return [p.strip() for p in row["bw"].split(",") if p.strip()]

    # Indexed columns: walk bw1, bw2, ... until the first absent key;
    # present-but-empty cells are skipped without stopping the walk.
    values: List[str] = []
    i = 1
    while f"bw{i}" in row:
        value = row[f"bw{i}"].strip()
        if value:
            values.append(value)
        i += 1
    return values
