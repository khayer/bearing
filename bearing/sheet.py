from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from .tsv import read_tsv_table


def load_rows(path: Path) -> List[Dict[str, str]]:
    """Load normalized sheet rows from TSV, preserving all columns."""
    _, rows = read_tsv_table(path)
    return rows


def require_columns(fields: List[str], required: List[str], context: str) -> None:
    """Raise ValueError when one or more required columns are missing."""
    missing = [c for c in required if c not in fields]
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
    values: List[str] = []
    if row.get("bw"):
        values.extend([p.strip() for p in row["bw"].split(",") if p.strip()])
        return values

    i = 1
    while True:
        key = f"bw{i}"
        if key not in row:
            break
        value = row.get(key, "").strip()
        if value:
            values.append(value)
        i += 1
    return values
