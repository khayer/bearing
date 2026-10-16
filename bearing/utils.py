from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple


def _split_key_value(item: str, missing_sep_msg: str, empty_msg: str) -> Tuple[str, str]:
    """Split ``KEY=VALUE`` on the first '=', stripping both sides.

    Raises ValueError with ``missing_sep_msg`` when there is no '=' and with
    ``empty_msg`` when either side is empty after stripping.
    """
    if "=" not in item:
        raise ValueError(missing_sep_msg)
    key, value = (part.strip() for part in item.split("=", 1))
    if not key or not value:
        raise ValueError(empty_msg)
    return key, value


def parse_key_value_items(items: List[str], arg_name: str) -> Dict[str, str]:
    """Parse repeatable CLI items of form KEY=VALUE into a dictionary.

    Later duplicates of a key overwrite earlier ones.
    """
    out: Dict[str, str] = {}
    for item in items:
        key, value = _split_key_value(
            item,
            missing_sep_msg=f"{arg_name} expects CONDITION=VALUE, got: {item}",
            empty_msg=f"Invalid {arg_name} item: {item}",
        )
        out[key] = value
    return out


def sanitize_token(text: str, allowed: str = "._-") -> str:
    """Return a filesystem-safe token, replacing disallowed characters with '_'"""
    return "".join(ch if ch.isalnum() or ch in allowed else "_" for ch in text)


def resolve_path(path_value: str, base_dir: Path) -> Path:
    """Resolve absolute or base-dir-relative paths to absolute paths."""
    path = Path(path_value)
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def iter_non_comment_lines(handle: Iterable[str]) -> Iterator[str]:
    """Yield non-empty lines that do not start with '#' after stripping."""
    for line in handle:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        yield line


def parse_ucsc_region(region_str: str) -> Tuple[str, int, int]:
    """Parse regions like 'chr6:50000000-52000000' into (chrom, start, end)."""
    chrom, coords = region_str.split(":", 1)
    start_s, end_s = coords.replace(",", "").split("-", 1)
    return chrom, int(start_s), int(end_s)


def parse_named_region_item(item: str, arg_name: str = "region") -> Dict[str, str]:
    """Parse a CLI item in NAME=chr:start-end format."""
    name, region = _split_key_value(
        item,
        missing_sep_msg=f"{arg_name} must be NAME=chr:start-end, got: {item}",
        empty_msg=f"Invalid {arg_name} entry: {item}",
    )
    return {"name": name, "region": region}


def missing_columns(fields: Iterable[str], required: Iterable[str]) -> List[str]:
    """Return the ``required`` columns absent from ``fields``, in required order."""
    present = set(fields)
    return [c for c in required if c not in present]
