from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Tuple


def parse_key_value_items(items: List[str], arg_name: str) -> Dict[str, str]:
    """Parse repeatable CLI items of form KEY=VALUE into a dictionary."""
    out: Dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"{arg_name} expects CONDITION=VALUE, got: {item}")
        key, value = item.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key or not value:
            raise ValueError(f"Invalid {arg_name} item: {item}")
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
    if "=" not in item:
        raise ValueError(f"{arg_name} must be NAME=chr:start-end, got: {item}")
    name, region = item.split("=", 1)
    name = name.strip()
    region = region.strip()
    if not name or not region:
        raise ValueError(f"Invalid {arg_name} entry: {item}")
    return {"name": name, "region": region}
