"""Shared utilities for Bearing scripts."""

from .utils import (
    iter_non_comment_lines,
    parse_key_value_items,
    parse_named_region_item,
    parse_ucsc_region,
    resolve_path,
    sanitize_token,
)
from .tsv import read_tsv_dict_rows
from .tsv import read_tsv_table
from .sheet import (
    extract_bw_values,
    load_rows,
    require_any_column,
    require_columns,
    sample_name_from_row,
)
from .validate import ensure_paths_exist
from .diagnostics import require_files, require_mapping_paths
from .runner import execute_command, execute_commands, format_command

__all__ = [
    "iter_non_comment_lines",
    "parse_key_value_items",
    "parse_named_region_item",
    "parse_ucsc_region",
    "resolve_path",
    "sanitize_token",
    "read_tsv_dict_rows",
    "read_tsv_table",
    "load_rows",
    "require_columns",
    "require_any_column",
    "sample_name_from_row",
    "extract_bw_values",
    "ensure_paths_exist",
    "require_files",
    "require_mapping_paths",
    "format_command",
    "execute_command",
    "execute_commands",
]
