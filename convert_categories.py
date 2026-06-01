#!/usr/bin/env python3
"""
Convert common categories file formats to the standard YAML format
expected by bearing_hic_plot.py and batch_bearing_hic_plots.py.

Usage:
    python convert_categories.py input.json output.yaml
    python convert_categories.py input.yaml output.yaml
"""

import sys
import json
from pathlib import Path


def load_input_file(path):
    """Try to load JSON or YAML, auto-detect format."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    
    try:
        import yaml
        has_yaml = True
    except ImportError:
        has_yaml = False
        print("[WARN] PyYAML not available; can only parse JSON", file=sys.stderr)
    
    with open(path) as fh:
        content = fh.read()
    
    # Try JSON first
    try:
        return json.loads(content), "json"
    except json.JSONDecodeError:
        pass
    
    # Try YAML
    if has_yaml:
        try:
            import yaml
            return yaml.safe_load(content), "yaml"
        except Exception:
            pass
    
    raise ValueError(f"Could not parse {path} as JSON or YAML")


def convert_to_standard(data):
    """
    Convert various formats to standard format:
      categories:
        - name: "Name1"
          color: "#color1"
        - name: "Name2"
          color: "#color2"
    
    Handles:
    - List of dicts (already standard)
    - List of strings (adds default grey color)
    - Numeric-key dict with [name, color] pairs (bearmon format)
    - Simple name->color dict
    """
    if isinstance(data, dict) and "categories" in data:
        categories = data["categories"]
        if isinstance(categories, list) and len(categories) > 0:
            # Check if it's already the right format (list of dicts)
            if isinstance(categories[0], dict) and "name" in categories[0]:
                return data  # Already correct
            # Check if it's a list of strings
            elif isinstance(categories[0], str):
                # Convert ["name1", "name2"] -> [{"name": "name1", "color": "#grey"}, ...]
                print(
                    "[WARN] Input has simple string categories; adding default grey color",
                    file=sys.stderr,
                )
                return {
                    "categories": [
                        {"name": name, "color": "#cccccc"}
                        for name in categories
                    ]
                }
        # Check if categories is a dict (numeric keys or simple mapping)
        elif isinstance(categories, dict):
            # Try to detect format
            has_numeric_keys = all(k.isdigit() for k in categories.keys())
            if has_numeric_keys:
                # Numeric-key format: {"1": ["ATAC", "#be92e0"], "2": [...]}
                print(
                    "[INFO] Converting numeric-key categories to standard list format",
                    file=sys.stderr,
                )
                sorted_keys = sorted(categories.keys(), key=int)
                converted_cats = []
                for key in sorted_keys:
                    entry = categories[key]
                    if isinstance(entry, (list, tuple)) and len(entry) >= 2:
                        # Format: ["ATAC", "#be92e0"]
                        converted_cats.append({"name": entry[0], "color": entry[1]})
                    elif isinstance(entry, dict):
                        # Format: {"name": "ATAC", "color": "#be92e0"}
                        converted_cats.append(entry)
                    else:
                        raise ValueError(f"Unexpected category entry for key {key}: {entry}")
                return {"categories": converted_cats}
            else:
                # Simple mapping: {"ATAC": "#color1", ...}
                print(
                    "[INFO] Converting name->color mapping to standard list format",
                    file=sys.stderr,
                )
                return {
                    "categories": [
                        {"name": name, "color": color}
                        for name, color in categories.items()
                    ]
                }
    
    # Check if input is a simple dict mapping names to colors (no "categories" key)
    if isinstance(data, dict) and all(isinstance(v, str) for v in data.values()):
        print(
            "[INFO] Converting simple name->color mapping to standard format",
            file=sys.stderr,
        )
        return {
            "categories": [
                {"name": name, "color": color}
                for name, color in data.items()
            ]
        }
    
    raise ValueError(
        f"Could not auto-convert categories format. "
        f"Input structure: {type(data).__name__}"
    )


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    try:
        print(f"Loading {input_path}...", file=sys.stderr)
        data, fmt = load_input_file(input_path)
        print(f"  Detected format: {fmt}", file=sys.stderr)
        
        print(f"Converting to standard format...", file=sys.stderr)
        converted = convert_to_standard(data)
        
        print(f"Writing to {output_path}...", file=sys.stderr)
        try:
            import yaml
            with open(output_path, "w") as fh:
                yaml.dump(converted, fh, default_flow_style=False, sort_keys=False)
            print(f"[OK] Converted and saved to {output_path}", file=sys.stderr)
        except ImportError:
            print(
                "[ERROR] PyYAML required to write YAML. Install: pip install pyyaml",
                file=sys.stderr,
            )
            sys.exit(1)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
