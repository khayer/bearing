from __future__ import annotations

import concurrent.futures
import shlex
import subprocess
from typing import List, Sequence, Tuple


def format_command(cmd: Sequence[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in cmd)


def emit_command(cmd: Sequence[str], label: str = "") -> None:
    if label:
        print(f"[{label}]")
    print(f"  $ {format_command(cmd)}")


def execute_command(cmd: Sequence[str], run: bool, label: str = "") -> int:
    emit_command(cmd, label=label)
    if not run:
        return 0
    return subprocess.run(list(cmd)).returncode


def execute_commands(commands: List[Sequence[str]], run: bool, threads: int = 1) -> Tuple[int, int]:
    """Execute command list. Returns (total, failures)."""
    total = len(commands)
    if not run:
        for cmd in commands:
            emit_command(cmd)
        return total, 0

    if threads <= 1:
        failures = 0
        for i, cmd in enumerate(commands, start=1):
            print(f"[{i}/{total}] running")
            print(f"  $ {format_command(cmd)}")
            rc = subprocess.run(list(cmd)).returncode
            if rc != 0:
                failures += 1
        return total, failures

    def _run_one(cmd: Sequence[str]) -> int:
        return subprocess.run(list(cmd)).returncode

    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
        for i, cmd in enumerate(commands, start=1):
            print(f"[{i}/{total}] queued")
            print(f"  $ {format_command(cmd)}")
        futures = [ex.submit(_run_one, c) for c in commands]
        for fut in concurrent.futures.as_completed(futures):
            rc = fut.result()
            if rc != 0:
                failures += 1
    return total, failures
