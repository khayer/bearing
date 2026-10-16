from __future__ import annotations

import concurrent.futures
import shlex
import subprocess
from typing import List, Sequence, Tuple


def format_command(cmd: Sequence[str]) -> str:
    """Render ``cmd`` as a shell-safe, copy-pasteable string."""
    return " ".join(shlex.quote(str(part)) for part in cmd)


def emit_command(cmd: Sequence[str], label: str = "") -> None:
    """Print an optional ``[label]`` header followed by the indented command."""
    if label:
        print(f"[{label}]")
    print(f"  $ {format_command(cmd)}")


def _run(cmd: Sequence[str]) -> int:
    return subprocess.run(list(cmd)).returncode


def execute_command(cmd: Sequence[str], run: bool, label: str = "") -> int:
    """Print ``cmd`` and, if ``run``, execute it. Returns its exit code (0 on dry run)."""
    emit_command(cmd, label=label)
    return _run(cmd) if run else 0


def execute_commands(commands: List[Sequence[str]], run: bool, threads: int = 1) -> Tuple[int, int]:
    """Execute command list. Returns (total, failures).

    With ``run=False`` commands are only printed. With ``threads > 1`` all
    commands are announced as queued, then run concurrently.
    """
    total = len(commands)
    if not run:
        for cmd in commands:
            emit_command(cmd)
        return total, 0

    status = "running" if threads <= 1 else "queued"

    def announce(i: int, cmd: Sequence[str]) -> None:
        print(f"[{i}/{total}] {status}")
        emit_command(cmd)

    if threads <= 1:
        return_codes = []
        for i, cmd in enumerate(commands, start=1):
            announce(i, cmd)
            return_codes.append(_run(cmd))
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
            for i, cmd in enumerate(commands, start=1):
                announce(i, cmd)
            return_codes = list(ex.map(_run, commands))

    return total, sum(1 for rc in return_codes if rc != 0)
