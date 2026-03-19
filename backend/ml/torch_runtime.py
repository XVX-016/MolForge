from __future__ import annotations

import subprocess
import sys
from functools import lru_cache
from typing import Iterable, Tuple


def _probe_code(imports: Tuple[str, ...]) -> str:
    lines = []
    for module_name in imports:
        lines.append(f"import {module_name}")
    lines.append("print('ok')")
    return "\n".join(lines)


@lru_cache(maxsize=None)
def probe_imports(imports: Tuple[str, ...]) -> tuple[bool, str]:
    """
    Probe optional imports in a subprocess so native DLL load failures do not
    crash the current Python process.
    """
    try:
        result = subprocess.run(
            [sys.executable, "-c", _probe_code(imports)],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
    except Exception as exc:  # pragma: no cover - defensive guard
        return False, str(exc)

    if result.returncode == 0:
        return True, ""

    stderr = (result.stderr or "").strip()
    stdout = (result.stdout or "").strip()
    message = stderr or stdout or f"subprocess exited with code {result.returncode}"
    return False, message


def torch_runtime_available(*extra_imports: str) -> tuple[bool, str]:
    imports = ("torch",) + tuple(extra_imports)
    return probe_imports(imports)
