from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from backend.ml.torch_runtime import torch_runtime_available


REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_ROOT = REPO_ROOT / "backend"

for candidate in (REPO_ROOT, BACKEND_ROOT):
    candidate_str = str(candidate)
    if candidate_str not in sys.path:
        sys.path.insert(0, candidate_str)


def _module_available(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def _safe_torch_runtime_available(*extra_imports: str) -> bool:
    available, _ = torch_runtime_available(*extra_imports)
    return available


def pytest_ignore_collect(collection_path: Path, config: pytest.Config) -> bool:
    filename = collection_path.name

    if filename == "test_studio_dashboard.py" and not _module_available("sqlmodel"):
        return True

    if filename == "test_attention_integration.py":
        if not _module_available("torch") or not _module_available("torch_geometric"):
            return True
        if not _safe_torch_runtime_available("torch_geometric.data"):
            return True

    return False
