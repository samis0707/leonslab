"""Pytest configuration: makes lib/ importable and sets up data dir for tests."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "lib"))


@pytest.fixture
def data_dir() -> Path:
    """Path to the bundled data directory."""
    return ROOT / "data" / "primer_design"


@pytest.fixture
def fixtures_dir() -> Path:
    """Path to test fixtures (mini genomes, mini gene FASTAs for unit tests)."""
    return ROOT / "tests" / "fixtures"


@pytest.fixture
def golden_dir() -> Path:
    """Path to golden artefacts for regression tests."""
    return ROOT / "tests" / "golden"
