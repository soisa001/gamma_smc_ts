from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Return the repository root, which run2 needs to find tracked resources."""
    return Path(__file__).resolve().parent.parent
