from __future__ import annotations

import runpy


def test_mechanical_fixture_verifier_loads_current_public_api() -> None:
    # Import/load only: main() requires the multi-gigabyte external fixture cache.
    namespace = runpy.run_path("scripts/verify_canFam3_canFam4_mechanical_fixture.py")

    assert "main" in namespace
