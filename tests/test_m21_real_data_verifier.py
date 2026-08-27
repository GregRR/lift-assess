from __future__ import annotations

import runpy


def test_m21_real_data_verifier_loads_recorded_cases() -> None:
    # Import/load only: main() requires the external canFam3->canFam4 fixture cache.
    namespace = runpy.run_path("scripts/verify_m21_b12_b14.py")

    assert "main" in namespace
    cases = namespace["CASES"]
    assert tuple(case.label for case in cases) == ("B12", "B13", "B14")
    assert tuple(case.locus for case in cases) == (
        "chrX:26956239-26956239",
        "chr28:1484906-1484906",
        "chr5:31705136-31705136",
    )
    assert cases[2].all_chain_count == 9
    assert cases[2].filtered_chain_count == 1
    assert cases[2].relationship == "FAVORS_ONE_PLACEMENT"
