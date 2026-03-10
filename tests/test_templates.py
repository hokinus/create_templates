#!/usr/bin/env python3
"""
Pytest test suite for create_templates_csv.py

This test discovers test cases from the results directory and validates
that get_template_labels produces the expected outputs.

Test case structure:
- {ID}.csv: input sequences file
- {ID}.txt: MMseqs alignment results file
- {ID}_results.csv: expected output (C1' coordinates only)
- {ID}_results.csv.allatom.csv: expected output (all atoms)
"""

from io import StringIO
import pytest
import pandas as pd
import pandas.testing as pd_testing
import numpy as np
import os
import sys
from pathlib import Path

# Add parent directory to path to import create_templates_csv
SCRIPT_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from create_templates_csv import get_template_labels_serial

# Import shared test utilities and discovery functions
from conftest import discover_test_cases, TemplateTestCase


def test_discovery(test_cases):
    """Test that test cases are discovered correctly."""
    assert len(test_cases) > 0, "No test cases discovered in results directory"
    test_ids = [tc.test_id for tc in test_cases]
    print(f"Discovered test cases: {test_ids}")


@pytest.mark.parametrize(
    "test_case",
    discover_test_cases(str(Path(__file__).parent / "results")),
    ids=lambda tc: tc.test_id,
)
def test_get_template_labels(test_case):
    """
    Test get_template_labels_serial against expected output.

    For each discovered test case, this test:
    1. Runs get_template_labels_serial with the input files
    2. Compares the generated output with expected output
    3. Validates both C1' coordinates and all-atom coordinates
    """

    # Determine CIF directory - should be in the project root
    test_dir = Path(test_case.sequences_file).parent
    project_root = test_dir.parent.parent  # tests/results -> tests -> project_root
    cif_dir = str(project_root / "PDB_RNA")

    # Get expected outputs
    expected_labels = pd.read_csv(test_case.expected_results)
    expected_allatom = pd.read_csv(test_case.expected_allatom)

    # Run the function
    try:
        output_labels, output_allatom_labels, targets, start_idx, end_idx = (
            get_template_labels_serial(
                sequences_file=test_case.sequences_file,
                mmseqs_results_file=test_case.mmseqs_file,
                skip_temporal_cutoff=True,
                MAX_TEMPLATES=40,
                cif_dir=cif_dir,
                id_map_file="",
                start_idx=0,
                end_idx=0,
            )
        )
    except Exception as e:
        pytest.fail(f"get_template_labels raised exception: {str(e)}")

    # Convert to DataFrame for comparison
    generated_labels = pd.DataFrame(output_labels)
    generated_allatom = pd.DataFrame(output_allatom_labels)

    # Convert to CSV for roundtrip comparison (deals with some weird '', ' ' string serializations and reading as NaN)
    tmp_io = StringIO()
    generated_labels.to_csv(tmp_io, index=False, float_format="%.3f")
    tmp_io.seek(0)

    generated_labels = pd.read_csv(tmp_io)
    tmp_io = StringIO()
    generated_allatom.to_csv(tmp_io, index=False, float_format="%.3f")
    tmp_io.seek(0)
    generated_allatom = pd.read_csv(tmp_io)

    # Compare C1' coordinates output
    pd_testing.assert_frame_equal(generated_labels, expected_labels, check_dtype=False)
    #
    # Compare all-atom coordinates output
    pd_testing.assert_frame_equal(
        generated_allatom, expected_allatom, check_dtype=False
    )


# ---------------------------------------------------------------------------
# Temporal-cutoff tests
#
# Uses the R1108 test case whose only template is 7qr3_C (released 2022-10-26).
# Per-target cutoff in the CSV is 2022-05-27 (before the release date).
# ---------------------------------------------------------------------------

CUTOFF_TEST_CASE_ID = "R1108"
TEMPLATE_RELEASE_DATE = "2022-10-26"  # release date of 7qr3_C
CUTOFF_BEFORE_RELEASE = "2022-05-01"  # before release → template should be excluded
CUTOFF_AFTER_RELEASE = "2023-01-01"  # after release  → template should be included


def _get_r1108(test_cases):
    """Return the R1108 TemplateTestCase, skipping if not present."""
    matches = [tc for tc in test_cases if tc.test_id == CUTOFF_TEST_CASE_ID]
    if not matches:
        pytest.skip(f"Test case {CUTOFF_TEST_CASE_ID} not found in results directory")
    return matches[0]


def _run_serial(tc, *, skip_temporal_cutoff=False, fixed_temporal_cutoff=""):
    """Helper: run get_template_labels_serial and return the labels DataFrame."""
    project_root = Path(tc.sequences_file).parent.parent.parent
    cif_dir = str(project_root / "PDB_RNA")
    output_labels, _, _, _, _ = get_template_labels_serial(
        sequences_file=tc.sequences_file,
        mmseqs_results_file=tc.mmseqs_file,
        skip_temporal_cutoff=skip_temporal_cutoff,
        MAX_TEMPLATES=40,
        cif_dir=cif_dir,
        fixed_temporal_cutoff=fixed_temporal_cutoff,
    )
    return pd.DataFrame(output_labels)


def test_temporal_cutoff_excludes_templates_before_release(test_cases):
    """With a cutoff before the template release date all coordinates must be NaN."""
    tc = _get_r1108(test_cases)
    df = _run_serial(tc, fixed_temporal_cutoff=CUTOFF_BEFORE_RELEASE)

    coord_cols = [c for c in df.columns if c.startswith(("x_", "y_", "z_"))]
    assert coord_cols, "Expected coordinate columns in output"
    assert df[coord_cols].isna().all().all(), (
        f"Expected all coordinates to be NaN when cutoff {CUTOFF_BEFORE_RELEASE} "
        f"is before template release {TEMPLATE_RELEASE_DATE}, but found non-NaN values"
    )


def test_temporal_cutoff_includes_templates_after_release(test_cases):
    """With a cutoff after the template release date coordinates must not all be NaN."""
    tc = _get_r1108(test_cases)
    df = _run_serial(tc, fixed_temporal_cutoff=CUTOFF_AFTER_RELEASE)

    coord_cols = [c for c in df.columns if c.startswith(("x_", "y_", "z_"))]
    assert coord_cols, "Expected coordinate columns in output"
    assert not df[coord_cols].isna().all().all(), (
        f"Expected at least some non-NaN coordinates when cutoff {CUTOFF_AFTER_RELEASE} "
        f"is after template release {TEMPLATE_RELEASE_DATE}"
    )


def test_fixed_temporal_cutoff_overrides_csv_cutoff(test_cases):
    """fixed_temporal_cutoff must override the per-target cutoff from the CSV.

    R1108's CSV cutoff (2022-05-27) is before the template release date, so without
    override the template would be excluded.  Passing a fixed cutoff *after* the
    release date should force the template to be included despite the CSV value.
    """
    tc = _get_r1108(test_cases)

    # Without override: CSV cutoff 2022-05-27 < release 2022-10-26 → excluded
    df_excluded = _run_serial(tc)
    coord_cols = [c for c in df_excluded.columns if c.startswith(("x_", "y_", "z_"))]
    assert df_excluded[coord_cols].isna().all().all(), (
        "Baseline: expected all NaN when using per-CSV cutoff that predates release"
    )

    # With fixed override after release → included
    df_included = _run_serial(tc, fixed_temporal_cutoff=CUTOFF_AFTER_RELEASE)
    assert not df_included[coord_cols].isna().all().all(), (
        "fixed_temporal_cutoff should override the CSV cutoff and include the template"
    )


def test_skip_temporal_cutoff_ignores_date(test_cases):
    """skip_temporal_cutoff=True must include the template regardless of dates."""
    tc = _get_r1108(test_cases)
    df = _run_serial(tc, skip_temporal_cutoff=True)

    coord_cols = [c for c in df.columns if c.startswith(("x_", "y_", "z_"))]
    assert not df[coord_cols].isna().all().all(), (
        "skip_temporal_cutoff=True should include all templates regardless of dates"
    )


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])
