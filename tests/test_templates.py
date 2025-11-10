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

from create_templates_csv import get_template_labels

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
    Test get_template_labels against expected output.

    For each discovered test case, this test:
    1. Runs get_template_labels with the input files
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
        output_labels, output_allatom_labels, targets = get_template_labels(
            sequences_file=test_case.sequences_file,
            mmseqs_results_file=test_case.mmseqs_file,
            skip_temporal_cutoff=True,
            MAX_TEMPLATES=40,
            cif_dir=cif_dir,
            id_map_file="",
            start_idx=0,
            end_idx=0,
        )
    except Exception as e:
        pytest.fail(f"get_template_labels raised exception: {str(e)}")

    # Convert to DataFrame for comparison
    generated_labels = pd.DataFrame(output_labels)
    generated_allatom = pd.DataFrame(output_allatom_labels)

    # Convert to CSV for roundtrip comparison (deals with some weird '', ' ' string serializations and reading as NaN)
    tmp_io = StringIO()
    generated_labels.to_csv(tmp_io, index=False)
    tmp_io.seek(0)

    generated_labels = pd.read_csv(tmp_io)
    tmp_io = StringIO()
    generated_allatom.to_csv(tmp_io, index=False)
    tmp_io.seek(0)
    generated_allatom = pd.read_csv(tmp_io)

    # Compare C1' coordinates output
    pd_testing.assert_frame_equal(generated_labels, expected_labels, check_dtype=False)
    #
    # Compare all-atom coordinates output
    pd_testing.assert_frame_equal(
        generated_allatom, expected_allatom, check_dtype=False
    )


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])
