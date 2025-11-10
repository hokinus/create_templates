#!/usr/bin/env python3
"""
Pytest-benchmark suite for create_templates_csv.py

Benchmarks the get_template_labels function using test cases
discovered from the results directory.
"""

import pytest
import pandas as pd
import os
import sys
from pathlib import Path

# Add parent directory to path to import create_templates_csv
SCRIPT_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from create_templates_csv import get_template_labels

# Import shared test utilities and discovery functions
from conftest import discover_test_cases, TemplateTestCase


@pytest.mark.parametrize("test_case", 
                         discover_test_cases(str(Path(__file__).parent / "results")),
                         ids=lambda tc: tc.test_id)
def test_benchmark_get_template_labels(benchmark, test_case):
    """
    Benchmark get_template_labels against discovered test cases.
    
    For each discovered test case, this benchmark:
    1. Measures the execution time of get_template_labels
    2. Reports timing statistics
    """
    
    # Determine CIF directory - should be in the project root
    test_dir = Path(test_case.sequences_file).parent
    project_root = test_dir.parent.parent  # tests/results -> tests -> project_root
    cif_dir = str(project_root / "PDB_RNA")
    
    # Define the function to benchmark
    def run_get_template_labels():
        return get_template_labels(
            sequences_file=test_case.sequences_file,
            mmseqs_results_file=test_case.mmseqs_file,
            skip_temporal_cutoff=True,
            MAX_TEMPLATES=40,
            cif_dir=cif_dir,
            id_map_file='',
            start_idx=0,
            end_idx=0
        )
    
    # Run the benchmark
    result = benchmark(run_get_template_labels)
    
    # Verify that function returns expected number of outputs
    assert len(result) == 3, "get_template_labels should return 3 items (output_labels, output_allatom_labels, targets)"


if __name__ == "__main__":
    # Run benchmarks with verbose output
    pytest.main([__file__, "-v", "--benchmark-only", "--benchmark-autosave"])
