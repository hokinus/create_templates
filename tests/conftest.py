#!/usr/bin/env python3
"""
Shared pytest configuration and utilities for test and benchmark suites.

Provides common test case discovery and fixtures used by multiple test modules.
"""

import pytest
from pathlib import Path
from collections import namedtuple

# Named tuple for test case data
TemplateTestCase = namedtuple('TemplateTestCase', ['test_id', 'sequences_file', 'mmseqs_file', 'expected_results', 'expected_allatom'])


def discover_test_cases(results_dir: str) -> list:
    """
    Discover test cases from the results directory.
    
    Returns list of TemplateTestCase namedtuples containing paths to:
    - test input files ({ID}.csv, {ID}.txt)
    - expected output files ({ID}_results.csv, {ID}_results.csv.allatom.csv)
    
    Args:
        results_dir: Path to the results directory containing test files
        
    Returns:
        List of TemplateTestCase namedtuples sorted by test_id
    """
    results_path = Path(results_dir)
    test_cases = []
    
    # Find all {ID}.txt files to identify test IDs
    for txt_file in results_path.glob("*.txt"):
        
        test_id = txt_file.stem  # filename without extension
        
        # Check for corresponding input files
        input_file = results_path / f"{test_id}.csv"
        mmseqs_file = txt_file  # The .txt file is the MMseqs results file
        expected_results = results_path / f"{test_id}_results.csv"
        expected_allatom = results_path / f"{test_id}_results.csv.allatom.csv"
        
        # Only add if all required files exist
        if input_file.exists() and expected_results.exists() and expected_allatom.exists():
            test_cases.append(TemplateTestCase(
                test_id=test_id,
                sequences_file=str(input_file),
                mmseqs_file=str(mmseqs_file),
                expected_results=str(expected_results),
                expected_allatom=str(expected_allatom)
            ))
    
    return sorted(test_cases, key=lambda x: x.test_id)


@pytest.fixture(scope="session")
def test_cases():
    """Discover and return all test cases."""
    results_dir = Path(__file__).parent / "results"
    cases = discover_test_cases(str(results_dir))
    return cases
