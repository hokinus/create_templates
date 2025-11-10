# Template Tests

This directory contains pytest tests for the `create_templates_csv.py` module.

## Test Structure

### Test Files

- `test_templates.py` - Main test suite for template generation

### Test Data Structure

Each test case requires the following files in the `results/` directory:

```
results/
├── {ID}.csv                      # Input sequences file
├── {ID}.txt                      # MMseqs alignment results
├── {ID}_results.csv              # Expected output (C1' coordinates)
└── {ID}_results.csv.allatom.csv  # Expected output (all-atom coordinates)
```

#### File Formats

**{ID}.csv** - Input sequences

- Format: CSV with columns `target_id`, `sequence`, `temporal_cutoff`, etc.
- Contains the RNA sequences to be processed

**{ID}.txt** - MMseqs alignment results

- Format: Space-separated values (one alignment per line)
- Columns: `query`, `template`, `evalue`, `qstart`, `qend`, `tstart`, `tend`, `qaln`, `taln`

**{ID}_results.csv** - Expected C1' coordinate output

- Format: CSV with columns:
  - `ID` - Unique identifier (`{target_id}_{resid}`)
  - `resname` - Residue name (A, C, G, U, or X)
  - `resid` - Residue ID (1-based position in sequence)
  - `x_N`, `y_N`, `z_N` - C1' coordinates for template N (N=1 to 40)

**{ID}_results.csv.allatom.csv** - Expected all-atom coordinate output

- Format: CSV with columns:
  - `ID`, `resname`, `resid`
  - For each template N (N=1 to 40):
    - 26 atoms × 3 coordinates: `{ATOM}_x_N`, `{ATOM}_y_N`, `{ATOM}_z_N`
    - Atoms: P, OP1, OP2, O5', O3', C1', C2', C3', C4', O4', C5', N1, C2, O2, N3, C4, N4, C5, C6, O4, N9, N7, C8, N6, N2, O6
    - PDB info: `pdb_id_N`, `pdb_seq_num_N`, `pdb_ins_code_N`, `pdb_resname_N`

## Running Tests

### Run all tests

```bash
pixi run -e dev python -m pytest tests/test_templates.py -v
```

### Run a specific test

```bash
pixi run -e dev python -m pytest tests/test_templates.py::test_get_template_labels[R1108] -v
```

### Run with detailed output

```bash
pixi run -e dev python -m pytest tests/test_templates.py -v -s
```

### Run with coverage

```bash
pixi run -e dev python -m pytest tests/test_templates.py --cov=create_templates_csv
```

## Test Discovery

Tests are automatically discovered from the `results/` directory by:

1. Finding all `{ID}.txt` files
2. Checking that corresponding `.csv`, `_results.csv`, and `_results.csv.allatom.csv` files exist
3. Creating parameterized tests for each valid test case

Current test cases:

- R1108 - simple match
- R1116 - gapped (template and query) match
- R1126 - partial match, on template and query
- R1149 - partial match, partial match to query, multiple templates

## Test Implementation

The `test_get_template_labels` test:

1. Discovers test cases from the results directory
2. For each test case:
   - Loads input sequences and MMseqs results
   - Calls `get_template_labels()` with the test data
   - Compares generated output with expected results using `pandas.testing.assert_frame_equal()`
   - Validates both C1' coordinates and all-atom coordinates

## Environment Setup

The tests require a development environment with pytest installed. This is configured in `pixi.toml`:

```toml
[feature.dev.dependencies]
pytest = ">=7.0,<9"

[environments]
default = {features = []}
dev = {features = ["dev"]}
```

Install the dev environment:

```bash
pixi install
pixi run -e dev python -m pytest tests/test_templates.py
```

## Running Benchmarks

To benchmark the `get_template_labels` function performance:

```bash
pixi run -e dev pytest tests/benchmark_templates.py --benchmark-only
```

For more detailed output with statistics:

```bash
pixi run -e dev pytest tests/benchmark_templates.py --benchmark-only --benchmark-verbose
```

To save benchmark results automatically:

```bash
pixi run -e dev pytest tests/benchmark_templates.py --benchmark-only --benchmark-autosave
```

The benchmark suite automatically discovers test cases from `results/` and runs the template label generation on each case, measuring execution time and providing detailed performance statistics. Benchmarks can be compared across runs by using the `--benchmark-autosave` option.

## Adding New Tests

To add a new test case:

1. Create `{NEW_ID}.csv` with your test sequences
2. Create `{NEW_ID}.txt` with the MMseqs alignment results
3. Generate the expected output files:
   - Run `create_templates_csv.py` with your test data
   - Save the output as `{NEW_ID}_results.csv` and `{NEW_ID}_results.csv.allatom.csv`
4. The test will be automatically discovered and run

Example:

```bash
cd tests/results
python ../../create_templates_csv.py \
  -s NEW_ID.csv \
  --mmseqs_results_file NEW_ID.txt \
  --cif_dir ../../PDB_RNA \
  -o ./ \
  --outfile NEW_ID_results.csv
```
