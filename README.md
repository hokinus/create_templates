# create_templates
Create 3D templates (or solution) file for RNA.

Currently based on MMseqs to carry out search and alignment.

See also [MMseqs2 3D RNA Template notebook](https://www.kaggle.com/code/rhijudas/mmseqs2-3d-rna-template-identification) on Kaggle, which has the full workflow! 

Requirements:   
1. `biopython`  
2. `PDB_RNA/` directory holding  
 - cif.gz or .cif files (get ID's from Advanced Search in PDB, searching for polymer entity RNA and then [batch download](https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script)  
 - `pdb_release_dates_NA.csv` (get from Advanced Search in PDB), and  
 - `pdb_seqres_NA.fasta` (available at [link](https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz)).  
 
 
 Note that a minimal `PDB_RNA` folder that works for the example is provided here. A version used for the Kaggle RNA folding competition in May 2025 is available [here](https://www.kaggle.com/competitions/stanford-rna-3d-folding/data). If that is not accessible, check out this [clone](https://www.kaggle.com/datasets/rhijudas/clone-of-stanford-rna-3d-modeling-competition-data).


## Example command line

```
cd example/
python3 ../create_templates_csv.py \
	 -s validation_sequences.csv \
	--mmseqs_results_file validation_Result.txt \
	--skip_temporal_cutoff \
	--outfile validation_templates.csv 
```

Output should match what is in `example/example_output/validation_templates.csv`

## All options

```
% python3 create_templates_csv.py -h

usage: create_templates_csv.py [-h] [-s SEQUENCES_FILE] [--mmseqs_results_file MMSEQS_RESULTS_FILE] [--outfile OUTFILE] [--dataset_name DATASET_NAME] [-o OUTDIR] [--max_templates MAX_TEMPLATES] [--cif_dir CIF_DIR]
                               [--skip_temporal_cutoff] [--start_idx START_IDX] [--end_idx END_IDX] [--id_map ID_MAP]

Prepare templates.csv file similar to labels.csv but with MMseqs2-identified templates

options:
  -h, --help            show this help message and exit
  -s, --sequences_file SEQUENCES_FILE
                        CSV file with columns including "target_id" and "sequence". Default is `test_sequences.csv`.
  --mmseqs_results_file MMSEQS_RESULTS_FILE
                        MMseqs output with query,target,evalue,qstart,qend,tstart,tend,qaln,taln.
  --outfile OUTFILE     Name of the output CSV file. Default is `templates.csv`.
  --dataset_name, --name DATASET_NAME
                        full dataset_name, tag for csvs
  -o, --outdir OUTDIR   Where to save output CSVs (Default ./)
  --max_templates MAX_TEMPLATES
                        Maximum number of templates for target. Default is 5. Use 40 to prepare solution
  --cif_dir CIF_DIR     Directory holding cif.gz files, pdb_release_dates_NA.csv, and pdb_seqres_NA.fasta
  --skip_temporal_cutoff
                        Disable tests of temporal cutoff
  --start_idx START_IDX
                        Start index (1,2,...) of test_sequences to work on, for parallelization. Default: 0 (do all sequences).
  --end_idx END_IDX     End index (1,2,...) of test_sequences to work on, for parallelization. Default: 0 (do all sequences).
  --id_map ID_MAP       CSV file with fields `orig` and `new` for mapping original target IDs to new target IDs. Default is `` (no mapping).
```

## How to run MMseqs2

To get the mmseqs file, install [MMseqs2](https://github.com/soedinglab/MMseqs2), and within `example/`, use command lines like

```
mmseqs createdb ../PDB_RNA/pdb_seqres_NA.fasta pdb_seqres_NA  --dbtype 2
mmseqs easy-search validation_sequences.fasta pdb_seqres_NA  validation_Result.txt tmp --search-type 3 --format-output "query,target,evalue,qstart,qend,tstart,tend,qaln,taln"
```




