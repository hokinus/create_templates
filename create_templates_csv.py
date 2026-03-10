#!/usr/bin/env python3
from collections import defaultdict

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.sequence as seq
import numpy as np
import pandas as pd
import os
import gzip
import sys
import csv
import warnings
import argparse
from copy import deepcopy
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm

parser = argparse.ArgumentParser(description="Prepare templates.csv file similar to labels.csv but with MMseqs2-identified templates")
parser.add_argument('-s', '--sequences_file',
                    default="validation_sequences.csv",
                    help='CSV file with columns including "target_id" and "sequence". Default is `test_sequences.csv`.')
parser.add_argument('--mmseqs_results_file',
                    default= '',
                    help='MMseqs output with query,target,evalue,qstart,qend,tstart,tend,qaln,taln.')
parser.add_argument('--outfile',
                    default='',
                    help='Name of the output CSV file. Default is `templates.csv`.')
parser.add_argument('--dataset_name','--name',
                    default='',
                    help='full dataset_name, tag for csvs')
parser.add_argument(
    "-o", "--outdir", default=None, help="Where to save output CSVs (Default ./)"
)
parser.add_argument('--max_templates',
                    default=40, type=int,
                    help='Maximum number of templates for target. Default is 5. Use 40 to prepare solution')
parser.add_argument('--cif_dir',
                    default='',
                    help='Directory holding cif.gz files, pdb_release_dates_NA.csv, and pdb_seqres_NA.fasta')
parser.add_argument('--skip_temporal_cutoff',
                    action='store_true',
                    help='Disable tests of temporal cutoff')
parser.add_argument(
    "--temporal_cutoff",
    default="",
    help="Fixed temporal cutoff date (YYYY-MM-DD) applied to all targets, "
    "overriding per-target cutoffs from the sequences file.",
)
parser.add_argument('--start_idx',
                    default=0, type=int,
                    help='Start index (1,2,...) of test_sequences to work on, for parallelization. Default: 0 (do all sequences).' )
parser.add_argument('--end_idx',
                    default=0, type=int,
                    help='End index (1,2,...) of test_sequences to work on, for parallelization. Default: 0 (do all sequences).' )
parser.add_argument('--id_map',
                    default='',
                    help='CSV file with fields `orig` and `new` for mapping original target IDs to new target IDs. Default is `` (no mapping).')
parser.add_argument(
    "--num_workers",
    default=1,
    type=int,
    help="Number of parallel worker processes for processing targets. Default is 1,sequential run.",
)


def clean_res_name( res_name ):
    if res_name in ['A', 'C', 'G', 'U']:
        return res_name
    else: # can be modified residue with 3-letter name.
        return 'X'

#all atoms
ALL_ATOMS=["P","OP1","OP2","O5'","O3'","C1'","C2'","O2'","C3'","C4'","O4'","C5'","N1","C2","O2","N3","C4","N4","C5","C6","O4","N9","N7","C8","N6","N2","O6"]
C1PRIME_KEY="C1\'"

def extract_title_release_date(pdbx_file):
    block = pdbx_file.block
    
    possible_title_fields = [
        'struct',
        'entry',
        'struct_keywords'
    ]

    pdb_title = None
    for field in possible_title_fields:
        if field in block:
            if field == 'struct' and 'title' in block[field]:
                pdb_title = block[field]['title'].as_item()
                break
            elif field == 'entry' and 'title' in block[field]:
                pdb_title = block[field]['title'].as_item()
                break
            elif field == 'struct_keywords' and 'pdbx_keywords' in block[field]:
                pdb_title = block[field]['pdbx_keywords'].as_item()
                break

    possible_date_fields = [
        ('pdbx_database_status', 'initial_release_date'),
        ('pdbx_database_status', 'recvd_initial_deposition_date'),
        ('database_PDB_rev', 'date')
    ]

    release_date = None
    for category, field in possible_date_fields:
        if category in block and field in block[category]:
            date_data = block[category][field]
            if hasattr(date_data, 'as_array'):
                release_date = date_data.as_array()[0]  # Take the first date if it's an array
            else:
                release_date = date_data.as_item()
            break

    return pdb_title, release_date


def extract_rna_sequence(pdbx_file, chain_id):
    block = pdbx_file.block

    # Extract _pdbx_poly_seq_scheme information
    if "pdbx_poly_seq_scheme" in block:
        scheme = block["pdbx_poly_seq_scheme"]

        strand_id = (
            scheme["pdb_strand_id"].as_array() if "pdb_strand_id" in scheme else []
        )

        mask = strand_id == chain_id

        # note use of auth_seq_num instead of pdb_seq_num since that is what Biotite uses for res_id
        pdb_mon_id = scheme["pdb_mon_id"].as_array()[mask]

        pdb_chain_seq_nums = scheme["seq_id"].as_array()[mask]
        pdb_chain_ins_codes = scheme["pdb_ins_code"].as_array()[mask]
        print(pdb_chain_seq_nums.dtype)
        if not np.issubdtype(pdb_chain_seq_nums.dtype, np.number):
            try:
                num_mask = np.strings.isnumeric(pdb_chain_seq_nums)
                pdb_chain_seq_nums[~num_mask] = 0
                pdb_chain_seq_nums = pdb_chain_seq_nums.astype(int)
            except AttributeError as e:
                raise RuntimeError(
                    f"Error processing sequence numbers for chain {chain_id}: {e}"
                )
        # Sanitize residue numbers
        pdb_chain_seq_nums[len(pdb_mon_id) :] = 0

        # Sanitize insertion codes
        pdb_chain_ins_codes = np.strings.replace(pdb_chain_ins_codes, ".", "")
        pdb_chain_ins_codes[len(pdb_mon_id) :] = ""

        pdb_chain_sequence = np.vectorize(clean_res_name)(pdb_mon_id)

        return pdb_chain_sequence, pdb_chain_seq_nums, pdb_chain_ins_codes
    else:
        return "", [], []


def get_coord_labels(
    pdbx_file, chain_id, chain_sequence, chain_seq_nums, chain_ins_codes
):
    """
    Extract coordinates for an RNA chain based on a reference sequence alignment.

    This function uses Biotite to parse a CIF file, finds the specified chain,
    and extracts coordinates for RNA residues if there is indeed a C1' (nan's otherwise).

    Parameters:
    cif_path (str): Path to the CIF file.
    chain_id (str): Chain identifier in the CIF file.
    chain_sequence (str): sequence of target derived from polyx_ fields, used for output and as sanity check.
    chain_seq_nums  (list of strings): numbers of residues in PDB (auth_seq_num)
    chain_ins_codes (list of strings): ins_codes in PDB (needed to ensure unique lookup!)

    Returns:
    list of tuples: Each tuple contains (resname, resid, xyz, pdb_seq_num), where:
                    resname: Residue name (A, C, G, or U) from the reference sequence
                    resid: Residue ID (1, 2, 3, ...) based on position in reference sequence
                    xyz: dictionary of xyz coords for all 26 (heavy) atom names,
                        P,OP1,OP2,O5',O3',C1',C2',C3',C4',O4',C5',N1,C2,O2,N3,C4,N4,C5,C6,O4,N9,N7,C8,N6,N2,O6
                    pdb_info: (author seq num, ins code, resname)

    The length of the returned list is equal to the length of input chain_sequence.
    """

    # Get structure using biotite, use label_ ids
    structure = pdbx.get_structure(pdbx_file, use_author_fields=False, model=1)

    # Filter for the specified chain
    chain_filter = (
        structure.chain_id == chain_id
    )  # TODO: make sure chain ids match label_asym_id
    chain_atoms = structure[chain_filter]

    # We need to match observed residues to requested sequence
    # We can match residues directly based on ids, both use label_ ids
    full_sequence = pd.DataFrame(
        {
            "resnum": chain_seq_nums,
            "ins_code": chain_ins_codes,
            "resname": chain_sequence,
        }
    )
    # Expand full sequence to include all_atoms, use long form (each atom is a row)
    all_atoms_df = pd.DataFrame({"atom_name": ALL_ATOMS, "temp_key": 1})
    full_sequence["temp_key"] = 1

    full_sequence_all_atom = full_sequence.merge(
        all_atoms_df, on="temp_key", how="outer"
    )

    # Create a DataFrame from chain_atoms to merge with
    chain_atoms_df = pd.DataFrame(
        {
            "resname": chain_atoms.res_name,
            "resnum": chain_atoms.res_id,
            "ins_code": chain_atoms.ins_code,
            "atom_name": chain_atoms.atom_name,
            "x": chain_atoms.coord[:, 0],
            "y": chain_atoms.coord[:, 1],
            "z": chain_atoms.coord[:, 2],
        }
    )

    # Deduplicate to handle alternate conformations (keep first occurrence per resnum+atom_name)
    chain_atoms_df = chain_atoms_df.drop_duplicates(subset=["resnum", "atom_name"], keep="first")

    full_sequence_all_atom = full_sequence_all_atom.merge(
        chain_atoms_df, on=["resnum", "atom_name"], how="left", suffixes=("", "_obs")
    )

    # Build result: list of (resname, resid, xyz, pdb_info) tuples — one per residue
    result = []
    n_atoms = len(ALL_ATOMS)
    for i in range(len(chain_sequence)):
        chunk = full_sequence_all_atom.iloc[i * n_atoms : (i + 1) * n_atoms]
        xyz = {atom: (np.nan, np.nan, np.nan) for atom in ALL_ATOMS}
        obs_resname = ""
        for _, row in chunk.iterrows():
            if not pd.isna(row["x"]):
                xyz[row["atom_name"]] = (float(row["x"]), float(row["y"]), float(row["z"]))
                if not obs_resname and "resname_obs" in chunk.columns and not pd.isna(row["resname_obs"]):
                    obs_resname = str(row["resname_obs"])
        ins_c = chain_ins_codes[i] if chain_ins_codes[i] else " "
        result.append((chain_sequence[i], i + 1, xyz, (chain_seq_nums[i], ins_c, obs_resname)))

    return result

def get_target_coord_data( chain_coord_data, alignment ):
    '''
    Inputs
      chain_coord_data = coordinates and other info, read out from PDB file for the chain
      alignment = two strings that map chain to target with gaps as '-'

    Output
      target_coord_data = coordinates for target sequence, with gaps filled with nan.
    '''
    target_coord_data = []
    chain_pos = -1
    target_pos = -1
    xyz_blank = { atom:(np.nan,np.nan,np.nan) for atom in ALL_ATOMS}
    for (chain_res,target_res) in zip(alignment[0],alignment[1]):
        if chain_res  != '-':
            chain_pos += 1
        if target_res != '-':
            target_pos += 1
            if chain_res != '-':
                coord_data=chain_coord_data[ chain_pos ]
                target_coord_data.append( (target_res, target_pos+1,coord_data[2],coord_data[3]) )
            else:
                target_coord_data.append( (target_res, target_pos+1,xyz_blank,(-1e18,' ','') ) )

    return target_coord_data

def is_before_or_on(d1, d2):
    date1 = pd.to_datetime(d1)
    date2 = pd.to_datetime(d2)
    return date1 <= date2

def read_id_map(id_map_file):
    if len(id_map_file)==0: return None
    id_map = {}
    try:
        with open(id_map_file, newline='') as f:
            reader = csv.DictReader(f)
            if 'orig' not in reader.fieldnames or 'new' not in reader.fieldnames:
                print("Warning: ID map file does not contain the fields 'orig' and 'new'. Using original IDs instead.")
                return id_map
            for row in reader:
                id_map[row['orig']] = row['new']
    except FileNotFoundError:
        print(f"Warning: ID map file {id_map_file} not found. Using original IDs instead.", file=sys.stderr)
    except Exception as exc:
        print(f"Error reading {id_map_file}: {exc}", file=sys.stderr)
    return id_map

def read_release_dates( release_data_file ):
    release_dates = {}
    # must have format Entry ID, Release Date
    with open(release_data_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            release_dates[row['Entry ID']] = row['Release Date']

    return release_dates

def read_cif_file(cif_path):
    if cif_path.endswith(".gz"):
        with gzip.open(cif_path, "rt") as cif_file:
            pdbx_file = pdbx.CIFFile.read(cif_file)
    else:
        pdbx_file = pdbx.CIFFile.read(cif_path)
    return pdbx_file

def process_target(target, sequence, temporal_cutoff, aln_lines, skip_temporal_cutoff,
                   MAX_TEMPLATES, cif_dir, release_dates, id_map):
    """Process a single target and return (output_labels, output_allatom_labels)."""
    output_labels = []
    output_allatom_labels = []

    templates = []
    template_coord_data = []
    for aln_line in aln_lines:
        if len(aln_line)!=9: continue # some kind of overflow in some alignments?

        query,template,eval,qstart,qend,tstart,tend,qaln,taln = aln_line

        if int(qend) < int(qstart):
            continue  # aligned to reverse complement!

        pdb_id, chain_id = template.split("_")

        cif_path = os.path.join(cif_dir, f"{pdb_id.upper()}.cif.gz")
        if not os.path.isfile(cif_path):
            cif_path = os.path.join(cif_dir, f"{pdb_id.lower()}.cif")  # kaggle style
            if not os.path.isfile(cif_path):
                continue  # occasional alignment to DNA, ignore!

        release_date = release_dates[pdb_id.upper()]  # pulled from PDB server

        if not skip_temporal_cutoff and is_before_or_on(temporal_cutoff, release_date):
            continue

        cif_file = read_cif_file(cif_path)
        # these release dates in the CIF files can be buggy!
        title, release_date_unreliable = extract_title_release_date(cif_file)

        print("\n", target, temporal_cutoff, "   ", template)
        if title:
            print(f"PDB Title: {title}")
        if release_date:
            print(f"PDB Release Date: {release_date}")

        chain_sequence, chain_seq_nums, chain_ins_codes = extract_rna_sequence(
            cif_file, chain_id
        )

        alignment = []
        qstart = int(qstart)
        qend = int(qend)
        tstart = int(tstart)
        tend = int(tend)
        alignment.append(
            sequence[: (qstart - 1)] + "-" * (tstart - 1) + qaln + sequence[qend:]
        )
        alignment.append(
            "-" * (qstart - 1)
            + "X" * (tstart - 1)
            + taln
            + "-" * (len(sequence) - qend)
        )
        print(alignment[0], "query")
        print(alignment[1], "template")
        chain_coord_data = get_coord_labels(
            cif_file, chain_id, chain_sequence, chain_seq_nums, chain_ins_codes
        )

        coord_data = get_target_coord_data(
            chain_coord_data, (alignment[1], alignment[0])
        )

        if len(coord_data) != len(sequence):
            print(
                "WARNING! len(coord_data) != len(sequence)",
                "len coord_data",
                len(coord_data),
                "len sequence",
                len(sequence),
                "qstart",
                qstart,
                "len qaln",
                len(qaln),
                "qend",
                qend,
            )
            continue

        templates.append(template)
        template_coord_data.append(coord_data)

        if len(templates) >= MAX_TEMPLATES:
            break

    print("Found", len(templates), "templates for", target, "\n")

    mapped_target = target
    if id_map is not None:
        mapped_target = id_map[target]

    for i in range(len(sequence)):
        output_label = {
            "ID": f"{mapped_target}_{i + 1}",
            "resname": sequence[i],
            "resid": i + 1,
        }
        output_allatom_label = deepcopy(output_label)

        for n in range(len(templates)):
            template = templates[n]
            res, resid, xyz, pdb_info = template_coord_data[n][i]
            assert resid == i + 1
            output_label[f"x_{n + 1}"] = xyz[C1PRIME_KEY][0]
            output_label[f"y_{n + 1}"] = xyz[C1PRIME_KEY][1]
            output_label[f"z_{n + 1}"] = xyz[C1PRIME_KEY][2]

            for atom in ALL_ATOMS:
                output_allatom_label.update(
                    {
                        f"{atom}_x_{n + 1}": xyz[atom][0],
                        f"{atom}_y_{n + 1}": xyz[atom][1],
                        f"{atom}_z_{n + 1}": xyz[atom][2],
                    }
                )
            output_allatom_label.update(
                {
                    f"pdb_id_{n + 1}": template,
                    f"pdb_seq_num_{n + 1}": int(pdb_info[0]),
                    f"pdb_ins_code_{n + 1}": pdb_info[1],
                    f"pdb_resname_{n + 1}": pdb_info[2],
                }
            )

        for n in range(len(templates), MAX_TEMPLATES):
            output_label[f"x_{n + 1}"] = np.nan
            output_label[f"y_{n + 1}"] = np.nan
            output_label[f"z_{n + 1}"] = np.nan

            for atom in ALL_ATOMS:
                output_allatom_label.update(
                    {
                        f"{atom}_x_{n + 1}": np.nan,
                        f"{atom}_y_{n + 1}": np.nan,
                        f"{atom}_z_{n + 1}": np.nan,
                    }
                )
            output_allatom_label.update(
                {
                    f"pdb_id_{n + 1}": "",
                    f"pdb_seq_num_{n + 1}": np.nan,
                    f"pdb_ins_code_{n + 1}": "",
                    f"pdb_resname_{n + 1}": "",
                }
            )

        output_labels.append(output_label)
        output_allatom_labels.append(output_allatom_label)

    return output_labels, output_allatom_labels

def process_and_save_target(
    target,
    sequence,
    temporal_cutoff,
    aln_lines,
    skip_temporal_cutoff,
    MAX_TEMPLATES,
    cif_dir,
    release_dates,
    id_map,
    outdir,
):
    output_labels, output_allatom_labels = process_target(
        target,
        sequence,
        temporal_cutoff,
        aln_lines,
        skip_temporal_cutoff,
        MAX_TEMPLATES,
        cif_dir,
        release_dates,
        id_map,
    )
    output_template_labels_to_csv(
        output_labels,
        output_allatom_labels,
        [target],
        outdir=outdir,
        dataset_name=target,
        start_idx=0,
        end_idx=0,
    )

def preprocess_inputs(
    sequences_file,
    mmseqs_results_file,
    cif_dir,
    id_map_file="",
    start_idx=0,
    end_idx=0,
    fixed_temporal_cutoff="",
):
    """Read sequences, MMseqs alignments, release dates, and build task arguments.

    Args:
        fixed_temporal_cutoff: If non-empty, overrides the per-target temporal cutoff
            from the sequences file with this single date (YYYY-MM-DD) for all targets.

    Returns:
        task_args: list of (target, sequence, cutoff, aln_lines) tuples
        targets: full list of target IDs from the sequences file
        cif_dir: resolved CIF directory path
        release_dates: dict mapping PDB ID -> release date
        id_map: dict mapping original -> new target IDs, or None
        start_idx, end_idx: resolved 1-based slice indices
    """
    if len(cif_dir) == 0:
        dir_name = os.path.dirname(os.path.abspath(sys.argv[0]))
        cif_dir = dir_name + "/PDB_RNA"

    df = pd.read_csv(sequences_file)
    targets = df["target_id"].to_list()
    sequences = df["sequence"].to_list()
    temporal_cutoffs = df["temporal_cutoff"].to_list()

    if start_idx == 0 and end_idx == 0:  # do all targets by default
        start_idx = 1
        end_idx = len(targets)

    selected_targets = targets[start_idx - 1 : end_idx]
    selected_sequences = sequences[start_idx - 1 : end_idx]
    selected_cutoffs = temporal_cutoffs[start_idx - 1 : end_idx]

    if fixed_temporal_cutoff:
        print(
            f"Overriding per-target temporal cutoffs with fixed cutoff: {fixed_temporal_cutoff}"
        )
        selected_cutoffs = [fixed_temporal_cutoff] * len(selected_targets)

    selected_targets_set = set(selected_targets)
    aln_lines = defaultdict(list)
    for line in tqdm(
        open(mmseqs_results_file), desc="Reading MMseqs results", file=sys.stderr
    ):
        # query,template,eval,qstart,qend,tstart,tend,qaln,taln
        parts = line.strip().split()
        if parts:
            query = parts[0]
            if query not in selected_targets_set:
                continue
            aln_lines[query].append(parts)

    id_map = read_id_map(id_map_file)
    release_dates = read_release_dates(cif_dir + "/pdb_release_dates_NA.csv")

    task_args = [
        (target, sequence, cutoff, aln_lines[target])
        for target, sequence, cutoff in zip(
            selected_targets, selected_sequences, selected_cutoffs
        )
    ]

    return task_args, targets, cif_dir, release_dates, id_map, start_idx, end_idx


def get_template_labels_serial(
    sequences_file,
    mmseqs_results_file,
    skip_temporal_cutoff,
    MAX_TEMPLATES,
    cif_dir,
    id_map_file="",
    start_idx=0,
    end_idx=0,
    fixed_temporal_cutoff="",
):
    """Preprocess inputs, process all targets sequentially, and return results.

    Args:
        fixed_temporal_cutoff: If non-empty, overrides per-target cutoffs for all targets.

    Returns:
        output_labels: list of per-residue C1'-only label dicts
        output_allatom_labels: list of per-residue all-atom label dicts
        targets: full list of target IDs from the sequences file
        start_idx, end_idx: resolved 1-based slice indices
    """
    task_args, targets, cif_dir, release_dates, id_map, start_idx, end_idx = (
        preprocess_inputs(
            sequences_file,
            mmseqs_results_file,
            cif_dir,
            id_map_file,
            start_idx,
            end_idx,
            fixed_temporal_cutoff,
        )
    )
    print("Running in sequential mode")
    worker_fn = partial(
        process_target,
        skip_temporal_cutoff=skip_temporal_cutoff,
        MAX_TEMPLATES=MAX_TEMPLATES,
        cif_dir=cif_dir,
        release_dates=release_dates,
        id_map=id_map,
    )
    results = [
        worker_fn(*args)
        for args in tqdm(task_args, desc="Processing targets", file=sys.stderr)
    ]

    output_labels = []
    output_allatom_labels = []
    for labels, allatom_labels in results:
        output_labels.extend(labels)
        output_allatom_labels.extend(allatom_labels)

    print(f'Completed {len(task_args)} targets\n')
    return output_labels, output_allatom_labels, targets, start_idx, end_idx


def get_template_labels_parallel(
    sequences_file,
    mmseqs_results_file,
    skip_temporal_cutoff,
    MAX_TEMPLATES,
    cif_dir,
    id_map_file="",
    start_idx=0,
    end_idx=0,
    num_workers=1,
    outdir="output",
    fixed_temporal_cutoff="",
):
    """Preprocess inputs and process targets in parallel, writing one output file per target."""
    task_args, targets, cif_dir, release_dates, id_map, start_idx, end_idx = (
        preprocess_inputs(
            sequences_file,
            mmseqs_results_file,
            cif_dir,
            id_map_file,
            start_idx,
            end_idx,
            fixed_temporal_cutoff,
        )
    )
    print("Running in parallel mode, output one file per target")
    worker_fn = partial(
        process_and_save_target,
        skip_temporal_cutoff=skip_temporal_cutoff,
        MAX_TEMPLATES=MAX_TEMPLATES,
        cif_dir=cif_dir,
        release_dates=release_dates,
        id_map=id_map,
        outdir=outdir,
    )
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(worker_fn, *args): i for i, args in enumerate(task_args)
        }
        with tqdm(
            total=len(task_args), desc="Processing targets", file=sys.stderr
        ) as pbar:
            for future in as_completed(futures):
                pbar.update(1)
    print(f"Completed {len(task_args)} targets\n")


# Create a DataFrame and write to CSV
def output_csv(output_data, outfile):
    df = pd.DataFrame(output_data)
    df.to_csv(outfile, index=False)
    print(f"Output written to {outfile}")


def output_template_labels_to_csv(
    output_labels,
    output_allatom_labels,
    targets,
    outdir="",
    outfile=None,
    outfile_allatom=None,
    dataset_name="",
    start_idx=0,
    end_idx=0,
):
    assert not (outfile is not None and len(dataset_name) > 0)

    os.makedirs(outdir, exist_ok=True)
    if outdir[-1] != "/":
        outdir += "/"
    split_tag = ""
    if start_idx > 0:
        num_digits = len(str(len(targets)))
        split_tag = f".{start_idx:0{num_digits}d}_{end_idx:0{num_digits}d}"

    if outfile is None:
        if len(dataset_name) == 0:
            dataset_name = "test"
        outfile = f"{outdir}{dataset_name}.templates{split_tag}.csv"
        outfile_allatom = f"{outdir}{dataset_name}.allatom_templates{split_tag}.csv"
    else:
        if outfile_allatom is None:
            if outfile.count("labels.csv") > 1:
                outfile_allatom = outfile.replace("labels.csv", "allatom.csv")
            elif outfile.endswith(".csv"):
                outfile_allatom = outfile.replace(".csv", ".allatom.csv")
            else:
                outfile_allatom = outfile + ".allatom.csv"

    output_csv(output_labels, outfile)
    output_csv(output_allatom_labels, outfile_allatom)


if __name__ == "__main__":
    args = parser.parse_args()

    sequences_file = args.sequences_file
    mmseqs_results_file = args.mmseqs_results_file
    MAX_TEMPLATES = args.max_templates
    cif_dir = args.cif_dir
    id_map_file = args.id_map
    start_idx = args.start_idx
    end_idx = args.end_idx
    skip_temporal_cutoff = args.skip_temporal_cutoff
    num_workers = args.num_workers
    parallel = num_workers > 1

    if parallel:
        if args.outdir is None:
            print(
                "Error: --outdir must be specified when using parallel workers",
                file=sys.stderr,
            )
            sys.exit(1)
        get_template_labels_parallel(
            sequences_file,
            mmseqs_results_file,
            skip_temporal_cutoff,
            MAX_TEMPLATES,
            cif_dir,
            id_map_file,
            start_idx,
            end_idx,
            num_workers,
            outdir=args.outdir,
            fixed_temporal_cutoff=args.temporal_cutoff,
        )
    else:
        output_labels, output_allatom_labels, targets, start_idx, end_idx = (
            get_template_labels_serial(
                sequences_file,
                mmseqs_results_file,
                skip_temporal_cutoff,
                MAX_TEMPLATES,
                cif_dir,
                id_map_file,
                start_idx,
                end_idx,
                fixed_temporal_cutoff=args.temporal_cutoff,
            )
        )
        output_template_labels_to_csv(
            output_labels,
            output_allatom_labels,
            targets,
            args.outdir,
            f"{args.outdir}/{args.outfile}",
            args.dataset_name,
            start_idx,
            end_idx,
        )
