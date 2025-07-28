#!/usr/bin/env python3
from Bio import SeqIO,PDB,BiopythonWarning
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.PDB import MMCIFParser
import numpy as np
import pandas as pd
import os
import gzip
import sys
import csv
import warnings
import argparse
from datetime import datetime

# Suppress warnings
warnings.simplefilter('ignore', BiopythonWarning)

parser = argparse.ArgumentParser(description="Prepare templates.csv file similar to labels.csv but with MMseqs2-identified templates")
parser.add_argument('-s', '--sequences_file',
                    default="validation_sequences.csv",
                    help='CSV file with columns including "target_id" and "sequence". Default is `test_sequences.csv`.')
parser.add_argument('--mmseqs_results_file',
                    default= '',
                    help='MMseqs output with query,target,evalue,qstart,qend,tstart,tend,qaln,taln.')
parser.add_argument('-o', '--outfile',
                    default='templates.csv',
                    help='Name of the output CSV file. Default is `templates.csv`.')
parser.add_argument('--max_templates',
                    default=40, type=int,
                    help='Maximum number of templates for target. Default is 5. Use 40 to prepare solution')
parser.add_argument('--cif_dir',
                    default='',
                    help='Directory holding cif.gz files, pdb_release_dates_NA.csv, and pdb_seqres_NA.fasta')
parser.add_argument('--skip_temporal_cutoff',
                    action='store_true',
                    help='Disable tests of temporal cutoff')
parser.add_argument('--start_idx',
                    default=0, type=int,
                    help='Start index (1,2,...) of test_sequences to work on, for parallelization. Default: 0 (do all sequences).' )
parser.add_argument('--end_idx',
                    default=0, type=int,
                    help='End index (1,2,...) of test_sequences to work on, for parallelization. Default: 0 (do all sequences).' )
parser.add_argument('--id_map',
                    default='',
                    help='CSV file with fields `orig` and `new` for mapping original target IDs to new target IDs. Default is `` (no mapping).')
args = parser.parse_args()

sequences_file = args.sequences_file
mmseqs_results_file = args.mmseqs_results_file
outfile = args.outfile
MAX_TEMPLATES = args.max_templates
cif_dir = args.cif_dir
id_map_file = args.id_map
start_idx = args.start_idx
end_idx = args.end_idx

if len(cif_dir) == 0:
    dir_name = os.path.dirname( os.path.abspath( sys.argv[0] ) )
    cif_dir = dir_name+'/PDB_RNA'

def clean_res_name( res_name ):
    if res_name in ['A', 'C', 'G', 'U']:
        return res_name
    else: # can be modified residue with 3-letter name.
        return 'X'

def extract_title_release_date( cif_path ):

    if cif_path.endswith('.gz'):
        with gzip.open(cif_path, 'rt') as cif_file:
            mmcif_dict = MMCIF2Dict(cif_file)
    else:
        mmcif_dict = MMCIF2Dict(cif_path)

    possible_title_fields = [
        '_struct.title',
        '_entry.title',
        '_struct_keywords.pdbx_keywords'
    ]

    pdb_title = None
    for field in possible_title_fields:
        if field in mmcif_dict:
            pdb_title = mmcif_dict[field]
            if isinstance(pdb_title, list):
                pdb_title = ' '.join(pdb_title)
            break

    possible_date_fields = [
        '_pdbx_database_status.initial_release_date',
        '_pdbx_database_status.recvd_initial_deposition_date',
        '_database_PDB_rev.date'
    ]

    release_date = None
    for field in possible_date_fields:
        if field in mmcif_dict:
            release_date = mmcif_dict[field]
            if isinstance(release_date, list):
                release_date = release_date[0]  # Take the first date if it's a list
            break

    return pdb_title, release_date


def extract_rna_sequence(cif_path,chain_id):

    if cif_path.endswith('.gz'):
        with gzip.open(cif_path, 'rt') as cif_file:
            mmcif_dict = MMCIF2Dict(cif_file)
    else:
        mmcif_dict = MMCIF2Dict(cif_path)

    pdb_sequence = None
    pdb_chain_id = None
    chain_seq_nums = None

    # Extract _pdbx_poly_seq_scheme information
    strand_id  = mmcif_dict.get('_pdbx_poly_seq_scheme.pdb_strand_id',[])
    mon_id     = mmcif_dict.get('_pdbx_poly_seq_scheme.mon_id',[])
    pdb_mon_id = mmcif_dict.get('_pdbx_poly_seq_scheme.pdb_mon_id',[])
    pdb_seq_num = mmcif_dict.get('_pdbx_poly_seq_scheme.pdb_seq_num',[])
    chain_ids = list(set(strand_id))
    seq_chains = []

    full_sequence = ''
    pdb_chain_sequence = ''
    pdb_chain_seq_nums = []
    for (strand,mon,pdb_mon,pdb_num) in zip(strand_id,mon_id,pdb_mon_id,pdb_seq_num):
        if strand==chain_id:
            full_sequence += clean_res_name( mon )
            pdb_chain_sequence += clean_res_name( pdb_mon )
            pdb_chain_seq_nums.append( pdb_num)

    #print(full_sequence)
    #print(pdb_chain_sequence)
    #print(pdb_chain_seq_nums)

    return full_sequence,pdb_chain_sequence,pdb_chain_seq_nums

def get_c1prime_labels(cif_path, chain_id, alignment, chain_seq_nums):
    """
    Extract C1' coordinates for an RNA chain based on a reference sequence alignment.

    This function uses Biopython to parse a CIF file, finds the specified chain,
    and extracts C1' coordinates for RNA residues. It aligns these coordinates
    with a reference sequence, handling gaps and missing residues.

    Parameters:
    cif_path (str): Path to the CIF file.
    chain_id (str): Chain identifier in the CIF file.
    alignment (list): A list containing two elements:
                      alignment[0]: List of residues for the reference sequence (A,C,G,U,-)
                      alignment[1]: List of residues for the chain sequence (A,C,G,U,X,-)
    chain_seq_nums (list): numbers of residues in PDB

    Returns:
    list of tuples: Each tuple contains (resname, resid, x, y, z), where:
                    resname: Residue name (A, C, G, or U) from the reference sequence
                    resid: Residue ID (1, 2, 3, ...) based on position in reference sequence
                    x, y, z: C1' coordinates (nan for missing residues/atoms)

    The length of the returned list is equal to the number of non-gap residues
    in the reference sequence.
    """
    # Parse the CIF file
    parser = MMCIFParser()
    if cif_path.endswith('.gz'):
        with gzip.open(cif_path, 'rt') as gz_file:
            structure = parser.get_structure('RNA', gz_file )
    else:
        structure = parser.get_structure('RNA', cif_path)

    # Get the specified chain
    chain = structure[0][chain_id]

    # getting residues out of chain is complex -- easier to get a list ahead of time.
    residues = {}
    for residue in chain: residues[ residue.id[1] ] = residue

    chain_seq = ''.join( [clean_res_name(residue.get_resname()) for residue in chain ] )
    #print(chain_seq)

    # Initialize the result list
    result = []

    # Counter for residue ID in reference sequence
    ref_resid = 0
    chain_idx = 0
    for ref_res, chain_res in zip(alignment[0], alignment[1]):
        if chain_res != '-': chain_idx += 1
        if ref_res != '-':
            ref_resid += 1
            if chain_res == '-': # or chain_res == 'X':
                # Missing residue in chain or unknown residue
                result.append((ref_res, ref_resid, np.nan, np.nan, np.nan, -1e18))
            else:
                # Find the corresponding residue in the chain
                try:
                    #chain_seq_num = int(chain_seq_nums[ref_resid-1])
                    chain_seq_num = int(chain_seq_nums[chain_idx-1])
                    residue = residues[chain_seq_num]
                    c1_prime = residue['C1\'']
                    coords = c1_prime.coord
                    if residue.get_resname() != chain_res:
                        print( 'WARNING!',ref_resid,chain_idx,chain_seq_num,residue.get_resname(),chain_res)
                    result.append((ref_res, ref_resid, coords[0], coords[1], coords[2], residue.id[1]))
                except KeyError:
                    # C1' atom not found
                    result.append((ref_res, ref_resid, np.nan, np.nan, np.nan, -1e18))
                except Exception as e:
                    # Any other error (e.g., residue not found)
                    result.append((ref_res, ref_resid, np.nan, np.nan, np.nan, -1e18))

    return result

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


# Prepare to collect output data
output_labels = []

# Read the FASTA file
df = pd.read_csv( sequences_file )
targets = df['target_id'].to_list()
sequences = df['sequence'].to_list()
temporal_cutoffs = df['temporal_cutoff'].to_list()

aln_lines = []
for line in open( mmseqs_results_file ).readlines():
    # query,template,eval,qstart,qend,tstart,tend,qaln,taln
    aln_lines.append( line.strip().split() )

id_map = read_id_map( id_map_file )

release_dates = read_release_dates( cif_dir + '/pdb_release_dates_NA.csv' )

if start_idx == 0 and end_idx == 0: # do all targets by default
    start_idx = 1
    end_idx = len(targets)

num_targets = 0
count = 0
for target,sequence,temporal_cutoff in zip(targets,sequences,temporal_cutoffs):
    count += 1
    if (count < start_idx) or (count > end_idx): continue


    # look for alignments and fill out C1' templates
    templates = []
    for aln_line in aln_lines:
        if len(aln_line)!=9: continue # some kind of overflow in some alignments?

        query,template,eval,qstart,qend,tstart,tend,qaln,taln = aln_line

        if query != target: continue

        if int(qend)<int(qstart): continue # aligned to reverse complement!

        pdb_id,chain_id = template.split('_')

        release_date = release_dates[pdb_id.upper()] # pulled from PDB server

        if not args.skip_temporal_cutoff and is_before_or_on(temporal_cutoff,release_date): continue

        # need to do alignment
        cif_path = os.path.join(cif_dir, f'{pdb_id.upper()}.cif.gz')
        if not os.path.isfile( cif_path ):
            cif_path = os.path.join(cif_dir, f'{pdb_id.lower()}.cif') # kaggle style
            if not os.path.isfile( cif_path ): continue # occasional alignment to DNA, ignore!


        # these release dates in the CIF files can be buggy!
        title,release_date_unreliable = extract_title_release_date( cif_path )

        print('\n',target,temporal_cutoff,"   ",template)
        if title: print(f"PDB Title: {title}")
        if release_date: print(f"PDB Release Date: {release_date}")

        # sometimes there is a mismatch between PDB's fasta files and what's actually stored in coordinates,
        # so best to get the actual residue numbers for the chain
        chain_full_sequence,chain_sequence,chain_seq_nums = extract_rna_sequence(cif_path,chain_id)

        # get 3d data
        alignment = []
        qstart=int(qstart)
        qend=int(qend)
        tstart=int(tstart)
        tend=int(tend)
        alignment.append( sequence[:(qstart-1)] + '-'*(tstart-1) + qaln + sequence[qend:]  )
        alignment.append( '-'*(qstart-1)        + 'X'*(tstart-1) + taln + '-'*(len(sequence)-qend) )
        print( alignment[0] )
        print( alignment[1] )
        c1prime_data = get_c1prime_labels( cif_path, chain_id, alignment, chain_seq_nums )

        # mismatch in FASTA sequence and the polyx info in the CIF file
        if len(c1prime_data) != len(sequence):
            print( 'WARNING! len(c1prime_data) != len(sequence)', 'len c1prime_data', len(c1prime_data), 'len sequence', len(sequence), 'qstart',qstart,'len qaln',len(qaln),'qend',qend)
            continue

        templates.append( c1prime_data )

        if len(templates) >= MAX_TEMPLATES: break

    print( "Found", len(templates), "templates for", target,'\n' )

    mapped_target = target
    if not id_map is None: mapped_target = id_map[target]

    for i in range(len(sequence)):
        output_label = {
            "ID": f'{mapped_target}_{i+1}',
            "resname": sequence[i],
            "resid": i+1,
        }

        # output templates
        for n in range(len(templates)):
            res,resid,x,y,z,pdb_seqnum = templates[n][i]
            assert( resid == i+1 )
            output_label[ f"x_{n+1}" ] = x
            output_label[ f"y_{n+1}" ] = y
            output_label[ f"z_{n+1}" ] = z

        # pad with blank models
        for n in range(len(templates),MAX_TEMPLATES):
            output_label[ f"x_{n+1}" ] = np.nan
            output_label[ f"y_{n+1}" ] = np.nan
            output_label[ f"z_{n+1}" ] = np.nan
        output_labels.append( output_label )

    num_targets += 1
    # if num_targets > 1: break # for debug!


print(f'Completed {num_targets} targets\n')

# Create a DataFrame and write to CSV
def output_csv( output_data, outfile ):
    df = pd.DataFrame(output_data)
    df.to_csv(outfile, index=False)
    print(f"Output written to {outfile}")

output_csv( output_labels, outfile )
