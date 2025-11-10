#!/usr/bin/env python3
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
parser.add_argument('-o','--outdir',
                    default='./',
                    help='Where to save output CSVs (Default ./)' )
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


def clean_res_name( res_name ):
    if res_name in ['A', 'C', 'G', 'U']:
        return res_name
    else: # can be modified residue with 3-letter name.
        return 'X'

#all atoms
ALL_ATOMS=["P","OP1","OP2","O5'","O3'","C1'","C2'","O2'","C3'","C4'","O4'","C5'","N1","C2","O2","N3","C4","N4","C5","C6","O4","N9","N7","C8","N6","N2","O6"]
C1PRIME_KEY="C1\'"

def extract_title_release_date( cif_path ):
    
    if cif_path.endswith('.gz'):
        with gzip.open(cif_path, 'rt') as cif_file:
            pdbx_file = pdbx.CIFFile.read(cif_file)
    else:
        pdbx_file = pdbx.CIFFile.read(cif_path)

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


def extract_rna_sequence(cif_path, chain_id):
    
    if cif_path.endswith('.gz'):
        with gzip.open(cif_path, 'rt') as cif_file:
            pdbx_file = pdbx.CIFFile.read(cif_file)
    else:
        pdbx_file = pdbx.CIFFile.read(cif_path)

    block = pdbx_file.block
    
    pdb_sequence = None
    pdb_chain_id = None
    chain_seq_nums = None

    # Extract _pdbx_poly_seq_scheme information
    if 'pdbx_poly_seq_scheme' in block:
        scheme = block['pdbx_poly_seq_scheme']
        
        strand_id = scheme['pdb_strand_id'].as_array() if 'pdb_strand_id' in scheme else []
        mon_id = scheme['mon_id'].as_array() if 'mon_id' in scheme else []
        pdb_mon_id = scheme['pdb_mon_id'].as_array() if 'pdb_mon_id' in scheme else []
        pdb_seq_num = scheme['pdb_seq_num'].as_array() if 'pdb_seq_num' in scheme else []
        auth_seq_num = scheme['auth_seq_num'].as_array() if 'auth_seq_num' in scheme else []
        pdb_ins_code = scheme['pdb_ins_code'].as_array() if 'pdb_ins_code' in scheme else []
        
        chain_ids = list(set(strand_id))
        seq_chains = []

        full_sequence = ''
        pdb_chain_sequence = ''
        pdb_chain_seq_nums = []
        pdb_chain_ins_codes = []

        for (strand,mon,pdb_mon,pdb_num,auth_num,ins_code) in zip(strand_id,mon_id,pdb_mon_id,pdb_seq_num,auth_seq_num,pdb_ins_code):
            if strand==chain_id:
                full_sequence += clean_res_name( mon )
                pdb_chain_sequence += clean_res_name( pdb_mon )
                # note use of auth_seq_num instead of pdb_seq_num since that is what Biotite uses for res_id
                pdb_chain_seq_nums.append( auth_num )
                pdb_chain_ins_codes.append( ins_code )

    return full_sequence,pdb_chain_sequence,pdb_chain_seq_nums,pdb_chain_ins_codes

def get_coord_labels(cif_path, chain_id, chain_sequence, chain_seq_nums, chain_ins_codes):
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
    # Parse the CIF file
    if cif_path.endswith('.gz'):
        with gzip.open(cif_path, 'rt') as cif_file:
            pdbx_file = pdbx.CIFFile.read(cif_file)
    else:
        pdbx_file = pdbx.CIFFile.read(cif_path)
    
    # Get structure using biotite
    structure = pdbx.get_structure(pdbx_file, model=1)
    
    # Filter for the specified chain
    chain_filter = structure.chain_id == chain_id
    chain_atoms = structure[chain_filter]
    
    # Group atoms by residue
    residue_dict = {}
    for i, atom in enumerate(chain_atoms):
        res_key = (int(chain_atoms.res_id[i]), str(chain_atoms.ins_code[i]))
        if res_key not in residue_dict:
            residue_dict[res_key] = {}
        residue_dict[res_key][str(chain_atoms.atom_name[i])] = {
            'coord': chain_atoms.coord[i],
            'res_name': chain_atoms.res_name[i]
        }

    # Initialize the result list
    result = []

    assert( len( chain_sequence) == len( chain_seq_nums ) )
    for i, chain_res in enumerate(chain_sequence):
        chain_resid = i+1
        chain_seq_num  = int(chain_seq_nums[i]) if (chain_seq_nums[i].isdigit() and i < len(chain_seq_nums)) else 0
        chain_ins_code = chain_ins_codes[i].replace('.','')
        res_key = (chain_seq_num, chain_ins_code)

        xyz = { atom:(np.nan,np.nan,np.nan) for atom in ALL_ATOMS}
        res_info = (chain_res, chain_resid, xyz, (-1e18,' ','') ) # blank
        if res_key in residue_dict:
            residue_atoms = residue_dict[res_key]
            if "C1'" in residue_atoms:
                resname = residue_atoms["C1'"]['res_name']
                if chain_res != clean_res_name( resname ):
                    print( f'Warning! mismatch residue at {chain_resid}: target {chain_res} pdb {resname} chain_seq_num {chain_seq_num}' )

                for atom in ALL_ATOMS:
                    if atom in residue_atoms:
                        xyz[atom] = tuple(residue_atoms[atom]['coord'])

                res_info = (chain_res, chain_resid, xyz, (res_key[0], res_key[1], resname) )
        result.append(res_info)
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


def get_template_labels( sequences_file, mmseqs_results_file, skip_temporal_cutoff,
                         MAX_TEMPLATES, cif_dir, id_map_file='', start_idx=0, end_idx=0 ):

    # Prepare to collect output data
    output_labels = []
    output_allatom_labels = []

    if len(cif_dir) == 0:
        dir_name = os.path.dirname( os.path.abspath( sys.argv[0] ) )
        cif_dir = dir_name+'/PDB_RNA'

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
        template_coord_data = []
        for aln_line in aln_lines:
            if len(aln_line)!=9: continue # some kind of overflow in some alignments?

            query,template,eval,qstart,qend,tstart,tend,qaln,taln = aln_line

            if query != target: continue

            if int(qend)<int(qstart): continue # aligned to reverse complement!

            pdb_id,chain_id = template.split('_')

            # need to do alignment
            cif_path = os.path.join(cif_dir, f'{pdb_id.upper()}.cif.gz')
            if not os.path.isfile( cif_path ):
                cif_path = os.path.join(cif_dir, f'{pdb_id.lower()}.cif') # kaggle style
                if not os.path.isfile( cif_path ): continue # occasional alignment to DNA, ignore!

            release_date = release_dates[pdb_id.upper()] # pulled from PDB server

            if not skip_temporal_cutoff and is_before_or_on(temporal_cutoff,release_date): continue

            # these release dates in the CIF files can be buggy!
            title,release_date_unreliable = extract_title_release_date( cif_path )

            print('\n',target,temporal_cutoff,"   ",template)
            if title: print(f"PDB Title: {title}")
            if release_date: print(f"PDB Release Date: {release_date}")

            # sometimes there is a mismatch between PDB's fasta files and what's actually stored in coordinates,
            # so best to get the actual residue numbers for the chain
            chain_full_sequence,chain_sequence,chain_seq_nums,chain_ins_codes = extract_rna_sequence(cif_path,chain_id)

            # get 3d data
            alignment = []
            qstart=int(qstart)
            qend=int(qend)
            tstart=int(tstart)
            tend=int(tend)
            alignment.append( sequence[:(qstart-1)] + '-'*(tstart-1) + qaln + sequence[qend:]  )
            alignment.append( '-'*(qstart-1)        + 'X'*(tstart-1) + taln + '-'*(len(sequence)-qend) )
            print( alignment[0],'query' )
            print( alignment[1],'template' )
            chain_coord_data = get_coord_labels( cif_path, chain_id, chain_sequence, chain_seq_nums, chain_ins_codes )

            coord_data = get_target_coord_data( chain_coord_data, (alignment[1],alignment[0]) )

            # mismatch in FASTA sequence and the polyx info in the CIF file
            if len(coord_data) != len(sequence):
                print( 'WARNING! len(coord_data) != len(sequence)', 'len coord_data', len(coord_data), 'len sequence', len(sequence), 'qstart',qstart,'len qaln',len(qaln),'qend',qend)
                continue

            templates.append( template )
            template_coord_data.append( coord_data )

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
            output_allatom_label = deepcopy(output_label)

            # output templates, C1'
            for n in range(len(templates)):
                template = templates[n]
                res,resid,xyz,pdb_info = template_coord_data[n][i]
                assert( resid == i+1 )
                output_label[ f"x_{n+1}" ] = xyz[C1PRIME_KEY][0]
                output_label[ f"y_{n+1}" ] = xyz[C1PRIME_KEY][1]
                output_label[ f"z_{n+1}" ] = xyz[C1PRIME_KEY][2]

                for atom in ALL_ATOMS:
                    output_allatom_label.update( {
                        f"{atom}_x_{n+1}": xyz[atom][0],
                        f"{atom}_y_{n+1}": xyz[atom][1],
                        f"{atom}_z_{n+1}": xyz[atom][2]
                    })
                output_allatom_label.update( {f"pdb_id_{n+1}": template,f"pdb_seq_num_{n+1}": int(pdb_info[0]), f"pdb_ins_code_{n+1}": pdb_info[1], f"pdb_resname_{n+1}": pdb_info[2]} )

            # pad with blank models
            for n in range(len(templates),MAX_TEMPLATES):
                output_label[ f"x_{n+1}" ] = np.nan
                output_label[ f"y_{n+1}" ] = np.nan
                output_label[ f"z_{n+1}" ] = np.nan

                for atom in ALL_ATOMS:
                    output_allatom_label.update( {
                        f"{atom}_x_{n+1}": np.nan,
                        f"{atom}_y_{n+1}": np.nan,
                        f"{atom}_z_{n+1}": np.nan
                    })
                output_allatom_label.update( {f"pdb_id_{n+1}": "",f"pdb_seq_num_{n+1}": np.nan, f"pdb_ins_code_{n+1}": '', f"pdb_resname_{n+1}": ''} )


            output_labels.append( output_label )
            output_allatom_labels.append( output_allatom_label)

        num_targets += 1
        # if num_targets > 1: break # for debug!

    print(f'Completed {num_targets} targets\n')

    return output_labels, output_allatom_labels, targets

# Create a DataFrame and write to CSV
def output_csv( output_data, outfile ):
    df = pd.DataFrame(output_data)
    df.to_csv(outfile, index=False)
    print(f"Output written to {outfile}")


def output_template_labels_to_csv( output_labels, output_allatom_labels, targets, outdir='', outfile='', dataset_name='', start_idx=0, end_idx=0 ):
    assert( not( len(outfile)>0 and len(dataset_name)>0 ) )

    os.makedirs(outdir, exist_ok=True)
    if outdir[-1] != '/': outdir += '/'
    split_tag = ''
    if start_idx > 0:
        num_digits = len(str(len(targets)))
        split_tag = f'.{start_idx:0{num_digits}d}_{end_idx:0{num_digits}d}'

    if len( outfile ) == 0:
        if len( dataset_name) == 0: dataset_name = 'test'
        outfile = f"{outdir}{dataset_name}.templates{split_tag}.csv"
        outfile_allatom = f"{outdir}{dataset_name}.allatom_templates{split_tag}.csv"
    else:
        outfile = f"{outdir}/{outfile}"
        if outfile.count('labels.csv')>1: outfile_allatom = outfile.replace('labels.csv','allatom.csv')
        elif outfile.endswith('.csv'):    outfile_allatom = outfile.replace('.csv','.allatom.csv')
        else: outfile_allatom = outfile + '.allatom.csv'

    output_csv( output_labels, outfile )
    output_csv( output_allatom_labels, outfile_allatom )


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

    output_labels,output_allatom_labels,targets = get_template_labels( sequences_file, mmseqs_results_file, skip_temporal_cutoff,
                                                       MAX_TEMPLATES, cif_dir, id_map_file, start_idx, end_idx )


    outdir = args.outdir
    outfile = args.outfile
    dataset_name = args.dataset_name
    output_template_labels_to_csv( output_labels, output_allatom_labels, targets, outdir, outfile, dataset_name, start_idx, end_idx )