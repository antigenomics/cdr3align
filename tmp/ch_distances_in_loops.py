import os
import numpy as np
import pandas as pd
import Bio.PDB
import argparse
from antigenomics_util import *
from geometry import center_of_mass


tto = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


pdb_parser = Bio.PDB.PDBParser(QUIET=True)
inppath = 'seqdata/HomoSapiens/cdrs/cdr3.annotation.txt'
outpath = 'seqdata/HomoSapiens/cdrs/cdr3.length_in_seq.txt'
inptable = pd.read_csv(inppath, sep='\t')
inptable['tcr_length'] = inptable['tcr_region_seq'].str.len()
inptable['loop_distance'] = 'NA'
inptable['loop_distance_var1'] = 'NA' #2 -2 for B, 1 -1 for A
inptable['loop_distance_var2'] = 'NA' #3 -3 for B, 2 -2 for A
for pdb in inptable.groupby('pdb_id'):
    work_structure = pdb_parser.get_structure("first", os.path.join('seqdata/HomoSapiens/allpdb', "pdb" + pdb[0] + ".ent"))
    work_model = work_structure[0]
    #for chain in list(pd.unique(pdb[1]['chain_tcr'])):
    for allele in pdb[1].groupby('chain_allele'):
        #print(inptable.iloc[(allele[1].index[0]), -1])
        chain = list(pd.unique(allele[1]['chain_tcr']))[0]
        work_range = range(int(allele[1]['tcr_region_start']), int(allele[1]['tcr_region_end']))
        work_residues = get_residues(work_model[chain], work_range)
        for i in work_residues:
            if i.get_id()[0] != ' ':
                work_residues.remove(i)
        work_atoms = []
        for work_res in work_residues:
            work_atoms.append([center_of_mass(work_res), work_res.get_id()[1], work_res.get_resname()])

        start, end = np.array(work_atoms[0][0]), np.array(work_atoms[-1][0])
        distance = np.linalg.norm(start - end)
        inptable.iloc[(allele[1].index[0]), -3] = distance

        start, end = np.array(work_atoms[2][0]), np.array(work_atoms[-2][0])
        distance = np.linalg.norm(start - end)
        inptable.iloc[(allele[1].index[0]), -2] = distance

        start, end = np.array(work_atoms[3][0]), np.array(work_atoms[-3][0])
        distance = np.linalg.norm(start - end)
        inptable.iloc[(allele[1].index[0]), -1] = distance

inptable.to_csv(outpath, sep='\t')