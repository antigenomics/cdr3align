import os
import numpy as np
import pandas as pd
import Bio.PDB
import argparse
from antigenomics_util import *
from geometry import center_of_mass
import ch_antigen_base as ch_w

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will get distances between CDR and antigen\
Current version is %s" % curr_version)

parser.add_argument("-i", help="Path to cdr3 annotation file. Default is '../seqdata/HomoSapiens/cdrs/cdr3.annotation.txt'",
                    required=False, dest="input_dir")
parser.add_argument("-o", help="Path to output directory. Default is '../seqdata/HomoSapiens/distances/'", required=False, dest='output_dir')
parser.add_argument("-pdb", help="path to PDB files. Default is ../seqdata/HomoSapiens/allpdb", required=False, dest="pdb")

myargs = parser.parse_args()

if myargs.input_dir is None:
    myargs.input_dir =  '../seqdata/HomoSapiens/cdrs/cdr3.annotation.txt'
if myargs.output_dir is None:
    myargs.output_dir = '../seqdata/HomoSapiens/distances/'
if myargs.pdb is None:
    myargs.pdb = '../seqdata/HomoSapiens/allpdb/'
ch_w.crdir(myargs.output_dir)

class Found(Exception): pass

tto = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def refine_seq(seq):
    for i in seq:
        if i.get_id()[0] != ' ':
            seq.remove(i)
    return(seq)

def list_atoms(seq):
    list_atoms = []
    for work_res in seq:
        list_atoms.append([work_res['CA'].get_vector(), work_res.get_id()[1], work_res.get_resname()])
    return list_atoms

def get_min_distance(amino, ant):
    distances = {}
    for k in range(len(ant)):
        distance = (amino[0] - ant[k][0]).norm()
        distances[distance] = [amino, ant[k], k]
    return(distances[min(distances)], min(distances))

def cdrantigen_distances(inppath = 'seqdata/HomoSapiens/cdrs/cdr3.annotation.txt',
                         outpath = 'seqdata/HomoSapiens/distances/cdr_and_antigen.txt',
                         inppdbs = 'seqdata/HomoSapiens/allpdb'):
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    inptable = pd.read_csv(inppath, sep='\t')
    inptable['tcr_length'] = inptable['tcr_region_seq'].str.len()
    with open(outpath, 'w') as out:
        out.write('pdb\tdistance\tCDR_chain\tAA_cdr\tAA_ant\tposcdr\tposant\tcdr_seq\tant_seq\tascordesc\n')
        def write_lines(pdb, cdrinfo, antinfo, distance, allele, citi, aiti, cdrseq, antseq, ascordesc):
            out.write('\t'.join(list(map(str,[pdb,
                                distance,
                                allele,
                                tto[cdrinfo[2]],
                                tto[antinfo[2]],
                                citi,
                                aiti,
                                cdrseq,
                                antseq,
                                ascordesc])))+'\n')

        for pdb in inptable.groupby('pdb_id'):
            work_structure = pdb_parser.get_structure("first",
                                                      os.path.join(inppdbs, "pdb" + pdb[0] + ".ent"))
            work_model = work_structure[0]
            for allele in pdb[1].groupby('chain_allele'):
                chain = list(pd.unique(allele[1]['chain_tcr']))[0]
                cdr_seq = list(pd.unique(allele[1]['tcr_region_seq']))[0]
                work_range = range(int(allele[1]['tcr_region_start']), int(allele[1]['tcr_region_end']))
                work_residues = get_residues(work_model[chain], work_range)
                work_residues = refine_seq(work_residues)
                work_atoms = list_atoms(work_residues)

                antigen = list(pd.unique(allele[1]['chain_antigen']))[0]
                antigen_seq = list(pd.unique(allele[1]['antigen_seq']))[0]
                antigen_range = range(int(0), int(len(antigen_seq)))
                antigen_residues = get_residues(work_model[antigen], antigen_range)
                antigen_residues = refine_seq(antigen_residues)
                antigen_atoms = list_atoms(antigen_residues)

                way = 'asc'
                iti = 0
                if len(work_atoms)%2 == 0:
                    lim = (len(work_atoms)//2) - 1
                    for i in range(len(work_atoms)):
                        atom, dist = get_min_distance(work_atoms[i], antigen_atoms)
                        write_lines(pdb[0], atom[0], atom[1], dist, allele[0], iti, atom[2], cdr_seq, antigen_seq, way)
                        if i < lim:
                            iti += 1
                        elif i == lim:
                            iti = iti
                            way = 'desc'
                        else:
                            iti -= 1
                else:
                    lim = (len(work_atoms)//2)
                    for i in range(len(work_atoms)):
                        atom, dist = get_min_distance(work_atoms[i], antigen_atoms)
                        write_lines(pdb[0], atom[0], atom[1], dist, allele[0], iti, atom[2], cdr_seq, antigen_seq, way)
                        if i < lim:
                            iti += 1
                        elif i == lim:
                            way = 'desc'
                            write_lines(pdb[0], atom[0], atom[1], dist, allele[0], iti, atom[2], cdr_seq, antigen_seq, way)
                            iti -= 1
                        else:
                            iti -= 1

def cdrantigen_min_distances(inppath='seqdata/HomoSapiens/distances/cdr_and_antigen.txt', outpath='seqdata/HomoSapiens/distances/min_cdr_and_antigen.txt'):
    inptable = pd.read_csv(inppath, sep='\t')
    inptable['tcr_length'] = inptable['cdr_seq'].str.len()
    pdbs = inptable.groupby('pdb')
    with open(outpath, 'w') as out:
        out.writelines('pdb\tdistance\tCDR_chain\tAA_cdr\tAA_ant\tposcdr\tposant\tcdr_seq\tant_seq\tascordesc\tcdr_len\n')
        for pdb in pdbs:
            for cdr in pdb[1].groupby('CDR_chain'):
                indx = inptable.iloc[cdr[1]['distance'].idxmin()]
                out.writelines('\t'.join(list(map(str, indx.iloc[:].values)))+'\n')

cdrantigen_distances(myargs.input_dir, os.path.join(myargs.output_dir, 'cdr_and_antigen.txt'), myargs.pdb)
cdrantigen_min_distances(os.path.join(myargs.output_dir, 'cdr_and_antigen.txt'), os.path.join(myargs.output_dir, 'min_cdr_and_antigen.txt'))