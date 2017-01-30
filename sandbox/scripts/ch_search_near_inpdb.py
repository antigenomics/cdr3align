import os
import numpy as np
import ch_get_superimpose_structures as ch_superimpose
import Bio.PDB
import argparse
from antigenomics_util import *
from geometry import center_of_mass

class Found(Exception): pass

tto = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H',
       'ILE':'I', 'LEU':'L', 'LYS':'L', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W',
       'TYR':'Y', 'VAL':'V'}

def get_pdbs_imposed(pdbs):
    pdbs_imposed = {}
    for antigen in sorted(pdbs):
        for idpdb1 in range(len(pdbs[antigen])):
            for idpdb2 in range(idpdb1 + 1, len((pdbs[antigen]))):
                pdb1 = sorted(pdbs[antigen])[idpdb1]
                pdb2 = sorted(pdbs[antigen])[idpdb2]
                if antigen in pdbs_imposed:
                    if pdb1 in pdbs_imposed[antigen]:
                        pdbs_imposed[antigen][pdb1][pdb2] = pdb2+"_align_"+pdb1+".pdb"
                    else:
                        pdbs_imposed[antigen][pdb1] = {pdb2:pdb2+"_align_"+pdb1+".pdb"}
                else:
                    pdbs_imposed[antigen] = {pdb1:{pdb2:pdb2+"_align_"+pdb1+".pdb"}}
    return pdbs_imposed

def getdistances(inppath, outpath, pdbs_imposed, pdbs, chainname='B'):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    try:
        with open(outpath, 'w') as out:
            out.write('pdb1\tpdb2\tpdb1_AA\tpdb2_AA\tdistance\tCDR_chain\tAA1\tAA2\t#1\t#2\tpdb1seq\tpdb2seq\n')
            for antigen in sorted(pdbs):                                            #select antigen 'antigen'
                for idpdb1 in range(len(pdbs[antigen])):                            #idpdb1 - id of pdbs
                    for idpdb2 in range(idpdb1 + 1, len((pdbs[antigen]))):          #idpdb2 - id of second pdb
                        pdb1 = sorted(pdbs[antigen])[idpdb1]
                        pdb2 = sorted(pdbs[antigen])[idpdb2]
                        first_structure = pdb_parser.get_structure("first", os.path.join(inppath, antigen, "pdb" + pdb1 + ".ent"))
                        second_structure = pdb_parser.get_structure("second", os.path.join(inppath, antigen, pdbs_imposed[antigen][pdb1][pdb2]))

                        first_model = first_structure[0]
                        second_model = second_structure[0]
                        for cdr_chain in pdbs[antigen][pdb1]:
                            if cdr_chain in pdbs[antigen][pdb2] and cdr_chain == chainname:
                                first_range = range(int(pdbs[antigen][pdb1][cdr_chain]['cdr_start']),int(pdbs[antigen][pdb1][cdr_chain]['cdr_end']))
                                sequence1 = pdbs[antigen][pdb1][cdr_chain]['CDR_seq']
                                second_range = range(int(pdbs[antigen][pdb2][cdr_chain]['cdr_start']),int(pdbs[antigen][pdb2][cdr_chain]['cdr_end']))
                                sequence2 = pdbs[antigen][pdb2][cdr_chain]['CDR_seq']
                                first_residues = get_residues(first_model[pdbs[antigen][pdb1][cdr_chain]['pdb_cdr']], first_range)
                                second_residues = get_residues(second_model[pdbs[antigen][pdb2][cdr_chain]['pdb_cdr']], second_range)
                                for i in first_residues:
                                    if i.get_id()[0] != ' ':
                                        first_residues.remove(i)
                                for i in second_residues:
                                    if i.get_id()[0] != ' ':
                                        second_residues.remove(i)

                                first_atoms = []
                                second_atoms = []

                                for first_res in first_residues:
                                    first_atoms.append([center_of_mass(first_res), first_res.get_id()[1], first_res.get_resname()])
                                for second_res in second_residues:
                                    second_atoms.append([center_of_mass(second_res), second_res.get_id()[1], second_res.get_resname()])
                                iti1 = 0
                                for k in first_atoms:
                                    iti2 = 0
                                    for i in second_atoms:
                                        aaa = np.array(k[0])
                                        bbb = np.array(i[0])
                                        distance = np.linalg.norm(aaa - bbb)
                                        out.write(pdb1 +
                                                  '\t' + pdb2+
                                                  '\t' + str(k[1]) +
                                                  '\t' + str(i[1]) +
                                                  '\t' + str(distance) +
                                                  '\t' + cdr_chain +
                                                  '\t' + tto[str(k[2])] +
                                                  '\t' + tto[str(i[2])] +
                                                  '\t' + str(iti1) +
                                                  '\t' + str(iti2) +
                                                  '\t' + sequence1 +
                                                  '\t' + sequence2 + '\n')
                                        iti2 += 1
                                    iti1 += 1
    except Found:
        print('')

