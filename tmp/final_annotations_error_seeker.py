import os
import glob
import pandas as pd
import copy
import Bio.PDB
from antigenomics_util import *

tto = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def crdir(path):
    if os.path.exists(path):
        return 0
    else:
        os.mkdir(path)

outp = 'annotation_errors'
crdir(outp)

crdir(os.path.join(outp,'pdb'))

inpdata = pd.read_csv('../final.annotations.txt', sep='\t')
ok_pdbs = pd.unique(inpdata['pdb_id'])


def get_pdbs(inplist, outpath):
    pdbl = Bio.PDB.PDBList()
    iti = 0
    print(len(inplist), iti)
    for i in inplist:
        iti += 1
        if not os.path.isfile(os.path.join(outpath,'pdb'+i+'.ent')):
            print(len(inplist), iti)
            pdbl.retrieve_pdb_file(i, pdir=os.path.join(outpath))

bad_pdbs = {'chain_problem':[], 'seqpos_problem':[], 'bad_structure_problem':[], 'loss_of_data':[]}
possible_problem = {'chain_problem':[], 'seqpos_problem':[], 'bad_structure_problem':[], 'loss_of_data':[]}

def intersect(a, b):
    return(sorted(list(set(sorted(a))&set(sorted(b)))))

def is_in_list(a, b):
    if sorted(intersect(a, b)) == sorted(a):
        return True
    else:
        return False

def refine(a):
    b = a
    to_pop = set()
    if len(a) == 1:
        return b
    for i in range(1, len(a)+1):
        for j in range(i+1, len(a)+1):
            if is_in_list(a[str(i)], a[str(j)]):
                to_pop.add(str(i))
            elif is_in_list(a[str(j)], a[str(i)]):
                to_pop.add(str(j))
    for i in to_pop:
       b.pop(i)
    return b

pdb_models = {}
for i in ok_pdbs:
    with open(os.path.join("../annotation_errors/pdb/pdb"+i+".ent"), 'r') as inp:
        yep = 0
        pdbinfo = {}
        pdb_chains_interch = []
        for line in inp:
            if line.startswith('COMPND'):
                if line.strip().split()[2] == 'CHAIN:':
                    pdb_chains_interch.append(list(map(lambda x: x.replace(',','').replace(';',''), line.strip().split()[3:])))
            if line.startswith("REMARK 350"):
                yep = 1
                for line in inp:
                    if line.startswith('REMARK 350 BIOMOLECULE:'):
                        biomol = line.strip().split()[3]
                        for line in inp:
                            if line.startswith('REMARK 350 APPLY THE FOLLOWING TO CHAINS:'):
                                pdbinfo[biomol] = list(map(lambda x: x.replace(',','').replace(';',''), line.strip().split()[7:]))
                                if line.strip().split()[-1].endswith(','):
                                    for line in inp:
                                        pdbinfo[biomol] += list(
                                            map(lambda x: x.replace(',', '').replace(';', ''),
                                                line.strip().split()[4:]))
                                        if line.strip().split()[-1].endswith(','):
                                            pass
                                        else:
                                            break
                                break
                    elif line.startswith('REMARK 350'):
                        pass
                    else:
                        break
            elif line.startswith("SEQRES"):
                if yep == 0:
                    possible_problem['chain_problem'].append(i)
                    possible_problem['seqpos_problem'].append(i)
                    possible_problem['bad_structure_problem'].append(i)
                    print('Possible problem found: ',i,'\n')
                else:
                    pdbinfo2 = refine(pdbinfo)
                    pdb_models[i] = [pdbinfo2, pdb_chains_interch]
                break

pdb_grouped = inpdata.groupby('pdb_id')
for ipdb in pdb_grouped:


    #PART I. CHAIN PROBLEM

    tcr = list(map(str, (list(pd.unique(ipdb[1]['chain_tcr'])))))
    antigen = list(map(str, list(pd.unique(ipdb[1]['chain_antigen']))))
    chains_in_seq = list(map(str, (list(pd.unique(ipdb[1]['chain_tcr'])) + list(pd.unique(ipdb[1]['chain_antigen'])))))
    #chain in seq = chains in annotations.
    #pdb_models[ipdb[0]] - pdb file
    #pdb_models[ipdb[0]][0] - pdbinfo2 - chains in each model. [1] - pdb_chains_interch - chains that can be changed into each other
    if 'nan' in chains_in_seq:
        bad_pdbs['loss_of_data'].append(ipdb[0])
    else:
        is_ok = 0
        #bad_pdbs = {'chain_problem':[], 'seqpos_problem':[], 'bad_structure_problem':[], 'loss_of_data':[]}
        for i in pdb_models[ipdb[0]][0]:
            if intersect(pdb_models[ipdb[0]][0][i], chains_in_seq) == sorted(chains_in_seq):
                is_ok = 1
                if len(pdb_models[ipdb[0]][1]) != 5:
                    bad_pdbs['bad_structure_problem'].append([ipdb[0], 'Found '+str(len(pdb_models[ipdb[0]][1]))+' chains, while should be 5.'])
            elif intersect(pdb_models[ipdb[0]][0][i], chains_in_seq) == []:
                pass
            else:
                if intersect(tcr, pdb_models[ipdb[0]][0][i]) == []:
                    pass
                elif intersect(tcr, pdb_models[ipdb[0]][0][i]) == sorted(tcr):
                    is_ok = 2
                    for chainclass in pdb_models[ipdb[0]][1]:
                        if is_in_list(antigen, chainclass):
                            is_ok = 2
                            for ant in chainclass:
                                if ant in pdb_models[ipdb[0]][0][i]:
                                    is_ok = 2
                                    outant = ant
                                    break
                                else:
                                    is_ok = 0
                            break
                        else:
                            is_ok = 0
                    if is_ok == 2:
                        bad_pdbs['chain_problem'].append([ipdb[0], 'Change antigen from ' +
                                                          str(antigen).replace('[', '').replace(']', '') +
                                                          ' to \'' + str(outant)+'\''])
                        if len(pdb_models[ipdb[0]][1]) != 5:
                            bad_pdbs['bad_structure_problem'].append([ipdb[0], 'Found ' + str(len(pdb_models[ipdb[0]][1])) + ' chains, while should be 5.'])
        if is_ok == 0:
            bad_pdbs['chain_problem'].append([ipdb[0], 'Chaos in tcr chains and/or in antigen'])
            is_ok = 0

    #PART II. ANNOTATION PROBLEM

    pdb_parser = Bio.PDB.PDBParser(QUIET=True)

    work_structure = pdb_parser.get_structure("first", os.path.join(outp,'pdb', "pdb" + ipdb[0] + ".ent"))
    work_model = work_structure[0]

    for seq in ipdb[1].groupby('tcr_region_seq'):
        chain = list(pd.unique(seq[1]['chain_tcr']))[0]
        work_range = range(int(seq[1]['tcr_region_start']), int(seq[1]['tcr_region_end']))
        work_residues = get_residues(work_model[chain], work_range)
        work_seq = ''
        for i in work_residues:
            if i.get_id()[0] != ' ':
                work_residues.remove(i)
            else:
                work_seq += tto[str(i.get_resname())]
        if seq[0] != work_seq:
            if ipdb[0] not in bad_pdbs['seqpos_problem']:
                bad_pdbs['seqpos_problem'] = {ipdb[0]:[[seq[0], work_seq]]}
            else:
                bad_pdbs['seqpos_problem'][ipdb[0]].append([seq[0], work_seq])

with open(os.path.join(outp, 'problems.txt'), 'w') as out:
    for i in sorted(bad_pdbs):
        if i == 'chain_problem':
            out.writelines('Problems with chain:\n')
            for k in bad_pdbs[i]:
                out.writelines('\t'.join(k)+'\n')
            out.writelines('//\n\n')
        elif i == 'loss_of_data':
            out.writelines('Loss of data:\n'+'\n'.join(bad_pdbs[i])+'\n//\n\n')
        elif i == 'seqpos_problem':
            out.writelines('Problem with data annotation:\npdb\tshould_be\twhat_is\n')
            for k in bad_pdbs[i]:
                out.writelines('\t'.join(k)+'\n')
            out.writelines('//\n\n')
        elif i == 'bad_structure_problem':
            out.writelines('Problem with structures:\n')
            for k in bad_pdbs[i]:
                out.writelines('\t'.join(k)+'\n')
            out.writelines('//\n\n')
print(possible_problem)