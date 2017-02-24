from Bio.SubsMat import MatrixInfo
import pandas as pd
import copy
import ch_antigen_base as ch_w
import glob
from shutil import copyfile
import random
import numpy as np
from Bio.PDB import *
from antigenomics_util import *
from datetime import datetime


def mindict(dic):
    i = None
    for obj in dic:
        if i == None:
            i = dic[obj]
        else:
            if i >= dic[obj]:
                i = dic[obj]
            for obj in dic:
                if i >= dic[obj]:
                    i = dic[obj]
    return i


# print(mindict(MatrixInfo.blosum30))
# aaa = 'test'
# mkdir(test)
# import numpy as np
# import matplotlib.pyplot as plt
# import glob
# from sklearn.neighbors import KernelDensity
#
# folders = glob.glob('group_pair_alignment/blosum*/')
# for folder in folders:
#    print(folder)
#    inps = glob.glob(folder+'trb_sequences*')
#    data = []
#    for i in inps:
#        inp = open(i,'r')
#        for line in inp:
#            if line.startswith('ins\tdel'):
#                for line in inp:
#                    if line.startswith('//'):
#                       break
#                    data.append(float(line.strip().split()[3]))
#       inp.close()
#    print('next_step_1')
#    X = np.array(data)[:, np.newaxis]
#    print('next_step_2')
#    X_plot = np.linspace(-10, 10, 2000)[:, np.newaxis]
#    fig, ax = plt.subplots()
#    print('next_step_3')
#    for kernel in ['gaussian', 'tophat', 'epanechnikov']:
#        kde = KernelDensity(kernel=kernel, bandwidth=0.2).fit(X)
#        log_dens = kde.score_samples(X_plot)
##        ax.plot(X_plot[:, 0], np.exp(log_dens), '-',
#               label="kernel = '{0}'".format(kernel))
#
#    ax.text(-1, -1, "N={0} points".format(len(data)))#
#
#    ax.legend(loc='upper left')
#    ax.plot(X[:, 0], -0.005 - 0.01 * np.random.random(X.shape[0]), '+k')#
#
#   ax.set_xlim(-2, 4)
#  ax.set_ylim(-2, 5)
# plt.savefig(folder.split('/')[1]+'.png')
#
aaa = pd.DataFrame(0, index=['A', 'B', 'C', 'D'], columns=['E', 'F', 'G'])
iti = 0
for i in aaa.columns:
    print(i)
    for k in aaa.index:
        print(k)
        iti += 1
        aaa[i][k] = iti
print(aaa)

# print(aaa.iloc[1:3].column)

# aaapandas = pd.DataFrame(aaascores, columns=['score', 'from', 'to', 'max/min'])
# print(aaapandas)
'''
aaamax = aaa.iloc[1:3].max(axis=1)
aaamin = aaa.iloc[1:3].min(axis=1)

print(aaamax)
for i in aaamax:
    print(i)

print(aaamin)
for i in aaamin:
    print(i)
print(aaa.index)

print('---')
print(aaa.iloc[1].idxmax())
len = 3
1,2
1,3
2,3
'''
'''
aaa = np.array([0.1*i for i in range(3)])
bbb = np.array([0.2*i for i in range(3)])
print(aaa, bbb)
dist = np.linalg.norm(aaa-bbb)
print(dist)

aaa = [1,2]
for i in range(len(aaa)-1):
    for k in range(i, len(aaa)):
        print(aaa[k])
'''
'''
input_file = 'final.annotations_3.txt'
table = pd.read_table(input_file)

bypdb = table.groupby("pdb_id")
print(bypdb)
for pdb_id, pdb_group in bypdb:
    if pdb_id == '1oga':
        print('\npdb_id:\n', pdb_id, '\npdb_group:\n', pdb_group)
        group = pdb_group
        break
pdb_id = '1oga'
pdb_parser = PDBParser()
pdb_list = PDBList()
pdb_file = pdb_list.retrieve_pdb_file(pdb_id, pdir='testdir')
print(pdb_id, "- preparing for computation")
model = pdb_parser.get_structure(pdb_id, pdb_file)[0]
pdb_annot = group.iloc[0]
print(pdb_annot)
antigen_chain = model[pdb_annot['chain_antigen']]
antigen_seq = pdb_annot['antigen_seq']
antigen_range = range(len(antigen_seq))
antigen_residues = get_residues(antigen_chain, antigen_range)
antigen_seq_obs = get_seq(antigen_residues)
print(antigen_residues)
print(antigen_seq_obs)
#antigen_seq_obs = get_seq(antigen_residues)
byregion = pdb_group.groupby(['chain_tcr', 'tcr_region'])
results_by_pdb = []
'''
'''for tcr_region_id, tcr_region_group in byregion:
    # Get and check tcr region residues
    tcr_annot = tcr_region_group.iloc[0]
    tcr_chain = model[tcr_annot['chain_tcr']]
    tcr_region_seq = tcr_annot['tcr_region_seq']
    tcr_region_range = range(tcr_annot['tcr_region_start'], tcr_annot['tcr_region_end'])
    tcr_region_residues = get_residues(tcr_chain, tcr_region_range)
    print(list(enumerate(tcr_chain)), tcr_chain)

    tcr_region_seq_obs = get_seq(tcr_region_residues)
    #print('-----\n', tcr_region_range, '\n', tcr_region_residues, '\n', tcr_region_seq_obs)

    if tcr_region_seq != tcr_region_seq_obs:
        warning("TCR:", tcr_region_id, " sequence mismatch (expected observed): ", tcr_region_seq,
                tcr_region_seq_obs, ". Replacing with one from PDB.")
        tcr_annot['tcr_region_seq'] = tcr_region_seq_obs
    #print(pdb_id, "- computing energies for", tcr_annot['tcr_v_allele'], ':', tcr_region_id[1])

    # Compute distances and add them to results
    distances = calc_distances(tcr_chain.get_id(), antigen_chain.get_id(),
                               tcr_annot['tcr_v_allele'], tcr_annot['tcr_region'],
                               tcr_region_residues, antigen_residues, tcr_region_range, antigen_range)
'''


# species, antigen_seq, pdb.
def fx(x):
    if isinstance(x, str):
        return (x[2])
    else:
        print(x)
        return x


'''
inptable = pd.read_csv('final.annotations_3.txt', sep='\t')
inptable_pdb = inptable.groupby('pdb_id')
inptable = inptable.dropna(axis=0)
outtable_species = inptable.sort_values(by=['species', 'antigen_seq', 'pdb_id'])
outtable_species['chain_allele']=outtable_species.tcr_v_allele.apply(fx)
outtable_species.to_csv('test.annotations_3.txt', sep='\t', columns=['species', 'antigen_seq', 'chain_allele',
                                                                     'tcr_region', 'pdb_id', 'tcr_region_seq', 'chain_tcr',
                                                                     'chain_antigen', 'tcr_region_start',
                                                                     'tcr_region_end'], index=False)
'''
'''
inptable = pd.read_csv('sorted.annotations_3.txt', sep='\t')
for species in inptable.groupby('species'):
    for tcr in species[1].groupby('tcr_region'):
        outinfo = tcr[1].sort_values(by=['antigen_seq', 'pdb_id'])
        outinfo.to_csv(os.path.join('testdir',str(species[0])+'_'+str(tcr[0])+'annotation.txt'), sep='\t',
                       columns=['species', 'antigen_seq', 'chain_allele', 'pdb_id', 'tcr_region_seq',
                                'chain_tcr', 'chain_antigen', 'tcr_region_start', 'tcr_region_end'],
                       index=False)
'''
'''
print('aaaa'
      'bbbb')
print('aaavvv'[::-1])

aaa = [1,2,3,4,5]
print(aaa.index(min(aaa)))

antigens = glob.glob('seqdata/HomoSapiens/pdbs/*')
for i in antigens:
    pdbs = glob.glob(i+'/*.ent')
    for k in pdbs:
        outpdb = k.split('/')[-1]
        copyfile(k, os.path.join('annotation_errors/pdb/', outpdb))

'''


aaa = {1:'a', 0:'b', 3:'c'}
print(aaa[min(aaa)])