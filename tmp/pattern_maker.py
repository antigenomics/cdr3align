

amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
amino_acids_volume = {'A':67, 'R':148, 'N':96, 'D':91, 'C':86, 'E':109, 'Q':114, 'G':48, 'H':118, 'I':124, 'L':124, 'K':135, 'M':124, 'F':135, 'P':90, 'S':73, 'T':93, 'W':163, 'Y':141, 'V':105}
amino_acids_hydrophobicity = {'A': 1, 'C': 1, 'D': -1, 'E': -1, 'F': 1, 'G': 1, 'H': -1, 'I': 1, 'K': -1, 'L': 1, 'M': 1, 'N': -1, 'P': 1, 'Q': -1, 'R': -1, 'S': -1, 'T': -1, 'V': 1, 'W': -1, 'Y': -1}
amino_acids_polarity = {'A': 0, 'C': 0, 'D': 1, 'E': 1, 'F': 1, 'G': 1, 'H': -1, 'I': 1, 'K': -1, 'L': 1, 'M': 1, 'N': -1, 'P': 1, 'Q': -1, 'R': -1, 'S': -1, 'T': -1, 'V': 1, 'W': -1, 'Y': -1}

import numpy as np
import matplotlib.pyplot as plt
import glob
from sklearn.neighbors import KernelDensity

'''
folders = glob.glob('group_pair_alignment/blosum*/')
for folder in folders:
    print(folder)
    inps = glob.glob(folder+'trb_sequences*')
    data = []
    for i in inps:
        inp = open(i,'r')
        for line in inp:
            if line.startswith('ins\tdel'):
                for line in inp:
                    if line.startswith('//'):
                       break
                    data.append(float(line.strip().split()[3]))
        inp.close()'''
inp = open('blast/trb_sequences_YPLHEQHGM_12','r')
seqdictvolume = {}
seqdicthydrophobicity = {}
seqdictpolarity = {}

for line in inp:
    if line.startswith('>'):
        pass
    else:
        sequence = line.strip()
        seqdictvolume[sequence] = []
        seqdicthydrophobicity[sequence] = []
        seqdictpolarity[sequence] = []
        for i in range(len(sequence)):
            seqdictvolume[sequence] += [i+1]*amino_acids_volume[sequence[i]]
            seqdicthydrophobicity[sequence] += [i+1]*int(str(amino_acids_hydrophobicity[sequence[i]]).replace('-1', '0'))
            seqdictpolarity[sequence] += [i+1]*amino_acids_polarity[sequence[i]]
        print(seqdictvolume[sequence])
        X_plot = np.linspace(-1, len(sequence), 100)[:, np.newaxis]
        fig, ax = plt.subplots()
        for dictt in [seqdicthydrophobicity, seqdictpolarity]:
            if dictt == seqdicthydrophobicity:
                name = 'seqdicthydrophobicity'
            elif dictt == seqdictpolarity:
                name = 'seqdictpolarity'
            X = np.array(dictt[sequence])[:, np.newaxis]
            #print('next_step_2')
            kernel = 'gaussian'
            kde = KernelDensity(kernel=kernel, bandwidth=0.5).fit(X)
            log_dens = kde.score_samples(X_plot)
            ax.plot(X_plot[:, 0], np.exp(log_dens), '-',
                    label="dictt = '{0}'".format(name))
            #print(np.exp(log_dens), log_dens)
        ax.legend(loc='upper left')
        #ax.plot(X[:, 0], -0.005 - 0.01 * np.random.random(X.shape[0]), '+k')#

        ax.set_xlim(-2, len(sequence)+3)
        ax.set_ylim(-1, 2)
        plt.savefig(sequence+'_volume.png')