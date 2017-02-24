import os
import pandas as pd
import ch_antigen_base as ch_w

class Found(Exception): pass


occurence = {'A': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
             'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

def read_scores(inppath, dataset = None):
    scores_datasets = {}
    if dataset is None:
        with open(inppath, 'r') as inp:
            for line in inp:
                for line in inp:
                    inf = line.strip().split()
                    if inf[2] in scores_datasets:
                        scores_datasets[inf[2]]['matrix'][inf[0]][inf[1]] = float(inf[3])
                    else:
                        scores_datasets[inf[2]] = {'matrix':pd.DataFrame(float(0), index=sorted(occurence), columns=sorted(occurence))}
                        scores_datasets[inf[2]]['matrix'][inf[0]][inf[1]] = float(inf[3])
    else:
        scores_datasets[dataset] = {'matrix':pd.DataFrame(float(0), index=sorted(occurence), columns=sorted(occurence))}
        with open(inppath, 'r') as inp:
            for line in inp:
                for line in inp:
                    inf = line.strip().split()
                    if inf[2] == dataset:
                        scores_datasets[inf[2]]['matrix'][inf[0]][inf[1]] = float(inf[3])
                    else:
                        pass
    for dset in scores_datasets:
        scores_datasets[dset]['maxscore'] = {}
        scores_datasets[dset]['minscore'] = {}
        for i in range(len(occurence)):
            scores_datasets[dset]['maxscore'][scores_datasets[dset]['matrix'].iloc[i].name] = {
                'score': scores_datasets[dset]['matrix'].iloc[i].max(),
                'to': scores_datasets[dset]['matrix'].iloc[i].idxmax(i)}
            scores_datasets[dset]['minscore'][scores_datasets[dset]['matrix'].iloc[i].name] = {
                'score': scores_datasets[dset]['matrix'].iloc[i].min(),
                'to': scores_datasets[dset]['matrix'].iloc[i].idxmin(i)}
    return scores_datasets

#same antigen, score, pair

def getmaxandmin_scores(scores):
    pass


def pair_alignment_to_scores(inppath, dataset, antig_seq, wrkdir, outname):
    ch_w.crdir(wrkdir)
    alignments = ch_w.get_alignments_var2(inppath)
    with open(os.path.join(wrkdir, outname), 'w') as out:
        out.write('sameantig\tscore\tpair\n')
        for pairset in alignments:
            p1, p2 = pairset.split('|')[0], pairset.split('|')[1]
            sameantig = 0
            try:
                for antig_seq1 in antig_seq[p1]:
                    if antig_seq1 in antig_seq[p2]:
                        raise Found
            except Found:
                sameantig = 1
            for pair in alignments[pairset]:
                iti = 0
                upscore = 0
                pair[0], pair[1] = pair[0].upper(), pair[1].upper()
                for i in range(len(pair[0])):
                    if pair[0][i] == 'C' or pair[1][i] == 'C':
                        pass
                    else:
                        iti += 1
                        upscore += dataset['matrix'][pair[0][i]][pair[1][i]]
                out.write(str(sameantig) + '\t' + str(upscore/iti) + '\t' + pair[0] + ' ' + pair[1] + '\n')
'''
sequences_with_antigen = ch_w.get_sequences_with_antigen('seqdata/HomoSapiens/trb_sequences_HomoSapiens_teach_antigen')
scores = read_scores('scores_modified.txt')

for dataset in scores:
    pair_alignment_to_scores('seqdata/HomoSapiens/sub_no_del_or_ins/out_teach_'+str(dataset),
                             scores[str(dataset)],
                             sequences_with_antigen,
                             'seqdata/HomoSapiens/scores',
                             'pairscores_'+str(dataset)+'.txt')
'''