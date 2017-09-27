from ch_base import *
from ch_read_blosum import *
from ch_imgt import *
import sys

input_file = '../../imgt_work/imgt_all.fasta'

dictimgt = get_specific_vdj(input_file, 'V') #organism, gene, seq

all_merged_loopseqs = {}
all_loopseq_representative = {}

to_remove = {}
for org in dictimgt:
    to_remove[org] = set()
    for gene in dictimgt[org]:
        dictimgt[org][gene] = ' '.join(list(str(dictimgt[org][gene][anchors[i][0]:anchors[i][1]]) for i in sorted(anchors)))

        if (gene.startswith('TRB') and len(dictimgt[org][gene]) != 30) or (gene.startswith('TRA') and len(dictimgt[org][gene]) !=30):
            #print(gene, dictimgt[org][gene])
            to_remove[org].add(gene)

    for gene in to_remove[org]:
        del dictimgt[org][gene]

    all_merged_loopseqs[org] = {}
    all_loopseq_representative[org] = {}
    all_loopseq_nbrs = {}
    for gene1 in dictimgt[org]:
        all_loopseq_nbrs[gene1] = set()
        for gene2 in dictimgt[org]:
            if dictimgt[org][gene1] == dictimgt[org][gene2]:
                all_loopseq_nbrs[gene1].add(gene2)
        all_loopseq_representative[org][gene1] = min(all_loopseq_nbrs[gene1])
        all_merged_loopseqs[org][min(all_loopseq_nbrs[gene1])] = dictimgt[org][gene1]
