import os
import glob
import pandas as pd
import csv
import json
import time
import copy
import math
import graphs_var1 as ch1_gr

files = ["out_teach_0_1_1_2",
         "out_teach_1_0_1_2",
         "out_teach_1_1_0_2"]

wrkdir = "seqdata/HomoSapiens"
inputdir = "seqdata/HomoSapiens/sub_del_ins_tot"
orgname = "HomoSapiens"

alignments = {}
for i in files:
    with open(os.path.join(inputdir, i), 'r') as inp:
        for line in inp:
            algn = line.strip().split()
            oseq, tseq = algn[0].replace('-', '').upper(), algn[1].replace('-', '').upper()
            if oseq + '|' + tseq in alignments:
                if [algn[0], algn[1]] in alignments[oseq + '|' + tseq]:
                    pass
                else:
                    alignments[oseq + '|' + tseq].append([algn[0], algn[1]])
            elif tseq + '|' + oseq in alignments:
                if [algn[1], algn[0]] in alignments[tseq + '|' + oseq]:
                    pass
                else:
                    alignments[tseq + '|' + oseq].append([algn[1], algn[0]])
            else:
                alignments[oseq + '|' + tseq] = [[algn[0], algn[1]]]

sequences_with_antigen = {}
with open(wrkdir+'/trb_sequences_'+orgname+'_teach_antigen', 'r') as inp:
    for line in inp:
        inf = line.strip().split()
        if inf[0] in sequences_with_antigen:
            sequences_with_antigen[inf[0]].append(inf[1])
        else:
            sequences_with_antigen[inf[0]] = [inf[1]]

#sequences, antigens = ch1_gr.get_seqs_antigens(sequences_with_antigen)
def get_gaps(alignments, indelsin, indelsout, indelsfrom, indelsto, indelsfromto, antig_seq, name, folder):
    windelsin = copy.deepcopy(indelsin)
    windelsout = copy.deepcopy(indelsout)
    windelsfromin = copy.deepcopy(indelsfrom)
    windelsfromout = copy.deepcopy(indelsfrom)
    windelstoin = copy.deepcopy(indelsto)
    windelssin = copy.deepcopy(indelsto)
    windelstoout = copy.deepcopy(indelsto)
    windelssout = copy.deepcopy(indelsto)
    windelsfromtoin = copy.deepcopy(indelsfromto)
    windelsfromtoout = copy.deepcopy(indelsfromto)
    def getinfo(pair, p1, p2, windels, windelsfrom, windelsto, windelsfromto, windelss):
        if len(p1) in windels:
            pass
        else:
            windels[len(p1)] = [0 for i in range(len(p1))]

        if len(p2) in windels:
            pass
        else:
            windels[len(p2)] = [0 for i in range(len(p2))]

        m1 = 0
        m2 = 0
        for m in range(len(pair[0])):
            a1, a2 = pair[0][m].upper(), pair[1][m].upper()
            if a1 == '-':
                windels[len(p2)][m2] += 1
                if (m1 != m or m1 != 0) and m1 < len(p1):
                    windelsfrom[p1[m1]] += 1
                if m1+1 < len(p1):
                    windelsto[p1[m1+1]] += 1
                if (m1 != m or m2 != 0) and m1+1 < len(p1):
                    windelsfromto[p1[m1]][p1[m1+1]] += 1
                windelss[p2[m2]] += 1
                m2 += 1
            elif a2 == '-':
                windels[len(p1)][m1] += 1
                if (m2 != m or m2 != 0) and m2 < len(p2):
                    windelsfrom[p2[m2]] += 1
                if m2 + 1 < len(p2):
                    windelsto[p2[m2 + 1]] += 1
                if (m2 != m or m2 != 0) and m2 + 1 < len(p2):
                    windelsfromto[p2[m2]][p2[m2+1]] += 1
                windelss[p1[m1]] += 1
                m1 += 1
            else:
                m1 += 1
                m2 += 1
        return windels, windelsfrom, windelsfromto, windelss

    for pairset in alignments:
        p1, p2 = pairset.split('|')[0], pairset.split('|')[1]
        sameantig = 0
        for antig_seq1 in antig_seq[p1]:
            for antig_seq2 in antig_seq[p2]:
                if antig_seq1 == antig_seq2:
                    sameantig = 1
        for pair in alignments[pairset]:
            if len(pair[0]) != len(pair[1]):
                print('different lens', pair[0], pair[1])
                print('ERRROROR')
                return 0
            elif sameantig == 1:
                windelsin, windelsfromin, windelsfromtoin, windelssin = getinfo(pair, p1, p2, windelsin, windelsfromin,
                                                                    windelstoin, windelsfromtoin, windelssin)
            else:
                windelsout, windelsfromout, windelsfromtoout, windelssout = getinfo(pair, p1, p2, windelsout, windelsfromout,
                                                                       windelstoout, windelsfromtoout, windelssout)

    with open(os.path.join(folder, name)+'_indels_all.txt', 'w') as out:
        #maxlen = (sorted(windelsin)[-1])
        #out.write('length:\t' + '\t'.join(map(str, [i for i in range(maxlen)]))+'\n')
        #for i in sorted(windelsin):
        #    out.write(str(i)+':\t'+'\t'.join(map(str, windelsin[i]))+'\n')
        out.write('len\tpos\tindels\ttype\n')
        for i in sorted(windelsin):
            for k in range(len(windelsin[i])):
                out.write(str(i) + '\t' + str(k) + '\t' + str(windelsin[i][k]) + '\tinner\n')
        for i in sorted(windelsout):
            for k in range(len(windelsout[i])):
                out.write(str(i) + '\t' + str(k) + '\t' + str(windelsout[i][k]) + '\touter\n')

    with open(os.path.join(folder, name) + '_indels_ft', 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['Amino', 'to_gap_in', 'to_gap_out', 'from_in','to_in', 'from_out', 'to_out'])
        for i in sorted(windelsfromin):
            writer.writerow([i, windelssin[i], windelssout[i], windelsfromin[i], windelstoin[i], windelsfromout[i], windelstoout[i]])

    with open(os.path.join(folder, name) + '_indels_fromto_all.txt', 'w') as out:
        out.write('from\tto\tcount\ttype\n')
        for i in indelsfrom:
            for k in indelsto:
                out.write(str(i)+'\t'+str(k)+'\t'+str(windelsfromtoin[i][k])+'\tinner\n')
        for i in indelsfrom:
            for k in indelsto:
                out.write(str(i) + '\t' + str(k) + '\t' + str(windelsfromtoout[i][k]) + '\touter\n')

    #windelsfromtoin.to_csv(os.path.join(folder, name)+'_fromto_in.csv')
    #windelsfromtoout.to_csv(os.path.join(folder, name) + '_fromto_out.csv')
indelsin = {}
indelsout = {}
indelsfrom = {}
indelsto = {}
amino_acids = {}

with open('aa_property_table.txt','r') as inp:
    for line in inp:
        if line.startswith('amino_acid'):
            properties = line.strip().split()[1:len(line.strip().split())-1]
            for line in inp:
                amino = line.strip().split()
                amino_acids[amino[0]]={}
                indelsto[amino[0]] = 0
                indelsfrom[amino[0]] = 0
                for i in range(len(properties)):
                    amino_acids[amino[0]].update({properties[i]: amino[i+1]})

indelsfromto = pd.DataFrame(0, index=indelsfrom, columns=indelsto)

ch1_gr.crdir(os.path.join(wrkdir, 'indels'))
get_gaps(alignments, indelsin, indelsout, indelsfrom, indelsto, indelsfromto, sequences_with_antigen, 'indels', os.path.join(wrkdir,'indels'))
#It would be interesting to somehow manage to concatenate position of indels with fromto..

