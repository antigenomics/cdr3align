import os
import glob
import pandas as pd
import csv
import json
import time
import copy
import math
import graphs_var1 as ch1
from datetime import datetime


def get_occurence_matrix(occurence, alignments, substitutions, work_occur=None):
    if work_occur != None:
        for i in alignments:
            loc_occurence = copy.deepcopy(occurence)
            localsub = pd.DataFrame(0, index=occurence, columns=occurence)
            for k in alignments[i]:
                if len(k[0]) != len(k[1]):
                    print('different lens', k[0], k[1])
                    print('ERRROROR')
                    return 0
                else:
                    for m in range(len(k[0])):
                        a1 = k[0][m].upper()
                        a2 = k[1][m].upper()
                        if a1 != '-' and a2 != '-':
                            if a1 != a2:
                                localsub[a1][a2] += 1
                                loc_occurence[a1] += 1
                                loc_occurence[a2] += 1
            substitutions += localsub / len(alignments[i])
            for k in work_occur:
                work_occur[k] += loc_occurence[k]/len(alignments[i])
        substitutions.to_csv('substitution_matrix_3_0_0_3.csv')
        out = open('aminoacid_occurence_3_0_0_3','w')
        for i in work_occur:
            out.write(str(i)+'\t'+str(work_occur[i])+'\n')
        return substitutions, work_occur
    else:
        for i in alignments:
            localsub = pd.DataFrame(0, index=occurence, columns=occurence)
            for k in alignments[i]:
                if len(k[0]) != len(k[1]):
                    print('different lens', k[0], k[1])
                    print('ERRROROR')
                    return 0
                else:
                    for m in range(len(k[0])):
                        a1 = k[0][m].upper()
                        a2 = k[1][m].upper()
                        if a1 != '-' and a2 != '-' and a1 != a2:
                            localsub[a1][a2] += 1
            substitutions+=localsub/len(alignments[i])
        substitutions.to_csv('substitution_matrix_3_0_0_3.csv')
        return substitutions

def get_occurence_matrix_default(occurence, alignments, substitutions, work_occur=None):
    if work_occur != None:
        for i in alignments:
            loc_occurence = copy.deepcopy(occurence)
            localsub = pd.DataFrame(0, index=occurence, columns=occurence)     #DEFAULT
            for k in alignments[i]:
                if len(k[0]) != len(k[1]):
                    print('different lens', k[0], k[1])
                    print('ERRROROR')
                    return 0
                else:
                    for m in range(len(k[0])):
                        a1 = k[0][m].upper()
                        a2 = k[1][m].upper()
                        if a1 != '-' and a2 != '-':
                            localsub[a1][a2] += 1
                            if a1 != a2:
                                loc_occurence[a1] += 1
                                loc_occurence[a2] += 1
            substitutions += localsub / len(alignments[i])                     #DEFAULT
            for k in work_occur:
                work_occur[k] += loc_occurence[k] / len(alignments[i])
        substitutions.to_csv('substitution_matrix_3_0_0_3_all.csv')
        out = open('aminoacid_occurence_3_0_0_3_all', 'w')
        for i in work_occur:
            out.write(str(i) + '\t' + str(work_occur[i]) + '\n')
        out.close()
        return substitutions, work_occur
    else:
        for i in alignments:
            localsub = pd.DataFrame(0, index=occurence, columns=occurence)
            for k in alignments[i]:
                if len(k[0]) != len(k[1]):
                    print('different lens', k[0], k[1])
                    print('ERRROROR')
                    return 0
                else:
                    for m in range(len(k[0])):
                        a1 = k[0][m].upper()                                    #DEFAULT
                        a2 = k[1][m].upper()
                        if a1 != '-' and a2 != '-':
                            localsub[a1][a2] += 1
            substitutions += localsub / len(alignments[i])
        substitutions.to_csv('substitution_matrix_3_0_0_3_all.csv')
        return substitutions

def get_occurence_matrix_2(occurence, alignments, substitutions, msubstitutions, antig_seq, name, folder, work_occur):
    mwork_occur = copy.deepcopy(work_occur)
    def getinfo(gisub, gioccurence):
        for m in range(len(pair[0])):
            a1, a2 = pair[0][m].upper(), pair[1][m].upper()
            if a1 != '-' and a2 != '-':
                # if a1 != a2:
                gisub[a1][a2] += 1
                gioccurence[a1] += 1
                gioccurence[a2] += 1
    for pairset in alignments:
        loc_occurence, loc_moccurence = copy.deepcopy(occurence), copy.deepcopy(occurence)
        localsub, localmsub = pd.DataFrame(0, index=occurence, columns=occurence), \
                              pd.DataFrame(0, index=occurence, columns=occurence)
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
                getinfo(localsub, loc_occurence)
            else:
                getinfo(localmsub, loc_moccurence)
        substitutions += localsub / len(alignments[pairset])
        msubstitutions += localmsub / len(alignments[pairset])
        for k in work_occur:
            work_occur[k] += loc_occurence[k]/len(alignments[pairset])
        for k in mwork_occur:
            mwork_occur[k] += loc_moccurence[k]/len(alignments[pairset])
    print(folder, name)
    substitutions.to_csv(folder+'substitution_matrix_'+name+'_2.csv')
    msubstitutions.to_csv(folder+'msubstitution_matrix_' + name + '_2.csv')
    out = open(folder+'aminoacid_occurence_'+name+'_2','w')
    for i in work_occur:
        out.write(str(i) + '\t' + str(work_occur[i]) + '\n')
    out.close()
    out = open(folder+'aminoacid_moccurence_' + name+'_2', 'w')
    for i in mwork_occur:
        out.write(str(i) + '\t' + str(mwork_occur[i]) + '\n')
    out.close()
    return substitutions, msubstitutions, work_occur, mwork_occur

def get_occurence_matrix_diagonal(occurence, alignments, substitutions, msubstitutions, antig_seq, name, folder, work_occur):
    mwork_occur = copy.deepcopy(work_occur)
    for i in alignments:
        loc_occurence, loc_moccurence = copy.deepcopy(occurence), copy.deepcopy(occurence)
        localsub, localmsub = pd.DataFrame(0, index=occurence, columns=occurence), pd.DataFrame(0, index=occurence, columns=occurence)      #VAR_2
        for k in alignments[i]:
            if len(k[0]) != len(k[1]):
                print('different lens', k[0], k[1])
                print('ERRROROR')
                return 0
            elif antig_seq[k[0].replace('-','').upper()] == antig_seq[k[1].replace('-','').upper()]:
                for m in range(len(k[0])):
                    a1, a2 = k[0][m].upper(), k[1][m].upper()
                    if a1 != '-' and a2 != '-':
                        if a1 == a2:
                            localsub[a1][a2] += 1
                            loc_occurence[a1] += 1
                            loc_occurence[a2] += 1
            else:
                if len(antig_seq[k[0].replace('-','').upper()]) > 1 or len(antig_seq[k[1].replace('-','').upper()]) > 1:
                    sameant = 0
                    for antig_seq1 in antig_seq[k[0].replace('-','').upper()]:
                        for antig_seq2 in antig_seq[k[1].replace('-','').upper()]:
                            if antig_seq1 == antig_seq2:
                                sameant = 1
                    if sameant == 0:
                        for m in range(len(k[0])):
                            a1, a2 = k[0][m].upper(), k[1][m].upper()
                            if a1 != '-' and a2 != '-':
                                if a1 == a2:
                                    localmsub[a1][a2] += 1
                                    loc_moccurence[a1] += 1
                                    loc_moccurence[a2] += 1
                    else:
                        for m in range(len(k[0])):
                            a1, a2 = k[0][m].upper(), k[1][m].upper()
                            if a1 != '-' and a2 != '-':
                                if a1 == a2:
                                    localsub[a1][a2] += 1
                                    loc_occurence[a1] += 1
                                    loc_occurence[a2] += 1
                else:
                    for m in range(len(k[0])):
                        a1, a2 = k[0][m].upper(), k[1][m].upper()
                        if a1 != '-' and a2 != '-':
                            if a1 == a2:
                                localmsub[a1][a2] += 1
                                loc_moccurence[a1] += 1
                                loc_moccurence[a2] += 1
        substitutions += localsub / len(alignments[i])
        msubstitutions += localmsub / len(alignments[i])
        for k in work_occur:
            work_occur[k] += loc_occurence[k]/len(alignments[i])
        for k in mwork_occur:
            mwork_occur[k] += loc_moccurence[k] / len(alignments[i])
    substitutions.to_csv(folder+'substitution_matrix_'+name+'.csv')
    msubstitutions.to_csv(folder+'msubstitution_matrix_' + name + '.csv')
    out = open(folder+'aminoacid_occurence_'+name,'w')
    for i in work_occur:
        out.write(str(i) + '\t' + str(work_occur[i]) + '\n')
    out.close()
    out = open(folder+'aminoacid_moccurence_' + name, 'w')
    for i in mwork_occur:
        out.write(str(i) + '\t' + str(mwork_occur[i]) + '\n')
    out.close()
    return substitutions, msubstitutions, work_occur, mwork_occur

def get_occurence_matrix_3(occurence, alignments, substitutions, msubstitutions, antig_seq, name, folder, work_occur):
    mwork_occur = copy.deepcopy(work_occur)
    def getinfo(gisub, gioccurence):
        for m in range(len(pair[0])):
            a1, a2 = pair[0][m].upper(), pair[1][m].upper()
            if a1 != '-' and a2 != '-':
                # if a1 != a2:
                gisub[a1][a2] += 1
                gioccurence[a1] += 1
                gioccurence[a2] += 1
        return(gisub, gioccurence)
    for pairset in alignments:
        loc_occurence, loc_moccurence = copy.deepcopy(occurence), copy.deepcopy(occurence)
        localsub, localmsub = pd.DataFrame(0, index=occurence, columns=occurence), \
                              pd.DataFrame(0, index=occurence, columns=occurence)
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
                localsub, loc_occurence = getinfo(localsub, loc_occurence)
            else:
                localmsub, loc_moccurence = getinfo(localmsub, loc_moccurence)
        substitutions += localsub / len(alignments[pairset])
        msubstitutions += localmsub / len(alignments[pairset])
        for k in work_occur:
            work_occur[k] += loc_occurence[k]/len(alignments[pairset])
        for k in mwork_occur:
            mwork_occur[k] += loc_moccurence[k]/len(alignments[pairset])
    substitutions.to_csv(folder+'substitution_matrix_'+name+'_3.csv')
    msubstitutions.to_csv(folder+'msubstitution_matrix_' + name + '_3.csv')
    out = open(folder+'aminoacid_occurence_'+name+'_3','w')
    for i in work_occur:
        out.write(str(i) + '\t' + str(work_occur[i]) + '\n')
    out.close()
    out = open(folder+'aminoacid_moccurence_' + name+'_3', 'w')
    for i in mwork_occur:
        out.write(str(i) + '\t' + str(mwork_occur[i]) + '\n')
    out.close()
    return substitutions, msubstitutions, work_occur, mwork_occur

def get_occurence(inpp, work_occur, occurence):
    seqs = set()
    with open(inpp, 'r') as inp:
        for line in inp:
            seqs.add(line.strip())
    for i in seqs:
        for k in i:
            work_occur[k.upper()]+=1
    sum_occur = 0
    for i in work_occur:
        sum_occur+=work_occur[i]
    log_sum_occur = math.log2(sum_occur)
    log_prob_occurence = copy.deepcopy(occurence)
    for i in log_prob_occurence:
        log_prob_occurence[i] = math.log2(work_occur[i])-log_sum_occur
    return log_prob_occurence

def get_log_prob_occur(work_occur, occurence):
    sum_occur = 0
    for i in work_occur:
        sum_occur+=work_occur[i]
    log_sum_occur = math.log2(sum_occur)
    log_prob_occurence = copy.deepcopy(occurence)
    for i in log_prob_occurence:
        log_prob_occurence[i] = math.log2(work_occur[i])-log_sum_occur
    return log_prob_occurence


#sequences_with_antigen, sequences_list = ch1.dict_sequences_with_antigen('trb_sequences_all_teach_antigen', var = 2)
#print(sequences_with_antigen)

alignments = {}
infname = 'out_teach_2_0_0_2'
with open('seqdata/HomoSapiens/sub_del_ins_tot/'+infname, 'r') as inp:
    for line in inp:
        algn = line.strip().split()
        oseq, tseq = algn[0].replace('-','').upper(), algn[1].replace('-','').upper()
        if oseq+'|'+tseq in alignments:
            alignments[oseq+'|'+tseq].append([algn[0], algn[1]])
        else:
            alignments[oseq + '|' + tseq]=[[algn[0], algn[1]]]

amino_acids = {}
occurence = {}
with open('aa_property_table.txt','r') as inp:
    for line in inp:
        if line.startswith('amino_acid'):
            properties = line.strip().split()[1:len(line.strip().split())-1]
            for line in inp:
                amino = line.strip().split()
                amino_acids[amino[0]]={}
                occurence[amino[0]] = 0
                for i in range(len(properties)): amino_acids[amino[0]].update({properties[i]: amino[i+1]})


orgname = "HomoSapiens"
wrkfolder = 'seqdata/'+orgname+'/'
sequences_with_antigen = {}
with open('seqdata/'+orgname+'/trb_sequences_'+orgname+'_teach_antigen', 'r') as inp:
    for line in inp:
        inf = line.strip().split()
        if inf[0] in sequences_with_antigen:
            sequences_with_antigen[inf[0]].append(inf[1])
        else:
            sequences_with_antigen[inf[0]] = [inf[1]]

working_occurence = copy.deepcopy(occurence)
substitutions = pd.DataFrame(0, index=sorted(occurence), columns=sorted(occurence))
msubstitutions = pd.DataFrame(0, index=sorted(occurence), columns=sorted(occurence))

print(substitutions, msubstitutions)
#substitutions, msubstitutions, working_occurence, mworking_occurence = get_occurence_matrix_2(occurence, alignments, substitutions, msubstitutions, antig_seq=sequences_with_antigen, name=infname, folder=wrkfolder, work_occur=working_occurence)

#substitutions = pd.read_csv('substitution_matrix_3_0_0_3_withm.csv', index_col=0)
#msubstitutions = pd.read_csv('msubstitution_matrix_3_0_0_3_withm.csv', index_col=0)

def other(substitutions, msubstitutions, occurence):
    sum_sub = substitutions.sum(axis=1).sum(axis=0)
    sum_msub = msubstitutions.sum(axis=1).sum(axis=0)

    sub_occur = pd.DataFrame(0, index=occurence, columns=occurence)
    msub_occur = pd.DataFrame(0, index=occurence, columns=occurence)

    sub_on_msub_norm = pd.DataFrame(0, index=occurence, columns=occurence)
    sub_on_msub = pd.DataFrame(0, index=occurence, columns=occurence)

    for i in occurence:
        for k in occurence:
            if substitutions[i][k] != 0:
                if i == k:
                    sub_occur[i][k] = 0
                else:
                    sub_occur[i][k] = float(substitutions[i][k])/float(sum_sub)
            else:
                if i == k:
                    sub_occur[i][k] = 0
                else:
                    sub_occur[i][k] = -15
            if msubstitutions[i][k] != 0:
                if i == k:
                    msub_occur[i][k] = 0
                else:
                    msub_occur[i][k] = float(msubstitutions[i][k])/float(sum_msub)
            else:
                if i == k:
                    msub_occur[i][k] = 0
                else:
                    msub_occur[i][k] = -15

            sub_on_msub_norm[i][k] = sub_occur[i][k]-msub_occur[i][k]
            if msubstitutions[i][k] != 0:
                sub_on_msub[i][k] = (float(substitutions[i][k]))/msubstitutions[i][k]
            else:
                sub_on_msub[i][k] = str(substitutions[i][k]) + '/0'

    sub_occur.to_csv('sub_occur.csv')
    msub_occur.to_csv('msub_occur.csv')
    sub_on_msub_norm.to_csv('sub_on_msub_norm.csv')
    sub_on_msub.to_csv('sub_on_msub.csv')