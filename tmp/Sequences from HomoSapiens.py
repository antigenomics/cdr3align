import os
import glob
import pandas as pd
import csv
import json
import time
import copy
import math
import substitution_matrix as ch2
import graphs_var1 as ch1


def crdir(path):
    if os.path.exists(path):
        return 0
    else:
        os.mkdir(path)

orgname = 'HomoSapiens'

def getdata(orgname):
    seqspecies = {}
    with open('vdjdb-2016-09-05/vdjdb.slim.txt', 'r') as inp:
        for line in inp:
            inf = line.strip().split()
            if inf[0] == 'TRB' and inf[2] == orgname:
                if inf[3] in seqspecies:
                    seqspecies[inf[3]].append(inf[1])
                else:
                    seqspecies[inf[3]]=[inf[1]]

    crdir('seqdata')
    with open('seqdata/trb_' + orgname + 'AntigSeq.txt', 'w') as out:
        for i in seqspecies:
            out.write('antigen.epitope:\t'+str(i)+'\t'+str(len(seqspecies[i]))+'\n')
            out.write('cdr3:\t'+'\t'.join(seqspecies[i]))
            out.write('\n//\n')

    teachseq = []
    with open('trb_sequences_all_teach', 'r') as inp:
        for line in inp:
            teachseq.append(line.strip().split()[0])

    firsttest = []
    with open('trb_sequences_half_first_testing', 'r') as inp:
        for line in inp:
            firsttest.append(line.strip().split()[0])

    sacredtest = []
    with open('trb_sequences_half_sacred_testing', 'r') as inp:
        for line in inp:
            sacredtest.append(line.strip().split()[0])

    crdir('seqdata/'+orgname)
    humanteachseq = {}
    humantestseq = {}
    for i in seqspecies: #i – antigen
        for k in seqspecies[i]: #k - sequence
            if k in teachseq:
                if i in humanteachseq:
                    humanteachseq[i].append(k)
                else:
                    humanteachseq[i] = [k]
            else:
                if i in humantestseq:
                    humantestseq[i].append(k)
                else:
                    humantestseq[i] = [k]

    with open('seqdata/'+orgname+'/trb_sequences_'+orgname+'_teach', 'w') as out:
        for i in humanteachseq:
            out.write('\n'.join(humanteachseq[i])+'\n')

    with open('seqdata/' + orgname + '/trb_sequences_' + orgname + '_teach_antigen', 'w') as out:
        for i in humanteachseq:
            out.write('\n'.join(list(map(lambda x: x+'\t'+i, humanteachseq[i])))+'\n')

    with open('seqdata/'+orgname+'/trb_sequences_'+orgname+'_test', 'w') as out:
        for i in humantestseq:
            out.write('\n'.join(humantestseq[i])+'\n')


    humanftestseq = {}
    humanstestseq = {}
    for i in humantestseq:
        for k in humantestseq[i]:
            if k in firsttest:
                if i in humanftestseq:
                    humanftestseq[i].append(k)
                else:
                    humanftestseq[i] = [k]
            elif k in sacredtest:
                if i in humanstestseq:
                    humanstestseq[i].append(k)
                else:
                    humanstestseq[i] = [k]
            else:
                print(i, k, 'something strange!')

    with open('seqdata/'+orgname+'/trb_sequences_'+orgname+'_half_first_testing', 'w') as out:
        for i in humanftestseq:
            out.write('\n'.join(humanftestseq[i])+'\n')

    with open('seqdata/' + orgname + '/trb_sequences_' + orgname + '_half_first_testing_antigen', 'w') as out:
        for i in humanftestseq:
            out.write('\n'.join(list(map(lambda x: x + '\t' + i, humanftestseq[i]))) + '\n')

    with open('seqdata/'+orgname+'/trb_sequences_'+orgname+'_half_sacred_testing', 'w') as out:
        for i in humanstestseq:
            out.write('\n'.join(humanstestseq[i])+'\n')

    with open('seqdata/' + orgname + '/trb_sequences_' + orgname + '_half_sacred_testing_antigen', 'w') as out:
        for i in humanstestseq:
            out.write('\n'.join(list(map(lambda x: x + '\t' + i, humanstestseq[i]))) + '\n')

#getdata('orgname')
def gethseq(inp):
    hseq = []
    with open(inp, 'r') as inp:
        inf = inp.read()
        hseq += inf.split('\n')
        #print(hseq)
        return hseq

def get_org_pairalignment_inrange(inp, r=4):
    for s in range(r):
        for i in range(r):
            for d in range(r):
                tot = s+i+d
                print(s,i,d,tot)
                with open('sub_del_ins_tot/out_teach_'+str(s)+'_'+str(i)+'_'+str(d)+'_'+str(tot), 'r') as inp:
                    with open('seqdata/'+orgname+'/sub_del_ins_tot/out_teach_'+str(s)+'_'+str(i)+'_'+str(d)+'_'+str(tot), 'w') as out:
                        for line in inp:
                            seqs = line.replace('-','').upper().strip().split()
                            if seqs[0] in hseq and seqs[1] in hseq:
                                out.write(line)


def groovy_treesearch(inp, opt, out):
    os.system("groovy TreeSearch.groovy " + inp + " " + opt + " " + out)

#groovy_treesearch('seqdata/'+orgname+'/trb_sequences_'+orgname+'_teach', '7,0,0,7', 'seqdata/'+orgname+'/sub_del_ins_tot/out_teach_7_0_0_7')
def get_alignments_set(wrkfolder, alignm):
    crdir(wrkfolder+'pairalignments')
    for i in alignm:
        algn = set()
        with open(i, 'r') as inp:
            for line in inp:
                algn.add(line.strip().replace('-','').upper())
        with open(wrkfolder+'pairalignments/'+i.split('/')[3]+'_set', 'w') as out:
            out.write('\n'.join(algn))


def get_information_from_alignment(alignmentset, sequences, antigens):
    sequences_2 = copy.deepcopy(sequences)
    antigens_2 = copy.deepcopy(antigens)
    with open(alignmentset, 'r') as inp:
        aaa = 0
        bbb = 0
        for line in inp:
            seq1, seq2 = line.strip().split()[0], line.strip().split()[1]
            sequences_2[seq1].add_link(sequences_2[seq2])
            sequences_2[seq2].add_link(sequences_2[seq1])
            aaa = seq1
            bbb = seq2
        for i in sequences_with_antigen:
            for k in sequences_with_antigen[i]:
                antigens_2[k].add_seq(sequences_2[i])
    return sequences_2, antigens_2

def write_information_from_matrices(alignmentset, outputfile, sequences, antigens):
    sequences_2, antigens_2 = get_information_from_alignment(alignmentset, sequences, antigens)

    listantigens = []
    for k in antigens_2:
        listantigens.append([k, antigens_2[k]])
    listantigens.sort(key = lambda x: len(x[1].sequences))
    with open(outputfile, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['antigen', 'number_of_sequences', 'possible_inlinks', 'real_inlinks', 'possible_outlinks',
                        'real_outlinks', 'number_of_boseq', 'inlinks_of_boseq', 'outlinks_of_boseq',
                        'possibe_boseq_links_to_other_antigens', 'real_boseq_links_to_other_antigens'])
        for k in listantigens:
            wantigen = k[0]
            wseqnumber = len(antigens_2[wantigen].sequences)+len(antigens_2[wantigen].boundarysequences)
            wposinlinks = ((wseqnumber**2)-wseqnumber)
            wrealinlinks = (antigens_2[wantigen].totinlinks)
            wposoutlinks = wseqnumber*(len(sequences_with_antigen)-wseqnumber)
            wrealoutlinks = antigens_2[wantigen].totoutlinks
            boseq = antigens_2[wantigen].boundarysequences
            if len(boseq) > 0:
                wboseqnumber = len(antigens_2[wantigen].boundarysequences)
                boseqin = 0
                boseqout = 0
                posboseqother = 0
                realboseqother = 0
                for i in boseq:
                    realboseqother += i.outlinked(wantigen)
                    boseqout += len(i.outlinks)+realboseqother
                    for antig in i.antigen:
                        if antig == wantigen:
                            pass
                        else:
                            #print(antig)
                            posboseqother += len(antigens_2[antig].sequences)
                            for seq in antigens_2[antig].boundarysequences:  #...
                                if wantigen in seq.antigen:
                                    #print(wantigen, seq.antigen, 'pass')
                                    pass
                                else:
                                    #print(wantigen, seq.antigen, 'not pass')
                                    posboseqother += 1
                    boseqin += len(i.inlinks[wantigen])
                wboseqinlinks = boseqin
                wboseqoutlinks = boseqout
                wposboseqother = posboseqother
                wrealboseqother = realboseqother
            else:
                wboseqnumber, wboseqinlinks, wboseqoutlinks, wposboseqother, wrealboseqother = 0, 0, 0, 0, 0
            writer.writerow([wantigen, wseqnumber, wposinlinks, wrealinlinks, wposoutlinks, wrealoutlinks, wboseqnumber, wboseqinlinks, wboseqoutlinks, wposboseqother, wrealboseqother])

def get_info_from_matrices(folder, outfolder, glword):
    alignms = glob.glob(folder+'/'+glword+'*')
    print(alignms)
    info = {}
    for name in alignms:
        with open(name, 'r') as inp:
            reader = csv.reader(inp, delimiter='\t')
            info[name] = {'realin':0,'realout':0, 'realboseqout':0, 'realboseqother':0, 'posin':0, 'posout':0, 'realboseqin':0, 'posboseqother':0}
            for row in reader:
                for row in reader:
                    info[name]['realin'] += int(row[3])
                    info[name]['realout'] += int(row[5])
                    info[name]['realboseqout'] += int(row[8])
                    info[name]['realboseqother'] += int(row[10])
                    info[name]['posin'] += int(row[2])
                    info[name]['posout'] += int(row[4])
                    info[name]['realboseqin'] += int(row[7])
                    info[name]['posboseqother'] += int(row[9])
    sortedinfo = []
    for name in info:
        sortedinfo.append([name, info[name]])
    sortedinfo.sort(key = lambda x: x[0])
    with open(outfolder+'sum_of_info_from_alignments', 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['name', 'rpin', 'rpout', 'rpin-rpout', 'realin', 'realout', 'realboseqin', 'realboseqout',
                         'realboseqother', 'posin', 'posout', 'posboseqother'])
        for inf in sortedinfo:
            name = inf[0]
            print(name)
            writer.writerow([name.split('/')[-1], "{0:.5f}".format(float(info[name]['realin'])/(info[name]['posin']+info[name]['posout'])),
                             "{0:.5f}".format(float(info[name]['realout'])/(info[name]['posin']+info[name]['posout'])),
                             "{0:.5f}".format((float(info[name]['realin']) / (info[name]['posin']+info[name]['posout']))-(float(info[name]['realout']) / (info[name]['posin']+info[name]['posout']))),
                             info[name]['realin'], info[name]['realout'],info[name]['realboseqin'], info[name]['realboseqout'],
                             info[name]['realboseqother'], info[name]['posin'], info[name]['posout'], info[name]['posboseqother']])

def get_info_from_matrices_2(folder, outfolder, glword):
    alignms = glob.glob(folder + '/' + glword + '*')
    print(alignms)
    info = {}
    for name in alignms:
        with open(name, 'r') as inp:
            reader = csv.reader(inp, delimiter='\t')
            info[name] = {'realin': 0, 'realout': 0, 'realboseqout': 0, 'realboseqother': 0, 'posin': 0,
                          'posout': 0, 'realboseqin': 0, 'posboseqother': 0}
            for row in reader:
                for row in reader:
                    info[name]['realin'] += int(row[3])
                    info[name]['realout'] += int(row[5])
                    info[name]['realboseqout'] += int(row[8])
                    info[name]['realboseqother'] += int(row[10])
                    info[name]['posin'] += int(row[2])
                    info[name]['posout'] += int(row[4])
                    info[name]['realboseqin'] += int(row[7])
                    info[name]['posboseqother'] += int(row[9])
    sortedinfo = []
    for name in info:
        sortedinfo.append([name, info[name]])
    sortedinfo.sort(key=lambda x: x[0])
    with open(outfolder + 'sum_of_info_from_alignments_2', 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(
            ['name', 'rpin', 'rpout', 'rpin-rpout', 'realin', 'realout', 'realboseqin', 'realboseqout',
             'realboseqother', 'posin', 'posout', 'posboseqother'])
        for inf in sortedinfo:
            name = inf[0]
            print(name)
            writer.writerow([name.split('/')[-1],
                             "{0:.5f}".format(float(info[name]['realin']) / (
                             info[name]['posin'] + info[name]['posout'] - info[name]['realout'] -
                             info[name]['realin'])),
                             "{0:.5f}".format(float(info[name]['realout']) / (
                             info[name]['posin'] + info[name]['posout'] - info[name]['realout'] -
                             info[name]['realin'])),
                             "{0:.5f}".format((float(info[name]['realin']) / (
                             info[name]['posin'] + info[name]['posout'])) - (
                                              float(info[name]['realout']) / (
                                              info[name]['posin'] + info[name]['posout'] - info[name][
                                                  'realout'] - info[name]['realin']))),
                             info[name]['realin'],
                             info[name]['realout'],
                             info[name]['realboseqin'],
                             info[name]['realboseqout'],
                             info[name]['realboseqother'],
                             info[name]['posin'],
                             info[name]['posout'],
                             info[name]['posboseqother']])


'''
hseq = gethseq('seqdata/'+orgname+'/trb_sequences_'+orgname+'_teach')
crdir('seqdata/'+orgname+'/sub_del_ins_tot')

sequences_with_antigen = {}
with open('seqdata/'+orgname+'/trb_sequences_'+orgname+'_teach_antigen', 'r') as inp:
    for line in inp:
        inf = line.strip().split()
        if inf[0] in sequences_with_antigen:
            sequences_with_antigen[inf[0]].append(inf[1])
        else:
            sequences_with_antigen[inf[0]] = [inf[1]]

wrkfolder = 'seqdata/'+orgname+'/'

alignm = glob.glob(wrkfolder+'sub_del_ins_tot/out_teach*')
sequences, antigens = ch1.get_seqs_antigens(sequences_with_antigen)
#alignmentsets = glob.glob(wrkfolder+'pairalignments/out_teach*')
#crdir(wrkfolder+'/information_from_alignments')
#for i in alignmentsets:
#    write_information_from_matrices(i, wrkfolder+'/information_from_alignments/alignments_'+i.split('out_teach')[1], sequences, antigens)
'''
'''
wrkfolder = 'seqdata/'+orgname+'/'
get_info_from_matrices_2(wrkfolder+'information_from_alignments', wrkfolder, 'alignments')
'''
'''
alignments = {}
inf = '4_0_0_4'
with open('seqdata/HomoSapiens/sub_del_ins_tot/out_teach_'+inf, 'r') as inp:
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

working_occurence = copy.deepcopy(occurence)
substitutions = pd.DataFrame(0, index=occurence, columns=occurence)
msubstitutions = pd.DataFrame(0, index=occurence, columns=occurence)

crdir(wrkfolder+'pictures')

#(occurence, alignments, substitutions, msubstitutions, antig_seq, name, work_occur=None):
#substitutions, msubstitutions, working_occurence, mworking_occurence = ch2.get_occurence_matrix_2(occurence, alignments, substitutions, msubstitutions, antig_seq=sequences_with_antigen, name=inf, folder=wrkfolder, work_occur=working_occurence)
substitutions, msubstitutions, working_occurence, mworking_occurence = ch2.get_occurence_matrix_diagonal(occurence, alignments, substitutions, msubstitutions, antig_seq=sequences_with_antigen, name=inf+'_diagonal', folder=wrkfolder, work_occur=working_occurence)
'''