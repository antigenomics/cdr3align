import os
import pandas as pd
import ch_antigen_base as ch_w

class Found(Exception): pass

def readannotations(inppath, outpath):
    dict = {'HomoSapiens':{}}

    with open(inppath, 'r') as inp:
        for line in inp:
            for line in inp:
                inf = line.strip().split()
                if inf[1] in dict:
                    if inf[8] in dict[inf[1]]:
                        if inf[10][2] in dict[inf[1]][inf[8]]:
                            if inf[12] in dict[inf[1]][inf[8]][inf[10][2]]:
                                dict[inf[1]][inf[8]][inf[10][2]][inf[12]] += [[inf[0], inf[15], [inf[9],inf[7]], [inf[13],inf[14]]]]

                            else:
                                dict[inf[1]][inf[8]][inf[10][2]][inf[12]]=[[inf[0], inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]]
                        else:
                            dict[inf[1]][inf[8]][inf[10][2]]={inf[12]:[[inf[0], inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]]}
                    else:
                        dict[inf[1]][inf[8]]={inf[10][2]:{inf[12]:[[inf[0], inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]]}}
                else:
                    dict[inf[1]]={inf[8]:{inf[10][2]:{inf[12]:[[inf[0], inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]]}}}

    #anitgen {A{CDR1[[pdb, seq, chain, region]],CDR2[[pdb, seq, chain, region]],CDR3[[pdb, seq, chain, region]]},B{CDR1,CDR2,CDR3}}

    with open(outpath, 'w') as out:
        out.write('species\tantigen\tchain_allele\tCDR\tpdb\tCDR_seq\tpdb_cdr\tpdb_antigen\tcdr_start\tcdr_end\n')
        for org in dict:
            for antigen in dict[org]:
                print(len(dict[org][antigen]))#for chain in sorted(dict[org][antigen]):
                chain = 'A'
                for seq in range(len(dict[org][antigen][chain]['CDR1'])):
                    for ch in ['A','B']:
                        info=dict[org][antigen][ch]
                        for i in range(1,4):
                            cdr = 'CDR'+str(i)
                            if cdr in info:
                                out.write(org+'\t'+antigen+'\t'+ch+'\t'+cdr+'\t'+info[cdr][seq][0]+'\t'
                                          +info[cdr][seq][1]+'\t'+info[cdr][seq][2][0]+'\t'+info[cdr][seq][2][1]+'\t'
                                          +info[cdr][seq][3][0]+'\t'+info[cdr][seq][3][1]+'\n')

#0 pdb_id 2f53
#1 species HomoSapiens
#2 chain_mhc_a A
#3 mhc_a_allele blabla
#4 chain_mhc_b B
#5 mhc_b_allele blabla
#6 mhc_type MHCI
#7 chain_antigen C
#8 antigen_seq !!!
#9 chain_tcr D
#10 tcr_v_allele TRAV
#11 tcr_j_allele TRAJ
#12 tcr_region CDR1
#13 tcr_region_start 27
#14 tcr_region_end 33
#15 tcr_region_seq DSAIYN
#anitgen {pdb{A{CDR1[[seq, chain, region]],CDR2[[seq, chain, region]],CDR3[[seq, chain, region]]},B{CDR1,CDR2,CDR3}}}
readannotations('final.annotations.txt', 'sorted.annotations.txt')
def writecdrs(inppath, outpath):
    dict = {'HomoSapiens':{}}
    with open(inppath, 'r') as inp:
        for line in inp:
            for line in inp:
                inf = line.strip().split()
                if inf[1] in dict:
                    if inf[8] in dict[inf[1]]:
                        if inf[0] in dict[inf[1]][inf[8]]:
                            if inf[10][2] in dict[inf[1]][inf[8]][inf[0]]:
                                dict[inf[1]][inf[8]][inf[0]][inf[10][2]][inf[12]]=[inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]
                            else:
                                dict[inf[1]][inf[8]][inf[0]][inf[10][2]]={inf[12]:[inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]}
                        else:
                            dict[inf[1]][inf[8]][inf[0]]={inf[10][2]:{inf[12]:[inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]}}
                    else:
                        dict[inf[1]][inf[8]]={inf[0]:{inf[10][2]:{inf[12]:[inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]}}}
                else:
                    dict[inf[1]]={inf[8]:{inf[0]:{inf[10][2]:{inf[12]:[inf[15], [inf[9],inf[7]],[inf[13],inf[14]]]}}}}

    cdr1 = open(os.path.join(outpath,'cdr1.annotation.txt'), 'w')
    cdr1.write('species\tantigen\tchain_allele\tpdb\tCDR_seq\tpdb_cdr\tpdb_antigen\tcdr_start\tcdr_end\n')
    cdr2 = open(os.path.join(outpath,'cdr2.annotation.txt'), 'w')
    cdr2.write('species\tantigen\tchain_allele\tpdb\tCDR_seq\tpdb_cdr\tpdb_antigen\tcdr_start\tcdr_end\n')
    cdr3 = open(os.path.join(outpath,'cdr3.annotation.txt'), 'w')
    cdr3.write('species\tantigen\tchain_allele\tpdb\tCDR_seq\tpdb_cdr\tpdb_antigen\tcdr_start\tcdr_end\n')
    for org in dict:
        for antigen in dict[org]:
            if len(dict[org][antigen]) > 1:
                #{cdr3:{chain:{sequence,pdb}}, cdr2:{sequence,pdb}, cdr1:{sequence,pdb}}
                dict2 = {}
                pdbs = {}
                for pdb in dict[org][antigen]:
                    for chain in dict[org][antigen][pdb]:
                        for cdr in dict[org][antigen][pdb][chain]:
                            if cdr in dict2:
                                if chain in dict2[cdr]:
                                    if dict[org][antigen][pdb][chain][cdr][0] in dict2[cdr][chain]:
                                        pass
                                    else:
                                        dict2[cdr][chain][dict[org][antigen][pdb][chain][cdr][0]]=pdb
                                else:
                                    dict2[cdr][chain] = {}
                            else:
                                dict2[cdr] = {chain:{}}
                print(dict2)
                for cdr in dict2:
                    for chain in sorted(dict2[cdr]):
                        if len(dict2[cdr][chain]) > 1:
                            for seq in dict2[cdr][chain]:
                                pdb = dict2[cdr][chain][seq]
                                chcdr = dict[org][antigen][pdb][chain][cdr][1][0]
                                chantigen = dict[org][antigen][pdb][chain][cdr][1][1]
                                cdrs = dict[org][antigen][pdb][chain][cdr][2][0]
                                cdre = dict[org][antigen][pdb][chain][cdr][2][1]
                                if cdr=='CDR1':
                                    cdr1.write(org + '\t' + antigen + '\t' + chain + '\t' + pdb + '\t'
                                        + seq + '\t' + chcdr + '\t' + chantigen + '\t'
                                        + cdrs + '\t' + cdre + '\n')
                                elif cdr=='CDR2':
                                    cdr2.write(org + '\t' + antigen + '\t' + chain + '\t' + pdb + '\t'
                                               + seq + '\t' + chcdr + '\t' + chantigen + '\t'
                                               + cdrs + '\t' + cdre + '\n')
                                elif cdr=='CDR3':
                                    cdr3.write(org + '\t' + antigen + '\t' + chain + '\t' + pdb + '\t'
                                               + seq + '\t' + chcdr + '\t' + chantigen + '\t'
                                               + cdrs + '\t' + cdre + '\n')

    cdr1.close()
    cdr2.close()
    cdr3.close()
writecdrs('final.annotations.txt', 'seqdata/HomoSapiens/cdrs/')
def writecdrs_2(inppath, outpath):
    all = []
    cdr = {}
    cdr['cdr1'] = set()
    cdr['cdr2'] = set()
    cdr['cdr3'] = set()
    cdr1_2 = []
    cdr2_3 = []
    cdr3_1 = []
    for i in range(1,4):
        with open(os.path.join(inppath, 'cdr'+str(i)+'.annotation.txt'), 'r') as inp:
            cd = 'cdr'+str(i)
            for line in inp:
                for line in inp:
                    cdr[cd].add(line.strip().split()[1]+'\t'+line.strip().split()[3])

    print(len(cdr['cdr1']))
    for i in cdr['cdr1']:
        if i in cdr['cdr2'] and i in cdr['cdr3']:
            all.append(i)
        elif i in cdr['cdr2']:
            cdr1_2.append(i)
        elif i in cdr['cdr3']:
            cdr3_1.append(i)
    for i in cdr['cdr2']:
        if i in cdr['cdr3'] and i not in cdr['cdr1']:
            cdr2_3.append(i)

    for i in all:
        cdr['cdr1'].remove(i)
        cdr['cdr2'].remove(i)
        cdr['cdr3'].remove(i)
    for i in cdr1_2:
        cdr['cdr1'].remove(i)
        cdr['cdr2'].remove(i)
    for i in cdr2_3:
        cdr['cdr2'].remove(i)
        cdr['cdr3'].remove(i)
    for i in cdr3_1:
        cdr['cdr3'].remove(i)
        cdr['cdr1'].remove(i)

    with open(outpath, 'w') as out:
        out.write('all\n'+'\n'.join(sorted(all))+'\n//\n')
        out.write('cdr1_2\n' + '\n'.join(sorted(cdr1_2)) + '\n//\n')
        out.write('cdr2_3\n' + '\n'.join(sorted(cdr2_3)) + '\n//\n')
        out.write('cdr3_2\n' + '\n'.join(sorted(cdr3_1)) + '\n//\n')
        out.write('cdr1\n' + '\n'.join(sorted(cdr['cdr1'])) + '\n//\n')
        out.write('cdr2\n' + '\n'.join(sorted(cdr['cdr2'])) + '\n//\n')
        out.write('cdr3\n' + '\n'.join(sorted(cdr['cdr3'])) + '\n//\n')

writecdrs_2('seqdata/HomoSapiens/cdrs/', 'seqdata/HomoSapiens/cdrs/cdrs.annotation.txt')
