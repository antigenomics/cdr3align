import os
import ch_search_near_inpdb as ch_s
import ch_get_superimpose_structures as ch_get
import ch_antigen_base as ch_w


#cdrs = ['CDR1', 'CDR2', 'CDR3']
#inpfiles = ['seqdata/HomoSapiens/cdrs/cdr1.annotation_unique.txt', 'seqdata/HomoSapiens/cdrs/cdr2.annotation_unique.txt',
#            'seqdata/HomoSapiens/cdrs/cdr3.annotation_unique.txt']
inppdb = 'seqdata/HomoSapiens/pdbs/'
chainname = ['A', 'B']
cdrs = ['CDR3']
inpfiles = ['seqdata/HomoSapiens/cdrs/cdr3.annotation_unique.txt']
for i in range(len(cdrs)):
    print(i)
    pdbs = ch_get.get_pdbs(inpfiles[i], inppdb)
    pdbs_imposed = ch_s.get_pdbs_imposed(pdbs)
    for k in chainname:
        print(k)
        ch_s.getdistances(inppath=inppdb,
                          outpath='seqdata/HomoSapiens/distances_from_pdb_'+str(k)+str(cdrs[i])+'.txt',
                          #outpath = 'seqdata/HomoSapiens/distances_test.txt',
                          pdbs=pdbs,
                          pdbs_imposed=pdbs_imposed,
                          chainname=k,
                          cdrname=cdrs[i])
