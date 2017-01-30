import os
import argparse
import ch_search_near_inpdb as ch_s
import ch_get_superimpose_structures as ch_get
import ch_antigen_base as ch_w

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script can superimpose structures and write down distances between "
                                             "superimposed cdr sequences. \
Current version is %s" % curr_version)

parser.add_argument("-i", help="Path to cdr annotation file. Default is '../seqdata/HomoSapiens/cdr/cdr3.annotation_unique.txt'",
                               required=False, dest="input_dir")
parser.add_argument("-superimpose", help="Do you need to superimpose files? Default is None",
                    required=False, dest="superimp")
parser.add_argument("-overwrite", help="if 'Y', then existing files can be overwritten during superimposing process. "
                                       "Default is False",
                    required=False, dest="overwrite")
parser.add_argument("-pdbs", help="Path tp pdb files. Default is '../seqdata/HomoSapiens/pdbs/'",
                    required=False, dest="pdbs")
parser.add_argument("-chain", help="Chain name. Default is 'B'", required=False, dest="chain")
parser.add_argument("-o", help="Path to output directory. '../seqdata/HomoSapiens/distances_from_pdb_test_Bcdr_3.txt'",
                    required=False, dest='output_dir')


myargs = parser.parse_args()

if myargs.input_dir is None:
    myargs.input_dir = '../seqdata/HomoSapiens/cdrs/cdr3.annotation_unique.txt'
if myargs.pdbs is None:
    myargs.pdbs = '../seqdata/HomoSapiens/pdbs/'
if myargs.chain is None:
    myargs.chain = 'B'
if myargs.output_dir is None:
    myargs.output_dir = '../seqdata/HomoSapiens/distances_from_pdb_test_Bcdr_3.txt'

ch_w.crdir(myargs.pdbs)

pdbs = ch_get.get_pdbs(myargs.input_dir, myargs.pdbs)
if myargs.superimp is not None and myargs.superimp != False:
    print('Start superimposing structures\n')
    if myargs.overwrite == 'Y':
        myargs.overwrite = True
    else:
        myargs.overwrite = False
    ch_get.superimpose(pdbs, myargs.pdbs, overwrite=myargs.overwrite)
    print('End superimposing structures\n')
pdbs_imposed = ch_s.get_pdbs_imposed(pdbs)
ch_s.getdistances(inppath=myargs.pdbs,
                  outpath=myargs.output_dir,
                  pdbs=pdbs,
                  pdbs_imposed=pdbs_imposed,
                  chainname=myargs.chain)
