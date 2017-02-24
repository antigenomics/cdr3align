import sys, os, argparse, glob
import ch_antigen_base as ch_w
import pandas as pd

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will make occurence matrices for pairwise alignments \
                                              from same antigen and from different antigens. \
Current version is %s" % curr_version)

parser.add_argument("-i", help="Path to pair alignments", nargs='+', required=True, dest="input_dir")
parser.add_argument("-iseq", help="Path to sequences with antigens used to create pairwise alignments",
                    required=True, dest="sequences")
parser.add_argument("-w", help="Path to working directory (use . to work where you run script)", required=True,
                    dest="work_dir")
parser.add_argument("-o", help="Name of output directory. Default is ''. ", required=False, dest='output_dir')
parser.add_argument("-use_all", help="if y, then it will recognize -i not as file, but as directory",
                    required=False, dest="use_all")
parser.add_argument("-q",
                    help="Would you like to be asked to do more, than this script should do \
                    (like running next scripts) (y/n)?", required=False, dest="q")

myargs = parser.parse_args()

occurence = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
             'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

if myargs.output_dir is None:
    myargs.output_dir = myargs.work_dir
else:
    myargs.output_dir = os.path.join(myargs.work_dir, myargs.output_dir)
    ch_w.crdir(myargs.output_dir)

sequences_with_antigen ={}
with open(myargs.sequences, 'r') as inp:
    for line in inp:
        inf = line.strip().split()
        if inf[0] in sequences_with_antigen:
            sequences_with_antigen[inf[0]].append(inf[1])
        else:
            sequences_with_antigen[inf[0]] = [inf[1]]

if myargs.use_all == 'y':
    files = glob.glob(os.path.join(myargs.input_dir[0], '') + '*')
    for i in files:
        alignments = ch_w.get_alignments(i)
        substitutions = pd.DataFrame(0, index=sorted(occurence), columns=sorted(occurence))
        msubstitutions = pd.DataFrame(0, index=sorted(occurence), columns=sorted(occurence))
        ch_w.get_occurence_matrix(occurence, alignments, substitutions, msubstitutions, antig_seq=sequences_with_antigen,
                               name=i.split('/')[-1], folder=myargs.output_dir)
else:
    for i in myargs.input_dir:
        alignments = ch_w.get_alignments(i)
        substitutions = pd.DataFrame(0, index=sorted(occurence), columns=sorted(occurence))
        msubstitutions = pd.DataFrame(0, index=sorted(occurence), columns=sorted(occurence))
        ch_w.get_occurence_matrix(occurence, alignments, substitutions, msubstitutions, antig_seq=sequences_with_antigen,
                                  name=i.split('/')[-1], folder=myargs.output_dir)

if myargs.q == 'y':
    print("We are sorry, right now this script can't do anything more.")