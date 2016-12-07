import sys, os, argparse, glob
import ch_antigen_base as ch_w

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will use pairwise alignments and get unique set of them. \
Current version is %s" % curr_version)

parser.add_argument("-i", help="Path to alignments", nargs='+', required=True, dest="input_dir")
parser.add_argument("-w", help="Path to working directory (use . to work where you run script)", required=True,
                    dest="work_dir")
parser.add_argument("-o", help="Name of output directory. Default is 'pairalignments'.", required=False, dest='output_dir')
parser.add_argument("-use_all", help="if y, then it will recognize -i not as file, but as directory",
                    required=False, dest="use_all")
parser.add_argument("-q",
                    help="Would you like to be asked to do more, than this script should do \
                    (like running next scripts) (y/n)?", required=False, dest="q")

myargs = parser.parse_args()

if myargs.output_dir is None:
    myargs.output_dir = 'pairalignments'

ch_w.crdir(os.path.join(myargs.work_dir, myargs.output_dir))

if myargs.use_all == 'y':
    files = glob.glob(os.path.join(myargs.input_dir[0],'')+'*')
    ch_w.get_alignments_set(myargs.work_dir, files, myargs.output_dir)
else:
    algn = set()
    for i in myargs.input_dir:
        with open(i, 'r') as inp:
            for line in inp:
                algn.add(line.strip().replace('-', '').upper())
        with open(os.path.join(myargs.work_dir, myargs.output_dir, i.split('/')[-1] + '_set'), 'w') as out:
            out.write('\n'.join(algn))

if myargs.q == 'y':
    print("We are sorry, right now this script can't do anything more.")