import sys, os, argparse
import ch_antigen_base as ch_w

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description =
"This script will use groovy_treesearch to find all possible pairwise alignments for current range. \
Current version is %s" % curr_version
)
parser.add_argument("-i", help="Path to sequences", required=True, dest="input_dir")
parser.add_argument("-w", help="Path to working directory (use . to work where you run script)", required=True,
                    dest="work_dir")
parser.add_argument("-o", help="Name of output directory. Default is 'sub_del_ins_tot'.", required=False,
                    dest='output_dir')
parser.add_argument("-range", help="Range for deletions, subtitutions and insertions", required=True,
                    dest="range")
parser.add_argument("-q",
                    help="Would you like to be asked to do more, than this script should do \
                    (like running next scripts) (y/n)?", required=False, dest="q")

myargs = parser.parse_args()

if myargs.output_dir is None:
    myargs.output_dir = 'sub_del_ins_tot'

ch_w.crdir(os.path.join(myargs.work_dir, myargs.output_dir))

for s in range(int(myargs.range)):
    for i in range(int(myargs.range)):
        for d in range(int(myargs.range)):
            tot = s+i+d
            opt = ','.join(map(str, [s,i,d,tot]))
            output_name = 'out_'+'_'.join(map(str, [s,i,d,tot]))
            ch_w.groovy_treesearch(myargs.input_dir, opt, os.path.join(myargs.work_dir, myargs.output_dir, output_name))

if myargs.q == 'y':
    print("We are sorry, right now this script can't do anything more.")