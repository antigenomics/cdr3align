import sys, os, argparse
import ch_antigen_base as ch_w

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description =
"This script will divide your sequences into teaching and testing samples. \
Current version is %s" % curr_version
)
parser.add_argument("-i", help="Path to your sequences", required=True,
                    dest="input_dir")
parser.add_argument("-w", help="Path to working directory (use . to work where you run script)", required=True,
                    dest="work_dir")
parser.add_argument("-o", help="Prefix for the output files. Default is ''.", required=False, dest='output')
parser.add_argument("-q", help="Would you like to be asked to do more, than this script should do\
                    (like running next scripts) (y/n)?", required=False, dest="q")

myargs = parser.parse_args()

if myargs.output is None:
    myargs.output = ""
else:
    myargs.output = myargs.output+"_"

seqspecies = {}
with open(myargs.input_dir, 'r') as inp:
    for line in inp:
        seqant = line.strip().split()
        if seqant[1] in seqspecies:
            seqspecies[seqant[1]].append(seqant[0])
        else:
            seqspecies[seqant[1]] = [seqant[0]]
output_name = myargs.input_dir.split('/')[-1].replace('_antigen.txt', '').replace('.txt', '')

answer = int(input(
    "What is the minimum quantity of sequences interacting with antigen (lower quantities won't go to testing sample)?"))

ch_w.test_and_teach(seqspecies, answer, 2, os.path.join(myargs.work_dir, output_name))

if myargs.q == 'y':
    print("We are sorry, right now this script can't do anything more.")