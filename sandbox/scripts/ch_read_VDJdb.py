import sys, os, argparse
import ch_antigen_base as ch_w

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description =
"This script will get info from VDJdb and write it in another file. \
Current version is %s" % curr_version
)
parser.add_argument("-i", help="Path to VDJdb", required=True, dest="input_dir")
parser.add_argument("-w", help="Path to working directory (use . to work where you run script)", required=True,
                    dest="work_dir")
parser.add_argument("-o", help="Prefix for the output files. Default is ''.", required=False,
                    dest='output')
parser.add_argument("-gene", help="Name of the gene from database", required=False,
                    dest="gene")
parser.add_argument("-org", help="Name of the organism from which you will get sequences", required=False,
                    dest="orgname")
parser.add_argument("-q", help="Would you like to be asked to do more, than this script should do \
                    (like running next scripts) (y/n)?", required=False, dest="q")

myargs = parser.parse_args()

if myargs.output is None:
    myargs.output = ""
else:
    myargs.output = myargs.output+"_"

ch_w_dict = ch_w.getdata(myargs.input_dir, myargs.gene, myargs.orgname)

output_name = myargs.output

if myargs.gene is not None:
    output_name += myargs.gene + "_"
else:
    output_name += "all_"
if myargs.orgname is not None:
    output_name += myargs.orgname
    ch_w.crdir(os.path.join(myargs.work_dir, myargs.orgname))
    output_name = os.path.join(myargs.orgname, output_name)
else:
    output_name += "all"
    ch_w.crdir(os.path.join(myargs.work_dir, 'all'))
    output_name = os.path.join('all', output_name)

with open(os.path.join(myargs.work_dir, output_name+"_seq_antigen.txt"), 'w') as out:
    for i in ch_w_dict["seqspecies"]:
        out.write('\n'.join(list(map(lambda x: str(x)+'\t'+str(i), ch_w_dict["seqspecies"][i])))+'\n')

with open(os.path.join(myargs.work_dir, output_name+"_seq"), 'w') as out:
    for i in ch_w_dict["seqspecies"]:
        out.write('\n'.join(ch_w_dict["seqspecies"][i]) + '\n')

if myargs.q == 'y':
    answer = input("Would you like to divide sequences into teaching and testing samples? (y/n)")
    if answer == 'y':
        answer = int(input("What is the minimum quantity of sequences interacting with antigen (lower quantities won't go to testing sample)?"))
        ch_w.test_and_teach(ch_w_dict["seqspecies"], answer, 2, os.path.join(myargs.work_dir, output_name))
