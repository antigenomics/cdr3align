import sys, os, argparse, glob
import ch_antigen_base as ch_w

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will get information from pairwise alignments. \
Current version is %s" % curr_version)

parser.add_argument("-i", help="Path to alignment sets", nargs='+', required=True, dest="input_dir")
parser.add_argument("-iseq", help="Path to sequences with antigens used to create pairwise alignments",
                    required=False, dest="sequences")
parser.add_argument("-w", help="Path to working directory (use . to work where you run script)", required=True,
                    dest="work_dir")
parser.add_argument("-o", help="Name of output directory. Default is 'information_from_alignments'.", required=False, dest='output_dir')
parser.add_argument("-use_all", help="if y, then it will recognize -i not as file, but as directory",
                    required=False, dest="use_all")
parser.add_argument("-q",
                    help="Would you like to be asked to do more, than this script should do \
                    (like running next scripts) (y/n)?", required=False, dest="q")
parser.add_argument("-sum_only",
                    help="If y, then you'll get only sum of information (separate files required in -i (-use_all still \
                    working))", required=False, dest="sum_only")

myargs = parser.parse_args()

if myargs.sum_only != 'y':
    if myargs.output_dir is None:
        myargs.output_dir = 'information_from_alignments'
        ch_w.crdir(os.path.join(myargs.work_dir, myargs.output_dir))

    if myargs.sequences is None:
        raise Exception('No path to sequences with antigens. Use -iseq to add path.')
    sequences_with_antigen = {}
    with open(myargs.sequences, 'r') as inp:
        for line in inp:
            inf = line.strip().split()
            if inf[0] in sequences_with_antigen:
                sequences_with_antigen[inf[0]].append(inf[1])
            else:
                sequences_with_antigen[inf[0]] = [inf[1]]

    sequences, antigens = ch_w.get_seqs_antigens(sequences_with_antigen)

    if myargs.use_all == 'y':
        files = glob.glob(os.path.join(myargs.input_dir[0], '') + '*')
        for i in files:
            ch_w.write_information_from_matrices(i, os.path.join(myargs.work_dir, myargs.output_dir,
                                                                 'alignments_'+str(i.split('/')[-1])),
                                                 sequences, antigens, sequences_with_antigen)
    else:
        for i in myargs.input_dir:
            ch_w.write_information_from_matrices(i, os.path.join(myargs.work_dir, myargs.output_dir,
                                                                 'alignments_' + str(i.split('/')[-1])),
                                                 sequences, antigens, sequences_with_antigen)

    if myargs.q == 'y':
        answer = input("Do you want to get sum of information from alignments (will be placed in workdir)?")
        if answer == 'y':
            if myargs.use_all == 'y':
                files = list(map(lambda x: os.path.join(myargs.work_dir, myargs.output_dir,
                                                        'alignments_' + str(x.split('/')[-1])), files))
                ch_w.get_info_from_matrices(myargs.work_dir, files)
            else:
                files = list(map(lambda x: os.path.join(myargs.work_dir, myargs.output_dir,
                                                        'alignments_' + str(x.split('/')[-1])), myargs.input_dir))
                ch_w.get_info_from_matrices(myargs.work_dir, files)
else:
    if myargs.output_dir is None:
        myargs.output_dir = ''

    if myargs.use_all == 'y':
        files = glob.glob(os.path.join(myargs.input_dir[0], '') + '*')
        ch_w.get_info_from_matrices(os.path.join(myargs.work_dir, myargs.output_dir), files)
    else:
        ch_w.get_info_from_matrices(os.path.join(myargs.work_dir, myargs.output_dir), myargs.input_dir)

    if myargs.q == 'y':
        print("We are sorry, right now this script can't do anything more.")