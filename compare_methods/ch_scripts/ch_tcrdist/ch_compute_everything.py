from ch_base import *
import argparse

#==============================

curr_version = 1.0
parser = argparse.ArgumentParser(description="Pipeline, that will compute everything for tcr_dist")

parser.add_argument("-vdjdb_path", type=str, default="../../imgt_work")
parser.add_argument("-distance_type", type=str, default="default")
parser.add_argument("-epitope_col", type=str, default="epitope.seq")
parser.add_argument("-multiproc", type=int, default=40)
parser.add_argument("-nrandom", type=int, default=1000)
parser.add_argument("-nlim", type=int, default=30)
parser.add_argument("-compare_methods_path", type=str, default="../ch_compare_methods")
parser.add_argument("-datasets", type=str, nargs="+", default="A.S;M.D;P.T;filtered".split(";"),
                    help="datasets.")
parser.add_argument("-organisms", type=str, nargs="+", default="HomoSapiens;MusMusculus".split(";"))
parser.add_argument("-chains", type=str, nargs="+", default="A;B;AB".split(";"))

args = parser.parse_args()

vdjdb_path = args.vdjdb_path
distance_type = args.distance_type
epitope_col = args.epitope_col
multiproc = args.multiproc
nrandom = args.nrandom
nlim = args.nlim
compare_methods_path = args.compare_methods_path
datasets = args.datasets
chains = args.chains
organisms = []
for i in args.organisms:
    if i in org_to_vdjdb_org:
        organisms.append(org_to_vdjdb_org[i])
    else:
        organisms.append(i)

crdir(os.path.join(compare_methods_path, 'compare'))

for organism in organisms:
    for dataset in datasets:
        if dataset == 'filtered':
            vdjdb_inp = '{}/vdjdb_imgted_filtered.txt'.format(vdjdb_path)
            output_dir = ('{}/compare/{}_{}_{}'.format(compare_methods_path, organism, 'filtered', distance_type))

        elif dataset == 'all':
            vdjdb_inp = '{}/vdjdb_imgted.txt'.format(vdjdb_path)
            output_dir = ('{}/compare/{}_{}_{}'.format(compare_methods_path, organism, all, distance_type))

        else:
            vdjdb_inp = '{}/vdjdb_imgted_filtered_{}.txt'.format(vdjdb_path, dataset)
            output_dir = ('{}/compare/{}_{}_{}'.format(compare_methods_path, organism, dataset, distance_type))

        for chain in chains:
            print('START ch_compute_distances_multiprocess.py; organism: {}; dataset: {}; chain: {}'.format(organism, dataset, chain))
            cmd = "python3 ch_compute_distances_multiprocess.py -clones_file {} -organism {} -chains {} -ch_distance_type {} -epitope_col {} -multiproc {}".format(vdjdb_inp, organism, chain, distance_type, epitope_col, multiproc)
            os.system(cmd)
            print('END ch_compute_distances_multiprocess.py\n')

            outtcrinfo = make_dirprefix(vdjdb_inp, 'dists') + '_{}_{}_nbrdists.tsv'.format(organism, chain)

            print('START ch_random_tcr_distances_multiprocess.py; organism: {}; dataset: {}; chain: {}'.format(organism, dataset, chain))
            cmd = "python3 ch_random_tcr_distances_multiprocess.py -clones_file {} -organism {} -nrandom {} -chains {} -ch_distance_type {} -epitope_col {} -multiproc {}".format(
                    vdjdb_inp, organism, nrandom, chain, distance_type, epitope_col, multiproc)
            os.system(cmd)
            print('END ch_random_tcr_distances_multiprocess.py\n')

            outrandominfo = make_dirprefix(vdjdb_inp, 'random_tcr_distances') + '_{}_{}_random_nbrdists.tsv'.format(organism, chain)

            crdir(output_dir)
            print('START roc_curves_tcrdist.R; organism: {}; dataset: {}; chain: {}'.format(organism, dataset, chain))
            cmd = "Rscript {}/roc_curves_tcrdist.R --tcr_distances {} --random_tcr_distances {} --organism {} --epitope_col {} --output_dir {} --output_file {} --script.dir {} --chains {} --nlim {}".format(
                compare_methods_path, outtcrinfo, outrandominfo, organism, epitope_col, output_dir, "{}_{}".format(dataset, chain), compare_methods_path, chain, nlim)
            os.system(cmd)
            print('END roc_curves_tcrdist.R\n\n')