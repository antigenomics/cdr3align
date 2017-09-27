"""
This file was exported from tcr-dist package with some modifications (https://github.com/phbradley/tcr-dist)

MIT License

Copyright (c) 2017 Philip Harlan Bradley and Jeremy Chase Crawford

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from ch_base import *
import ch_read_blosum as blosum
import argparse
import random
import numpy as np
import csv
import pandas as pd
import ch_read_files
from multiprocessing import Pool
from datetime import datetime

#==============================

curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will compute random distances between all "
                                             "the tcrs from db and randomly created tcrs and write out a distances matrix. "
                                             "The order of rows and columns"
                                             "in the distance matrix will match that in the input file."
                                             "It will also compute the rank scores for each tcr with respect to all the "
                                             "epitopes present")

parser.add_argument("-clones_file", type=str, default="../../imgt_work/vdjdb_imgted_filtered_P.T.txt",
                    help="input file path.")
parser.add_argument("-organism", type=str, default='MusMusculus')
parser.add_argument("-nrandom", type=int, default=50)
parser.add_argument("-chains", type=str, default='A',
                    help="chains.")
parser.add_argument("-nbrdist_percentiles", type=str, default='5;10;25',
                    help="nbrdist_percentiles.")
parser.add_argument("-ch_distance_type", type=str, default='default')
parser.add_argument("-matrix", type=str, help='path to matrix')
parser.add_argument("-distance_params", type=str,
                    default='gap_penalty_v_region:4,gap_penalty_cdr3_region:8,weight_cdr3_region:3,align_cdr3s:False,'
                            'trim_cdr3s:True,scale_factor:1.0439137134052388')
parser.add_argument("-epitope_col", type=str, default='epitope')
parser.add_argument("-constant_seed", type=bool)
parser.add_argument("-multiproc", type=int, default=2)

#work with args
args = parser.parse_args()
clones_file = args.clones_file
organism = args.organism
nrandom = args.nrandom
nbrdist_percentiles = list(map(int, args.nbrdist_percentiles.split(';')))
chains = args.chains
epitope_col = args.epitope_col
multiproc = args.multiproc

ch_distance_type = args.ch_distance_type
if ch_distance_type == 'default':
    sign = 1
    import ch_tcr_distances_default as ch_tcr_distances
else:
    sign = -1 #just little hack for distance computation. Biopython is looking for alignments with higheest score. So I made score <0 and dist = -score
    import ch_tcr_distances_custom as ch_tcr_distances

if not args.matrix:
    matrix = blosum.bsd4
else:
    matrix = blosum.compute_cdrmatrix(blosum.read_blosum(args.matrix), sign=sign)
distance_params = ch_tcr_distances.DistanceParams( distance_matrix = matrix, config_string = args.distance_params )

if args.constant_seed:
    random.seed(1)

#


rep_dists = ch_tcr_distances.compute_all_v_region_distances(organism, distance_params)

tcr_col = ['cdr1.alpha', 'cdr2.alpha', 'cdr2.5.alpha', 'cdr3.alpha',
           'cdr1.beta', 'cdr2.beta', 'cdr2.5.beta', 'cdr3.beta',
           'v.alpha', 'v.beta', epitope_col, 'species']

all_tcrs = ch_read_files.read_tcr(clones_file, organism, chains, epitope_col)

random_tcrs = ch_read_files.read_random_tcr('../db/', organism, chains, nrandom, seed=args.constant_seed)

random_tcrs.to_csv('{}_{}_{}_random_tcrs.tsv'.format(make_dirprefix(clones_file, 'all_tcrs'), organism, chains),
                   sep='\t', index=None)

outfields = ['va_reps','vb_reps','cdr3a','cdr3b']
print('compute_dists')

def compute_dists(rindex):
    random_tcr = random_tcrs.loc[rindex, ]
    edists = []
    for tindex, tcr in epigroup.iterrows():
        edists.append(ch_tcr_distances.compute_distance(tcr['tcr_info'], random_tcr['tcr_info'],
                                                        chains, distance_params, rep_dists=rep_dists))
    edists.sort()
    for nbrdist_percentile in nbrdist_percentiles:
        nbrdist = ch_tcr_distances.sort_and_compute_nbrdist_from_distances(edists, nbrdist_percentile,
                                                                           dont_sort=True)
        wtd_nbrdist = ch_tcr_distances.sort_and_compute_weighted_nbrdist_from_distances(edists,
                                                                                        nbrdist_percentile,
                                                                                        dont_sort=True)
        random_tcr['{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = '{:.3f}'.format(
            nbrdist)
        random_tcr['{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = '{:.3f}'.format(
            wtd_nbrdist)
    return random_tcr


print('Start computing distances')
numepitopes = len(pd.unique(all_tcrs[epitope_col]))
iti = 0
for epitope, epigroup in all_tcrs.groupby(epitope_col):
    if len(epigroup) < 2:
        iti += 1
        continue
    iti += 1
    print('epitope {}, {}/{}'.format(epitope, iti, numepitopes))
    for nbrdist_percentile in nbrdist_percentiles: #create columns
        random_tcrs['{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = ''
        random_tcrs['{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = ''
        outfields.append('{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile))
        outfields.append('{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile))

    pool = Pool(multiproc)
    results = pool.map(compute_dists, [i for i in range(len(random_tcrs.index))])
    pool.close()
    pool.join()
    random_tcrs = pd.concat(results, axis=1).T

print('Computation of distances is completed')
random_tcrs.to_csv('{}_{}_{}_random_nbrdists.tsv'.format(make_dirprefix(clones_file, 'random_tcr_distances'), organism, chains),
                   sep='\t', index=None, columns=outfields)