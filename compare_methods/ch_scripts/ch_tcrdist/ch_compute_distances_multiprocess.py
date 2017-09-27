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
## the BLOSUM62 scoring matrix
##
## Amino acid substitution matrices from protein blocks.
## Henikoff S, Henikoff JG.
## Proc Natl Acad Sci U S A. 1992 Nov 15;89(22):10915-9.
## PMID: 1438297

from ch_base import *
import pandas as pd
import argparse
import ch_read_blosum as blosum
import ch_read_files
from multiprocessing import Pool, Manager


#==============================

curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will compute distances between all "
                                             "the tcrs and write out a distances matrix. The order of rows and columns"
                                             "in the distance matrix will match that in the input file."
                                             "It will also compute the rank scores for each tcr with respect to all the "
                                             "epitopes present")

parser.add_argument("-clones_file", type=str, default="../../imgt_work/vdjdb_imgted_filtered_P.T.txt",
                    help="input file path.")
parser.add_argument("-distfile_prefix", type=str)
parser.add_argument("-outfile", type=str)
parser.add_argument("-organism", type=str, default='MusMusculus')
parser.add_argument("-intrasubject_nbrdists", type=bool, default=2)
parser.add_argument("-clones_files", type=str)
parser.add_argument("-epitope_prefixes", type=str)
parser.add_argument("-nbrdist_percentiles", type=str, default='5;10;25',
                    help="nbrdist_percentiles.")
parser.add_argument("-chains", type=str, default='B',
                    help="chains.")
parser.add_argument("-ch_distance_type", type=str, default='default')
parser.add_argument("-matrix", type=str, help='path to matrix')
parser.add_argument("-distance_params", type=str,
                    default='gap_penalty_v_region:4,gap_penalty_cdr3_region:8,weight_cdr3_region:3,align_cdr3s:False,'
                            'trim_cdr3s:True,scale_factor:1.0439137134052388')
parser.add_argument("-epitope_col", type=str, default='epitope')
parser.add_argument("-multiproc", type=int, default=2)

#work with args
args = parser.parse_args()
clones_file = args.clones_file
organism = args.organism
intrasubject_nbrdists = args.intrasubject_nbrdists
new_nbrdists = not intrasubject_nbrdists
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

if not args.clones_files:
    clones_files=[clones_file]
else:
    clones_files=args.clones_files

if not args.epitope_prefixes:
    epitope_prefixes = ['']*len(clones_files)
else:
    epitope_prefixes = args.epitope_prefixes

if not args.outfile:
    outfile = '{}_{}_{}_nbrdists.tsv'.format(make_dirprefix(clones_files[0], 'dists'), organism, chains)
else:
    outfile = args.outfile

if not args.distfile_prefix:
    distfile_prefix = distfile_prefix = make_dirprefix(clones_files[0], 'dists')
else:
    distfile_prefix = args.distfile_prefix
#

rep_dists = ch_tcr_distances.compute_all_v_region_distances(organism, distance_params)

tcr_col = ['cdr1.alpha', 'cdr2.alpha', 'cdr2.5.alpha', 'cdr3.alpha',
           'cdr1.beta', 'cdr2.beta', 'cdr2.5.beta', 'cdr3.beta',
           'v.alpha', 'v.beta', epitope_col, 'species']

all_tcrs = ch_read_files.read_tcr(clones_file, organism, chains, epitope_col)

all_tcrs.to_csv('{}_{}_{}_all_tcrs.tsv'.format(make_dirprefix(clones_files[0], 'all_tcrs'), organism, chains))

print(all_tcrs)

def compute_dists(index):
    ntcr = all_tcrs.loc[index,]
    edists = []
    for tindex, tcr in epigroup.iterrows():
        if tcr['tcr_info'] == ntcr['tcr_info']:
           continue
        edists.append(ch_tcr_distances.compute_distance(tcr['tcr_info'], ntcr['tcr_info'],
                                                        chains, distance_params, rep_dists=rep_dists))
    edists.sort()
    for nbrdist_percentile in nbrdist_percentiles:
        nbrdist = ch_tcr_distances.sort_and_compute_nbrdist_from_distances(edists, nbrdist_percentile,
                                                                           dont_sort=True)
        wtd_nbrdist = ch_tcr_distances.sort_and_compute_weighted_nbrdist_from_distances(edists,
                                                                                        nbrdist_percentile,
                                                                                        dont_sort=True)
        ntcr['{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = '{:.3f}'.format(
            nbrdist)
        ntcr['{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = '{:.3f}'.format(
            wtd_nbrdist)
    return ntcr


print('Start computing distances')
numepitopes = len(pd.unique(all_tcrs[epitope_col]))
iti = 0
for epitope, epigroup in all_tcrs.groupby(epitope_col):
    if len(epigroup) < 2:
        iti += 1
        continue
    iti += 1
    print('epitope {}, {}/{}'.format(epitope, iti, numepitopes))
    for nbrdist_percentile in nbrdist_percentiles:
        all_tcrs['{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = ''
        all_tcrs['{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = ''
        tcr_col.append('{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile))
        tcr_col.append('{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile))

    pool = Pool(multiproc)
    results = pool.map(compute_dists, [i for i in range(len(all_tcrs.index))])
    pool.close()
    pool.join()
    all_tcrs = pd.concat(results, axis=1).T


print('Computation of distances is completed')
all_tcrs.to_csv(outfile, sep='\t', index=None, columns=tcr_col)
