import pandas as pd
import numpy as np
import ch_cdr3s_human
from ch_base import *


def read_tcr(filename, organism, chains, epitope_col):
    tcr_col = ['cdr1.alpha', 'cdr2.alpha', 'cdr2.5.alpha', 'cdr3.alpha',
               'cdr1.beta', 'cdr2.beta', 'cdr2.5.beta', 'cdr3.beta',
               'v.alpha', 'v.beta', epitope_col, 'species']

    all_tcrs = pd.read_table(filename, usecols=tcr_col)
    all_tcrs = all_tcrs[all_tcrs['species'] == organism]

    for chain in chains:
        if chain == 'A':
            all_tcrs = all_tcrs[pd.notnull(all_tcrs['v.alpha'])]
            all_tcrs = all_tcrs[pd.notnull(all_tcrs['cdr3.alpha'])]
            all_tcrs = all_tcrs[all_tcrs['cdr3.alpha'].str.len() > 5]
        elif chain == 'B':
            all_tcrs = all_tcrs[pd.notnull(all_tcrs['v.beta'])]
            all_tcrs = all_tcrs[pd.notnull(all_tcrs['cdr3.beta'])]
            all_tcrs = all_tcrs[all_tcrs['cdr3.beta'].str.len() > 5]

    all_tcrs['v_alpha_rep'] = all_tcrs.loc[:, 'v.alpha'].map(
        ch_cdr3s_human.all_loopseq_representative[organism.lower()])
    all_tcrs['v_beta_rep'] = all_tcrs.loc[:, 'v.beta'].map(ch_cdr3s_human.all_loopseq_representative[organism.lower()])
    all_tcrs['tcr_info'] = list(zip(all_tcrs.v_alpha_rep.str.split(','), all_tcrs.v_beta_rep.str.split(','),
                                    all_tcrs['cdr3.alpha'], all_tcrs['cdr3.beta']))

    all_tcrs = all_tcrs.drop_duplicates(subset=['v_alpha_rep', 'v_beta_rep', 'cdr3.alpha', 'cdr3.beta', epitope_col],
                                        keep='first') #remove duplicates
    all_tcrs = all_tcrs.drop_duplicates(subset=['v_alpha_rep', 'v_beta_rep', 'cdr3.alpha', 'cdr3.beta'], keep=False) #remove crossreactivity

    all_tcrs = all_tcrs.reset_index(drop=True)

    return all_tcrs


def read_random_tcr(path_to_files, organism, chains, nrandom, seed=None):

    def all_loopseq_on_list(item):
        return ','.join(list(map(ch_cdr3s_human.all_loopseq_representative[organism.lower()].get, item.split(','))))


    random_tcrs = pd.DataFrame()
    for ab in 'AB':
        if ab not in chains:
            random_chain = pd.DataFrame(np.nan, index=range(0, nrandom),
                                        columns=['v_{}_rep'.format(ab.lower()), 'cdr3{}'.format(ab.lower())],
                                        dtype='str')
        else:
            random_chains_file = os.path.join(path_to_files, 'new_nextgen_chains_{}_{}.tsv'.format(orgtoorg[organism], ab))
            random_chains = pd.read_table(random_chains_file, usecols=['v_reps', 'j_reps', 'cdr3'])
            random_chains = random_chains[random_chains['cdr3'].str.len() > 5]

            random_chains.columns = ['v{}_reps'.format(ab.lower()), 'j{}_reps'.format(ab.lower()),
                                     'cdr3{}'.format(ab.lower())]
            random_chains['v_{}_rep'.format(ab.lower())] = random_chains.loc[:, 'v{}_reps'.format(ab.lower())].apply(
                all_loopseq_on_list, 1)
            random_chains = random_chains.drop_duplicates(subset=['v_{}_rep'.format(ab.lower()), 'cdr3{}'.format(ab.lower())], keep='first')

            random_chains = random_chains.reset_index(drop=True)

            if not seed is None:
                random_chain = random_chains.sample(frac=1, random_state=1).reset_index(drop=True)
            else:
                random_chain = random_chains.sample(frac=1).reset_index(drop=True)

        random_tcrs = pd.concat([random_tcrs, random_chain[:nrandom]], axis=1)

    random_tcrs['tcr_info'] = list(zip(random_tcrs.v_a_rep.str.split(','), random_tcrs.v_b_rep.str.split(','),
                                       random_tcrs['cdr3a'], random_tcrs['cdr3b']))

    return random_tcrs