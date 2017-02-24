import os
import pandas as pd
import numpy as np
import ch_antigen_base as ch_w
import ch_get_superimpose_structures as ch_superimpose
import Bio.PDB
import math

class Found(Exception): pass

tto = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'L', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}
#0 pdb_id 2f53
#1 species HomoSapiens
#2 chain_mhc_a A
#3 mhc_a_allele blabla
#4 chain_mhc_b B
#5 mhc_b_allele blabla
#6 mhc_type MHCI
#7 chain_antigen C
#8 antigen_seq !!!
#9 chain_tcr D
#10 tcr_v_allele TRAV
#11 tcr_j_allele TRAJ
#12 tcr_region CDR1
#13 tcr_region_start 27
#14 tcr_region_end 33
#15 tcr_region_seq DSAIYN
