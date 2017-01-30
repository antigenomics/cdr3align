import os, argparse
import ch_antigen_base as ch_w
import pandas as pd

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will get information from tcr annotation file and then write "
                                             "down pdbs with unique cdrs. \
Current version is %s" % curr_version)

parser.add_argument("-i", help="Path to tcr annotation file. Default is '../final.annotations.txt'",
                    required=False, dest="input_dir")
parser.add_argument("-osort", help="Output path for sorted annotation file. Default is '../sorted.annotations.txt'",
                    required=False, dest="sort_dir")
parser.add_argument("-o", help="Path to output directory. Default is '../seqdata/HomoSapiens/cdrs/'", required=False, dest='output_dir')
parser.add_argument("-org", help="Species name. Default is Homo_sapiens", required=False, dest="organism")

myargs = parser.parse_args()

if myargs.input_dir is None:
    myargs.input_dir = '../final.annotations.txt'
if myargs.sort_dir is None:
    myargs.sort_dir = '../sorted.annotations.txt'
if myargs.output_dir is None:
    myargs.output_dir = '../seqdata/HomoSapiens/cdrs/'
ch_w.crdir(myargs.output_dir)
if myargs.organism is None:
    myargs.organism = 'Homo_sapiens'

class Found(Exception): pass


def readannotations(inppath, outpath):
    def fx(x):
        if isinstance(x, str):
            return (x[2])
        else:
            print(x)
            return x

    inptable = pd.read_csv(inppath, sep='\t')
    inptable_pdb = inptable.groupby('pdb_id')
    inptable = inptable.dropna(axis=0)
    outtable_species = inptable.sort_values(by=['species', 'antigen_seq', 'pdb_id'])
    outtable_species['chain_allele'] = outtable_species.tcr_v_allele.apply(fx)
    outtable_species.to_csv(outpath, sep='\t', columns=['species', 'antigen_seq', 'chain_allele',
                                                                         'tcr_region', 'pdb_id', 'tcr_region_seq',
                                                                         'chain_tcr',
                                                                         'chain_antigen', 'tcr_region_start',
                                                                         'tcr_region_end'], index=False)

def writecdrs(inppath, outpath, organism=None):
    inptable = pd.read_csv(inppath, sep='\t')
    for species in inptable.groupby('species'):
        if species[0] == organism or None:
            for tcr in species[1].groupby('tcr_region'):
                chains = []
                iti = 0
                for chain in tcr[1].groupby('chain_allele'):
                    chains.append(chain[1].groupby('antigen_seq').filter(lambda x: len(x)>1))
                preoutinfo = pd.concat(chains)
                outinfo = preoutinfo.sort_values(by=['antigen_seq', 'pdb_id'])
                outinfo.to_csv(os.path.join(outpath, str(tcr[0]).lower() + '.annotation.txt'), sep='\t',
                           columns=['species', 'antigen_seq', 'chain_allele', 'pdb_id', 'tcr_region_seq',
                                    'chain_tcr', 'chain_antigen', 'tcr_region_start', 'tcr_region_end'],
                           index=False)


def writecdrs_2(inppath, outpath):
    all = []
    cdr = {}
    cdr['cdr1'] = set()
    cdr['cdr2'] = set()
    cdr['cdr3'] = set()
    cdr1_2 = []
    cdr2_3 = []
    cdr3_1 = []
    for i in range(1,4):
        with open(os.path.join(inppath, 'cdr'+str(i)+'.annotation.txt'), 'r') as inp:
            with open(os.path.join(inppath, 'cdr'+str(i)+'.annotation_unique.txt'), 'w') as out:
                cd = 'cdr'+str(i)
                antig_set = set()
                tcr_regions = []
                for line in inp:
                    out.writelines(line)
                    for line in inp:
                        lineinfo = line.strip().split()
                        if lineinfo[1] in antig_set:
                            if str(lineinfo[2])+'\t'+str(lineinfo[4]) not in tcr_regions:
                                out.writelines(line)
                                tcr_regions.append(str(lineinfo[2])+'\t'+str(lineinfo[4]))
                                cdr[cd].add(lineinfo[1] + '\t' + lineinfo[3])
                        else:
                            out.writelines(line)
                            antig_set.add(lineinfo[1])
                            tcr_regions = [str(lineinfo[2])+'\t'+str(lineinfo[4])]
                            cdr[cd].add(lineinfo[1]+'\t'+lineinfo[3])

    for i in cdr['cdr1']:
        if i in cdr['cdr2'] and i in cdr['cdr3']:
            all.append(i)
        elif i in cdr['cdr2']:
            cdr1_2.append(i)
        elif i in cdr['cdr3']:
            cdr3_1.append(i)
    for i in cdr['cdr2']:
        if i in cdr['cdr3'] and i not in cdr['cdr1']:
            cdr2_3.append(i)

    for i in all:
        cdr['cdr1'].remove(i)
        cdr['cdr2'].remove(i)
        cdr['cdr3'].remove(i)
    for i in cdr1_2:
        cdr['cdr1'].remove(i)
        cdr['cdr2'].remove(i)
    for i in cdr2_3:
        cdr['cdr2'].remove(i)
        cdr['cdr3'].remove(i)
    for i in cdr3_1:
        cdr['cdr3'].remove(i)
        cdr['cdr1'].remove(i)

    with open(outpath, 'w') as out:
        out.write('all\n'+'\n'.join(sorted(all))+'\n//\n')
        out.write('cdr1_2\n' + '\n'.join(sorted(cdr1_2)) + '\n//\n')
        out.write('cdr2_3\n' + '\n'.join(sorted(cdr2_3)) + '\n//\n')
        out.write('cdr3_2\n' + '\n'.join(sorted(cdr3_1)) + '\n//\n')
        out.write('cdr1\n' + '\n'.join(sorted(cdr['cdr1'])) + '\n//\n')
        out.write('cdr2\n' + '\n'.join(sorted(cdr['cdr2'])) + '\n//\n')
        out.write('cdr3\n' + '\n'.join(sorted(cdr['cdr3'])) + '\n//\n')


readannotations(myargs.input_dir, myargs.sort_dir)
writecdrs(myargs.sort_dir, myargs.output_dir, organism=myargs.organism)
writecdrs_2(myargs.output_dir, os.path.join(myargs.output_dir, "cdrs.annotation.txt"))

