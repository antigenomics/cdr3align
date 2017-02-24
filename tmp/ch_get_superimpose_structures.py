import os
import ch_antigen_base as ch_w
import Bio.PDB

class Found(Exception): pass


class NotDisordered(Bio.PDB.Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == 'A'

#{antigen -> structure -> chain -> other info}

def get_pdbs(inppath, outpath):
    pdbl = Bio.PDB.PDBList()
    pdbdict = {}
    with open(inppath, 'r') as inp:
        for line in inp:
            for line in inp:
                info = line.strip().split()
                if info[1] in pdbdict:
                    if info[3] in pdbdict[info[1]]:
                        pdbdict[info[1]][info[3]][info[2]] = {'CDR_seq': info[4], 'pdb_cdr': info[5],
                                                              'pdb_antigen': info[6], 'cdr_start': info[7],
                                                              'cdr_end': info[8]}
                    else:
                        pdbl.retrieve_pdb_file(info[3], pdir=os.path.join(outpath, info[1]))
                        pdbdict[info[1]][info[3]] = {info[2]: {'CDR_seq': info[4], 'pdb_cdr': info[5],
                                                               'pdb_antigen': info[6], 'cdr_start': info[7],
                                                               'cdr_end': info[8]}}
                else:
                    ch_w.crdir(os.path.join(outpath, info[1]))
                    pdbl.retrieve_pdb_file(info[3], pdir=os.path.join(outpath, info[1]))
                    pdbdict[info[1]] = {info[3]: {info[2]: {'CDR_seq': info[4], 'pdb_cdr': info[5],
                                                            'pdb_antigen': info[6], 'cdr_start': info[7],
                                                            'cdr_end': info[8]}}}
    return pdbdict
def superimpose(pdbdict, outpath, overwrite=True):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    for antigen in pdbdict:
        for idpdb1 in range(len(pdbdict[antigen])):
            for idpdb2 in range(idpdb1+1, len((pdbdict[antigen]))):
                pdb1 = sorted(pdbdict[antigen])[idpdb1]
                pdb2 = sorted(pdbdict[antigen])[idpdb2]
                if overwrite == False:
                    if os.path.exists(os.path.join(outpath, antigen, pdb2+"_align_"+pdb1+".pdb")):
                        continue
                for chain in pdbdict[antigen][pdb1]:
                    chain1 = pdbdict[antigen][pdb1][chain]['pdb_antigen']
                    break
                for chain in pdbdict[antigen][pdb2]:
                    chain2 = pdbdict[antigen][pdb2][chain]['pdb_antigen']
                    break
                ref_structure = pdb_parser.get_structure("reference", os.path.join(outpath, antigen, "/pdb"+pdb1+".ent"))
                sample_structure = pdb_parser.get_structure("sample", os.path.join(outpath, antigen, "/pdb"+pdb2+".ent"))
                ref_model    = ref_structure[0]
                sample_model = sample_structure[0]
                ref_atoms = []
                sample_atoms = []
                for ref_res in ref_model[chain1]:
                    if ref_res.get_id()[0] == ' ':
                        ref_atoms.append(ref_res['CA'])
                for sample_res in sample_model[chain2]:
                    if sample_res.get_id()[0] == ' ':
                        sample_atoms.append(sample_res['CA'])
                super_imposer = Bio.PDB.Superimposer()
                super_imposer.set_atoms(ref_atoms, sample_atoms)
                super_imposer.apply(sample_model.get_atoms())
                io = Bio.PDB.PDBIO()
                io.set_structure(sample_structure)
                io.save(os.path.join(outpath, antigen, pdb2+"_align_"+pdb1+".pdb"), select=NotDisordered())

def superimpose_two(pdb1, pdb2, inppath, outpath):
    pdbdict = {}
    with open(inppath, 'r') as inp:
        for line in inp:
            for line in inp:
                info = line.strip().split()
                if info[3] == pdb1:
                    if info[1] not in pdbdict:
                        pdbdict[info[1]] = {info[3]: {info[2]: {'CDR_seq': info[4], 'pdb_cdr': info[5],
                                                                'pdb_antigen': info[6], 'cdr_start': info[7],
                                                                'cdr_end': info[8]}}}
                    elif info[3] not in pdbdict[info[1]]:
                        pdbdict[info[1]][info[3]] = {info[2]: {'CDR_seq': info[4], 'pdb_cdr': info[5],
                                                               'pdb_antigen': info[6], 'cdr_start': info[7],
                                                               'cdr_end': info[8]}}
                elif info[3] == pdb2:
                    if info[1] not in pdbdict:
                        pdbdict[info[1]] = {info[3]: {info[2]: {'CDR_seq': info[4], 'pdb_cdr': info[5],
                                                                'pdb_antigen': info[6], 'cdr_start': info[7],
                                                                'cdr_end': info[8]}}}
                    elif info[3] not in pdbdict[info[1]]:
                        pdbdict[info[1]][info[3]] = {info[2]: {'CDR_seq': info[4], 'pdb_cdr': info[5],
                                                               'pdb_antigen': info[6], 'cdr_start': info[7],
                                                               'cdr_end': info[8]}}
                if len(pdbdict) == 2:
                    break
            break
    superimpose(pdbdict, outpath)
    return pdbdict