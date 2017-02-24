import os
import random
import copy

def crdir(path):
    if os.path.exists(path):
        return 0
    else:
        os.mkdir(path)

filename = 'trb_sequences_all'

seq_dict = {}
with open(filename, 'r') as inp:
    for line in inp:
        if line.startswith('>'):
            antigen = line.strip().split('_')[0].replace('>', '')
            if antigen in seq_dict:
                pass
            else:
                seq_dict[antigen] = []
        else:
            seq_dict[antigen].append(line.strip())
print(seq_dict)
seqs = []
for i in seq_dict:
    seqs.append([i, len(seq_dict[i])])
seqs.sort(key=lambda x: x[1])
seqs = seqs[::-1]
print(seqs)
'''
with open(filename+'_for_gelfand', 'w') as out:
    for i in seqs:
        out.write(i[0]+'\t'+str(i[1])+'\n')
'''
def test_and_teach(seq_dict, limit, times):
    """seq_dict – dictionary. key = antigen. limit – how many sequences should interact with specific antigen. times –
how size of test sample is lower than quantity of all sequences"""
    seq_dict_teach = copy.deepcopy(seq_dict)
    seq_dict_test = {}
    print('1. len(seq_dict_teach) = ', len(seq_dict_teach))
    for i in seq_dict_teach:
        if len(seq_dict_teach[i]) <= limit:
            pass
        else:
            print('2. len(seq_dict_teach[', i, ']) = ', len(seq_dict_teach[i]))
            seq_dict_test[i] = []
            for k in range(int(len(seq_dict_teach[i])/times)):
                number = random.randint(0, len(seq_dict_teach[i])-1)
                seq_dict_test[i].append(seq_dict_teach[i][number])
                seq_dict_teach[i].pop(number)
            print('3. len(seq_dict_teach[', i, ']) = ', len(seq_dict_teach[i]),'\n')
    with open(filename+'_test', 'w') as out:
        for i in seq_dict_test:
            pass
            out.write('\n'.join(seq_dict_test[i])+'\n')
    with open(filename+'_teach', 'w') as out:
        for i in seq_dict_teach:
            out.write('\n'.join(seq_dict_teach[i])+'\n')
    #for i in seq_dict_teach:
    #    print('5. len(seq_dict_teach[', i, ']) = ', len(seq_dict_teach[i]))

def dict_sequences_with_antigen(inpp, outt=None):
    """for each sequence we will get all antigens with what sequence can interact. In outt we can write it."""
    with open(inpp, 'r') as inp:
        sequences_with_antigen = {}
        for line in inp:
            if line.startswith('>'):
                antigen = line.strip().split('_')[0].replace('>','')
            else:
                if line.strip() in sequences_with_antigen:
                    sequences_with_antigen[line.strip()].append(antigen)
                else:
                    sequences_with_antigen[line.strip()] = [antigen]
    if outt == None:
        print('outt == None; output file won\'t be created')
        return sequences_with_antigen
    else:
        with open(outt, 'w') as out:
            for i in sequences_with_antigen:
                out.write(i+'\n')
    return sequences_with_antigen

def get_antigen_for_sequence(dict_with_all_seq, our_seqs):
    """for your antigen you can get sequences from dict with all seqs."""
    seq_with_antigen = {}
    for i in our_seqs:
        seq_with_antigen[i] = dict_with_all_seq[i]
    return seq_with_antigen


def get_seqs(inputt):
    """from input you can get sequences."""
    with open(inputt, 'r') as inp:
        our_seqs = []
        for line in inp:
            our_seqs.append(line.strip())
    return our_seqs

def write_seqs_and_antigen(inputt, output):
    """write >seq \t seq1, seq2, seqn"""
    with open(output, 'w') as out:
        for i in inputt:
            out.write('>'+i+'\t'+'\t'.join(inputt[i])+'\n')
            return 0

def write_halfs(inputt, outt):
    """for partition of test sample in halves."""
    if type(inputt) != list:
        raise TypeError(str(inputt) + ' not a list')
    half_1 = copy.deepcopy(inputt)
    half_2 = []
    for i in range(int(len(half_1)/2)):
        k = random.randint(0, len(half_1)-1)
        half_2.append(half_1[k])
        half_1.pop(k)
    print('len(half_1) =', len(half_1))
    print('len(half_2) =', len(half_2))
    with open(outt + '_first_testing', 'w') as out:
        out.write('\n'.join(half_1))
    with open(outt + '_sacred_testing', 'w') as out:
        out.write('\n'.join(half_2))


#test_and_teach(seq_dict, 16, 2)

#sequences_with_antigen = dict_sequences_with_antigen('trb_sequences_all')
#test_sample = get_seqs('trb_sequences_all_test')
#write_halfs(test_sample, 'trb_sequences_half')