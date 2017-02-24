import os
import glob
import networkx as nx
import pandas as pd
import csv
import json
import time
import copy
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as mcls

def crdir(path):
    if os.path.exists(path):
        return 0
    else:
        os.mkdir(path)

print(os.getcwd())




def dict_sequences_with_antigen(inpp, outt=None, var=1):
    if var == 1:
        sequences_with_antigen = {}
        with open(inpp, 'r') as inp:
            for line in inp:
                if line.startswith('>'):
                    antigen = line.strip().split('_')[0].replace('>', '')
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
                    out.write(i + '\n')
        return sequences_with_antigen
    elif var == 2:
        sequences_with_antigen = {}
        sequences_list = []
        with open(inpp, 'r') as inp:
            for line in inp:
                data = line.strip().split()
                sequences_with_antigen[data[0].replace('>', '')] = data[1:]
                sequences_list.append(data[0].replace('>', ''))
        return sequences_with_antigen, sequences_list
#sequences_with_antigen, sequences_list = dict_sequences_with_antigen('trb_sequences_all_teach_antigen', var = 2)
#crdir('matrices_from_alignment')
#for subs in range(4):
#    for dels in range(4):
#        for ins in range(4):
#            total = [subs+dels+ins]
#            for tot in total:
#                break
#subs, dels, ins, tot = 3, 3, 3, 9
#print(subs, dels, ins, tot)
#bigmatrix = pd.DataFrame(0, index = sequences_list, columns=sequences_list)
#bigmatrix.to_csv('test_bigmatrix')
#data = set()
#data.add(('CASSLVGQGARQPQHF','CASSLADRVNTEAFF'))
#print(data)
#for aaa in data:
#    print(bigmatrix[aaa[0]][aaa[1]])
#    print(bigmatrix['CASSLVGQGARQPQHF']['CASSLADRVNTEAFF'])

def bigmatrix_creator(folder, outfolder, subs, dels, ins, tot):
    bigmatrix = pd.DataFrame(0, index=sequences_list, columns=sequences_list)
    bigmatrix.to_csv(outfolder + '/bigmatrix_' + str(subs) + '_' + str(dels) + '_' + str(ins) + '_' + str(tot))
    #subs, dels, ins, tot = 3, 3, 3, 9
    with open(folder+'/out_teach_'+str(subs)+'_'+str(dels)+'_'+str(ins)+'_'+str(tot), 'r') as inp:
        pairs = set()
        i = 0
        for line in inp:
            aaa = line.strip().split()
            pairs.add((aaa[0].replace('-','').upper(),aaa[1].replace('-','').upper()))
        print(len(pairs))
        for aaa in pairs:
            i += 1
            bigmatrix[aaa[0]][aaa[1]] = 1
            #print(i)
        bigmatrix.to_csv(outfolder+'/bigmatrix_'+str(subs)+'_'+str(dels)+'_'+str(ins)+'_'+str(tot))

#print(len(sequences_with_antigen))
#print(len(sequences_list))

def sequences_edges_antigens(inp_list, inp_dict, outp):
    with open(outp, 'w') as out:
        out.write('id,antigen\n')
        for i in inp_list:
            out.write(i+','+'_'.join(inp_dict[i])+'\n')
#sequences_edges_antigens(sequences_list, sequences_with_antigen, 'trb_sequences_teach_antigens.csv')

def renamefiles(usertemplate, format):
    files = glob.glob(usertemplate)
    for i in files:
        os.rename(i, i+'.'+format)

def antigens_to_jgroups(sequences_with_antigen):
    antigens = {}
    l = 1
    for i in sequences_with_antigen:
        if len(sequences_with_antigen[i]) < 2:
            for k in sequences_with_antigen[i]:
                if k in antigens:
                    pass
                else:
                    antigens[k] = l
                    l+=1
        else:
            for k in sequences_with_antigen[i]:
                if k in antigens:
                    pass
                else:
                    antigens[k] = l
            l+=1
    return antigens
#antigens = antigens_to_jgroups(sequences_with_antigen)
#print(len(antigens))

def json_antigen(sequences_with_antigen, antigens, csv_file, outfile):
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            nodes = row[1:]
            for i in range(len(nodes)):
                for k in sequences_with_antigen[nodes[i]]:
                    antigen = k
                    break
                nodes[i] = {"id":nodes[i], "group":antigens[antigen]}
            links = []
            for row in reader:
                linkf = row[1:]
                source = row[0]
                for i in range(len(linkf)):
                    if int(linkf[i]) == 1:
                        links.append({"source": source, "target": nodes[i]["id"]})
    data = {"nodes":nodes, "links":links}
    with open(outfile, 'w') as out:
        json.dump(data, out)

#json_antigen(sequences_with_antigen, antigens, 'matrices_from_alignment/bigmatrix_2_3_3_8.txt', 'matrices_from_alignment/bigmatrix_2_3_3_8.json')
#jgraph.draw('matrices_from_alignment/bigmatrix_2_3_3_8.json', directed = False)

def crgraphs(subs, dels, ins, total):
    print(subs, dels, ins, total)
    os.system("Rscript AdaptiveImm.R bigmatrix_"+str(subs)+"_"+str(dels)+"_"+str(ins)+"_"+str(total))
    print('done')
    time.sleep(5)


class Sequence:
    """New type that has following attributes: seq (sequence), antigen (all antigens, to whom this sequence belongs),
inlinks and outlinks (links to other sequences). You can use following commands: add_link(new_sequence) and
outlinked(antigen, opt)"""
    def __init__(self, seq, antigen):
        self.seq = seq
        self.antigen = antigen
        self.inlinks = {}
        self.outlinks = set()
        for i in self.antigen:
            self.inlinks[i]=[]
    def add_link(self, sequence):
        """
Add inlink or outlink to sequence (depending of antigen for new sequence). Type of added sequence is Sequence.
However, by self.inlink or self.outlink you'll get only string.
        """
        inantigen = 0
        for i in sequence.antigen:
            if i in self.antigen:
                self.inlinks[i].append(sequence)
                inantigen += 1
        if inantigen == 0:
            self.outlinks.update([sequence])

    def outlinked(self, antig):
        """
This command will get all links to sequences, that do not belong to chosen antigen (by default), but do belong to other
sequence's antigens.
        """
        outlink = 0
        inlink = 0
        not_links = set(self.inlinks[antig])
        inlink += len(not_links)
        for i in self.inlinks:
            if i == antig:
                pass
            else:
                for k in self.inlinks[i]:
                    if k in not_links:
                        inlink += 1
                        pass
                    else:
                        outlink += 1
        return outlink

class AntigenGroup:
    """New type that has been created only to concentrate sequences belonging to one antigen. Following attributes are:
antigen, sequences, boundarysequences (they belong to more than one group), outlinks, inlinks, totoutlinks, totinlinks.
You can use following command: add_seq(sequence), if you want to add new sequence to group.
"""
    def __init__(self, antigen):
        self.antigen = antigen
        self.sequences = []
        self.boundarysequences = [] #they belong to more than one group
        self.outlinks = set()
        self.inlinks = set()
        self.totoutlinks = 0
        self.totinlinks = 0
    def add_seq(self, sequence):
        """Add new sequence to Antigengroup.
"""
        if len(sequence.antigen) == 1:
            if sequence not in self.sequences:
                self.sequences.append(sequence)
                self.outlinks.update(sequence.outlinks)
                self.totoutlinks += len(sequence.outlinks)
                self.inlinks.update(sequence.inlinks[self.antigen])#not all inlinks are from one antigen
                self.totinlinks += len(sequence.inlinks[self.antigen])
        elif len(sequence.antigen) > 1:
            if sequence not in self.sequences:
                self.boundarysequences.append(sequence)
                self.outlinks.update(sequence.outlinks)
                self.totoutlinks += len(sequence.outlinks)
                out = set()
                inn = set()
                inn.update(sequence.inlinks[self.antigen])
                self.totinlinks += len(sequence.inlinks[self.antigen])
                for i in sequence.antigen:
                    if i == self.antigen:
                        pass
                    else:
                        for k in sequence.inlinks[i]:
                            if k in inn:
                                pass
                                #self.totinlinks += 1
                            else:
                                out.update([k])
                                self.totoutlinks += 1
                self.outlinks.update(out)
                self.inlinks.update(inn)


def get_information_from_matrices(matrice, sequences, antigens):
    sequences_2 = copy.deepcopy(sequences)
    antigens_2 = copy.deepcopy(antigens)
    print(matrice)
    with open(matrice, 'r') as inp:
        reader = csv.reader(inp)
        for row in reader:
            nodes = row[1:]
            for row in reader:
                node = row[0]
                links = row[1:]
                for i in range(len(links)):
                    if int(links[i]) == 1:
                        sequences_2[node].add_link(sequences_2[nodes[i]])
                for k in sequences_with_antigen[node]:
                    antigens_2[k].add_seq(sequences_2[node])
    return sequences_2, antigens_2

def write_information_from_matrices(matrice, outputfile, sequences, antigens):
    sequences_2, antigens_2 = get_information_from_matrices(matrice, sequences, antigens)

    listantigens = []
    for k in antigens_2:
        listantigens.append([k, antigens_2[k]])
    listantigens.sort(key = lambda x: len(x[1].sequences))
    with open(outputfile, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['antigen', 'number of sequences', 'possible inlinks', 'real inlinks', 'possible outlinks',
                        'real outlinks', 'number of boseq', 'inlinks of boseq', 'outlinks of boseq',
                        'possibe boseq links to other antigens', 'real boseq links to other antigens'])
        for k in listantigens:
            wantigen = k[0]
            wseqnumber = len(antigens_2[wantigen].sequences)+len(antigens_2[wantigen].boundarysequences)
            wposinlinks = ((wseqnumber**2)-wseqnumber)
            wrealinlinks = (antigens_2[wantigen].totinlinks)
            wposoutlinks = wseqnumber*(len(sequences_list)-wseqnumber)
            wrealoutlinks = antigens_2[wantigen].totoutlinks
            boseq = antigens_2[wantigen].boundarysequences
            if len(boseq) > 0:
                wboseqnumber = len(antigens_2[wantigen].boundarysequences)
                boseqin = 0
                boseqout = 0
                posboseqother = 0
                realboseqother = 0
                for i in boseq:
                    realboseqother += i.outlinked(wantigen)
                    boseqout += len(i.outlinks)+realboseqother
                    for antig in i.antigen:
                        if antig == wantigen:
                            pass
                        else:
                            #print(antig)
                            posboseqother += len(antigens_2[antig].sequences)
                            for seq in antigens_2[antig].boundarysequences:  #...
                                if wantigen in seq.antigen:
                                    #print(wantigen, seq.antigen, 'pass')
                                    pass
                                else:
                                    #print(wantigen, seq.antigen, 'not pass')
                                    posboseqother += 1
                    boseqin += len(i.inlinks[wantigen])
                wboseqinlinks = boseqin
                wboseqoutlinks = boseqout
                wposboseqother = posboseqother
                wrealboseqother = realboseqother
            else:
                wboseqnumber, wboseqinlinks, wboseqoutlinks, wposboseqother, wrealboseqother = 0, 0, 0, 0, 0
            writer.writerow([wantigen, wseqnumber, wposinlinks, wrealinlinks, wposoutlinks, wrealoutlinks, wboseqnumber, wboseqinlinks, wboseqoutlinks, wposboseqother, wrealboseqother])

def get_seqs_antigens(sequences_with_antigen):
    sequences = {}
    antigens = {}
    for i in sequences_with_antigen:
        sequences[i] = Sequence(i, sequences_with_antigen[i])
        for k in sequences_with_antigen[i]:
            if k in antigens:
                pass
            else:
                antigens[k] = AntigenGroup(k)
    return(sequences, antigens)

folder = 'matrices_from_alignment'
folder2 = 'information_from_matrices'

def get_info_from_matrices(folder, folder2, glword):
    matrices = glob.glob(folder+'/'+glword+'*')
    matrices = list(map(lambda x: x.replace(folder+'/', '').replace('.txt', ''), matrices))
    crdir(folder2)
    info = {}
    for name in matrices:
        with open(folder2 + '/' + name + '.csv', 'r') as inp:
            reader = csv.reader(inp, delimiter='\t')
            info[name] = {'realin':0,'realout':0, 'realboseqout':0, 'realboseqother':0, 'posin':0, 'posout':0, 'realboseqin':0, 'posboseqother':0}
            for row in reader:
                for row in reader:
                    info[name]['realin'] += int(row[3])
                    info[name]['realout'] += int(row[5])
                    info[name]['realboseqout'] += int(row[8])
                    info[name]['realboseqother'] += int(row[10])
                    info[name]['posin'] += int(row[2])
                    info[name]['posout'] += int(row[4])
                    info[name]['realboseqin'] += int(row[7])
                    info[name]['posboseqother'] += int(row[9])
    sortedinfo = []
    for name in info:
        sortedinfo.append([name, info[name]])
    sortedinfo.sort(key = lambda x: x[0])
    with open('sum_of_info_from_matrices', 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['name', 'rpin', 'rpout', 'rpin-rpout', 'realin', 'realout', 'realboseqin', 'realboseqout',
                         'realboseqother', 'posin', 'posout', 'posboseqother'])
        for inf in sortedinfo:
            name = inf[0]
            writer.writerow([name, "{0:.5f}".format(float(info[name]['realin'])/(info[name]['posin']+info[name]['posout'])),
                             "{0:.5f}".format(float(info[name]['realout'])/(info[name]['posin']+info[name]['posout'])),
                             "{0:.5f}".format((float(info[name]['realin']) / (info[name]['posin']+info[name]['posout']))-(float(info[name]['realout']) / (info[name]['posin']+info[name]['posout']))),
                             info[name]['realin'], info[name]['realout'],info[name]['realboseqin'], info[name]['realboseqout'],
                             info[name]['realboseqother'], info[name]['posin'], info[name]['posout'], info[name]['posboseqother']])
#name = 'bigmatrix_3_3_3_9'
#sequences_2, antigens_2 = get_information_from_matrices(folder+'/'+name+'.txt', sequences, antigens)
#print(antigens_2['SGEGSFQPSQENP'].antigen)
#print(antigens_2['SGEGSFQPSQENP'].sequences)
#print(antigens_2['SGEGSFQPSQENP'].inlinks)
#for seq in antigens_2['SGEGSFQPSQENP'].sequences:
#    print('new_seq')
#    print(seq.seq)
#    print(seq.inlinks)

#sub = [4, 5, 6, 7]
#for i in sub:
#    bigmatrix_creator('sub_del_ins_tot', 'matrices_from_alignment', i, 0, 0, i)
#    write_information_from_matrices(folder+'/bigmatrix_'+str(i)+'_'+str(0)+'_'+str(0)+'_'+str(i), folder2+'/bigmatrix_'+str(i)+'_'+str(0)+'_'+str(0)+'_'+str(i), sequences, antigens)

def crgraph(nodess, edgess, seqwantigen):
    nodes = []
    with open(nodess, 'r') as inp:
        for line in inp:
            nodes.append(line.strip())

    edges = []
    with open(edgess, 'r') as inp:
        for line in inp:
            connodes = line.strip().split()
            edges.append((connodes[0], connodes[1]))

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    antigen_with_sequences = {}
    for i in seqwantigen:
        for k in seqwantigen[i]:
            if k in antigen_with_sequences:
                antigen_with_sequences[k].append(i)
            else:
                antigen_with_sequences[k] = [i]
    print(len(antigen_with_sequences))

    colors = {}
    iti = 0
    for i in antigen_with_sequences:
        colors[i] = iti
        iti+=10
    node_colors = []
    for i in G.nodes():
        if len(seqwantigen[i]) == 1:
            node_colors.append(colors[sequences_with_antigen[i][0]])
        else:
            node_colors.append(iti+50)
    #for i in range(len(node_colors)):
    #    node_colors[i] = mcls.colorConverter.to_rgb(node_colors[i])

    edge_colors = []
    edge_weights = []
    for i in G.edges():
        if len(seqwantigen[i[0]]) == 1 and len(seqwantigen[i[1]]) == 1:
            if seqwantigen[i[0]][0] == seqwantigen[i[1]][0]:
                edge_colors.append(300)
                edge_weights.append(1)
            else:
                edge_colors.append(400)
                edge_weights.append(5)
        else:
            sameantigen = 0
            for k in seqwantigen[i[0]]:
                for l in seqwantigen[i[1]]:
                    if k == l:
                        sameantigen = 1
            if sameantigen == 1:
                edge_colors.append(300)
                edge_weights.append(1)
            else:
                edge_colors.append(400)
                edge_weights.append(5)

    #pos = nx.fruchterman_reingold_layout(G)
    print(len(G.nodes()))
    print(len(node_colors))
    print(len(G.edges()))
    print(len(edge_colors))
    plt.figure(1, figsize=(10, 10))
    nx.draw(G, pos=nx.fruchterman_reingold_layout(G, k=0.1),
            node_color=node_colors,
            cmap=plt.cm.jet,
            node_size=10,
            edge_color=edge_colors,
            edge_cmap=plt.cm.RdBu,
            edge_weight=edge_weights)
    plt.savefig('testaaa.png', dpi=300)

def groovy_treesearch(inp, opt, out):
    os.system("groovy TreeSearch.groovy " + inp + " " + opt + " " + out)