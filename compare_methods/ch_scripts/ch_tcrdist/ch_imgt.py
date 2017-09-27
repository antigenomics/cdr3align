from Bio import SeqIO

#http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
anchors = {'CDR1': [26, 38], 'CDR2': [55, 65], 'CDR2.5': [80, 86]}

def get_vdj(input_file):
    dictimgt = {}

    with open(input_file, "r") as inp:
        for record in SeqIO.parse(inp, "fasta"):
            info = record.description.strip().split('|')
            specinfo = info[2].replace(' ', '').lower()
            if specinfo not in dictimgt:
                dictimgt[specinfo] = {}
            if info[1] in dictimgt[specinfo] and record.seq != dictimgt[specinfo][info[1]]:
                print(record.seq, info[1], dictimgt[specinfo][info[1]])
            if info[1].startswith('TRGC') or info[1].startswith('TRAC') or info[1].startswith('TRBC') or info[1].startswith('TRDC'):
                continue
            dictimgt[specinfo][info[1]] = record.seq
    return dictimgt

def get_specific_vdj(input_file, gene='V'):
    dictimgt = {}

    with open(input_file, "r") as inp:
        for record in SeqIO.parse(inp, "fasta"):
            info = record.description.strip().split('|')
            specinfo = info[2].replace(' ', '').lower()
            if specinfo not in dictimgt:
                dictimgt[specinfo] = {}
            if info[1] in dictimgt[specinfo] and record.seq != dictimgt[specinfo][info[1]]:
                print(record.seq, info[1], dictimgt[specinfo][info[1]])
            if info[1].startswith('TRB{}'.format(gene)) or info[1].startswith('TRA{}'.format(gene)):
                dictimgt[specinfo][info[1]] = record.seq
    return dictimgt

def changeseq(inpseq, positions):
    outdict, outdict2 = {}, {}
    for i in positions:
        outseq = str(inpseq[positions[i][0]:positions[i][1]])
        outdict[i] = outseq
        outdict2[i] = outseq.replace('.', '')
    return outdict, outdict2