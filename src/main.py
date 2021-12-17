from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import pandas as pd

def transform_file(name):
    f = open(name, 'r+')
    lines = f.readlines()
    f.seek(0)
    f.write('locus_tag' + '\t')
    f.write('name' + '\t')
    for i in range(8):
        f.write(str(i+1) + '\t')
    f.write('e_value' + '\t')
    f.write(str(9) + '\n')
    for line in lines:
        f.write(line)

    f.close()

def first_part():

    list_record = []

    for record in SeqIO.parse("scaffolds.fasta", "fasta"):
      list_record.append(record)

    for i in range(len(list_record)):
      list_record[i].annotations['molecule_type'] = 'DNA'

    SeqIO.write(list_record, "GENOME.gbk", "genbank")

    def protocol_fantom(s):
        i = 0
        locus_tag = ''
        start = ''
        finish = ''
        while s[i] == ' ':
            i += 1
        while s[i] != ' ':
            locus_tag += s[i]
            i += 1
        while s[i] == ' ':
            i += 1
        strand = s[i]
        i += 1
        while s[i] == ' ' or s[i] == '<' or s[i] == '>':
            i += 1
        while s[i] != ' ':
            start += s[i]
            i += 1
        while s[i] == ' ' or s[i] == '<' or s[i] == '>':
            i += 1
        while s[i] != ' ':
            finish += s[i]
            i += 1
        return [locus_tag, strand, start, finish]



    f = open('gms2.lst', 'r')
    mas_of_s = []
    for i in f:
        mas_of_s.append(i[:-2])
    flag = False
    fitst_flag = True
    big_list = []
    result_in_one_scaf = []
    ctr_s = 'SequenceID'
    ctr_f = 'total_logodd'
    for i in mas_of_s:
        if ctr_f in i:
            flag = False
            big_list.append({scaf_index : result_in_one_scaf})
            result_in_one_scaf = []
        if flag:
            result_in_one_scaf.append(protocol_fantom(i))
        if ctr_s in i:
            flag = True
            fitst_flag = False
            k = i.find('_')
            try:
                scaf_index = int(i[k-2:k])
            except ValueError:
                scaf_index = int(i[k - 1])
            #print(scaf_index)


    #print(big_list[0][1])
    #print(len(big_list))
    #print(len(big_list[1][2]))
    #

    big_list_of_features = []

    for i in range(len(big_list)):
        list_of_features = []
        for j in range(len(big_list[i][i+1])):
            locus_tag = big_list[i][i+1][j][0]
            if big_list[i][i+1][j][1] == '-':
                strand = -1
            if big_list[i][i+1][j][1] == '+':
                strand = 1
            start = int(big_list[i][i+1][j][2])
            finish = int(big_list[i][i+1][j][3])
            #print(locus_tag, strand, start, finish)
            sc = SeqFeature(FeatureLocation(start, finish, strand=strand), type="CDS")
            sc.qualifiers['locus_tag'] = [locus_tag]
            list_of_features.append(sc)
        big_list_of_features.append({i + 1: list_of_features})


    for i in range(len(list_record)):
        list_record[i].features = big_list_of_features[i][i+1]

    SeqIO.write(list_record, "GENOME.gbk", "genbank")

    for record in SeqIO.parse("proteins.fasta", 'fasta'):
        for i in range(len(big_list_of_features)):
            for j in big_list_of_features[i][i+1]:
                if record.id == j.qualifiers['locus_tag'][0]:
                    j.qualifiers['translation'] = [record.seq]



    return list_record, big_list_of_features

list_record, big_list_of_features = first_part()
SeqIO.write(list_record, "GENOME.gbk", "genbank")

def second_part(list_record, big_list_of_features):
    df = pd.read_csv('scaffolds.hits_from_MIL_1.txt', sep='\t')

    list_of_names = {}
    for i in range(max(df['locus_tag'])):
        gene = df.loc[df['locus_tag'] == i+1]
        try:
            tru_gene = gene.loc[df['e_value'] == min(gene['e_value'])]
            name = list(tru_gene['name'])
            name = name[0]
            name = name.split('_')
            name = name[2]
            list_of_names[i+1] = name
        except ValueError:
            pass

    #print(list_of_names)

    #print(list_of_names[1])

    for i in range(len(big_list_of_features)):
        for j in big_list_of_features[i][i+1]:
            key = int(j.qualifiers['locus_tag'][0])
            try:
                j.qualifiers['note'] = list_of_names[key]
            except KeyError:
                pass
    #
    #
    #for j in big_list_of_features[1][2]:
    #    try:
    #        print(j.qualifiers['note'])
    #    except KeyError:
    #        pass

    mil_1_genome = SeqIO.read("T_oleivorans_MIL_1.gbk", "genbank")


    for mil_f in mil_1_genome.features:
        if 'protein_id' not in mil_f.qualifiers:
            continue

        for i in range(len(big_list_of_features)):
            for j in big_list_of_features[i][i + 1]:
                try:
                    a = j.qualifiers['note']
                except KeyError:
                    continue
                b = mil_f.qualifiers['protein_id'][0]
                if j.qualifiers['note'] == mil_f.qualifiers['protein_id'][0]:
                    j.qualifiers['product'] = mil_f.qualifiers['product']

    SeqIO.write(list_record, "GENOME.gbk", "genbank")
    return list_record, big_list_of_features

name_f = 'scaffolds.hits_from_MIL_1.txt'
transform_file(name_f)

list_record, big_list_of_features = second_part(list_record, big_list_of_features)

name_f = 'scaffolds.hits_from_SwissProt.txt'
transform_file(name_f)

df = pd.read_csv(name_f, sep='\t')

list_of_names = {}
for i in range(max(df['locus_tag'])):
    gene = df.loc[df['locus_tag'] == i+1]
    try:
        tru_gene = gene.loc[df['e_value'] == min(gene['e_value'])]
        name = list(tru_gene['name'])
        name = name[0]
        name = name.split('|')
        name = name[2]
        list_of_names[i+1] = name
    except ValueError:
        pass

print(list_of_names)

for i in range(len(big_list_of_features)):
    for j in big_list_of_features[i][i + 1]:
        key = int(j.qualifiers['locus_tag'][0])
        try:
            j.qualifiers['note'] = list_of_names[key]
        except KeyError:
            pass

SeqIO.write(list_record, "GENOME.gbk", "genbank")