import csv
import sys
from Bio import SeqIO
from Bio.Seq import Seq

species = 'BaAt'
genome_file = 'BaAt_mit.fasta'
gff_file = 'tmp_BaAt.gff3'

gff = []
with open(gff_file) as list:
        for line in csv.reader(list, delimiter="\t"):
                gff.append(line)

allowed_value = ['tRNA']
gff_pass = []
for line in gff:
        if line[2] in allowed_value:
                gff_pass.append(line)

for i in SeqIO.parse(genome_file, 'fasta'):
        genome = i.seq

for line in gff_pass:
	name = line[8].split('product=')[1]
	name_file = name + '_' + line[3] + '_' + line[4]	
	start = int(line[3]) - 1
	end = int(line[4]) - 1 
	strand = line[6]
	sequence = genome[start:end]	
	if strand == '-':
		sequence = Seq.reverse_complement(sequence)
	with open(f'{name_file}.fa', "w") as pre:
        	pre.write('>' + name_file '\n')
                pre.write(sequence + '\n')
		pre.close()






