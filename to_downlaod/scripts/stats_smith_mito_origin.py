#find influence of type of region on smith origin 
import csv
from Bio import SeqIO

genome_file = 'DaRe_mit.fasta'
gff_file = 'DaRe_mit.gff3'
smith_file = 'forward_dare_tsa_22.txt'

tab = []
with open(smith_file) as list:
	for line in csv.reader(list, delimiter="\t"):
	        tab.append(line)

for i in SeqIO.parse(genome_file, 'fasta'):
	genome = i.seq

gff = []
with open(gff_file) as list:
	for line in csv.reader(list, delimiter="\t"):
	        gff.append(line)

allowed_value = ['CDS', 'rRNA', 'tRNA', 'D_loop']
gff_pass = []
for line in gff:
	if line[2] in allowed_value:
		gff_pass.append(line)


conv_dict = {'tRNA':'T', 'CDS':'C', 'rRNA':'R', 'D_loop':'D'}
anno_genome = []
for nuc in range(0,len(genome)):
	anno_genome.append('U')

for line in gff_pass:
	region_type = line[2]
	start = int(line[3])
	end = int(line[4])
	anno_genome[start:end+1] = conv_dict[region_type]*((end-start)+1)

anno_genome = ''.join(anno_genome)

#keep only names of the sequences
name = tab[::2]
type_smith = []

for i in name:
	type_smith.append(str(anno_genome[int(str(i).split('_')[2]):int(str(i).split('_')[3])]))

type_smith_join = ''.join(type_smith)

#bases_dict={}
bases_type = ['T','C','R','U']

bases_tot = [{},{},{}]

for base in bases_type: 
	bases_tot[0][base] = type_smith_join.count(base)
	bases_tot[1][base] = bases_tot[0][base]/(len(type_smith_join))
	bases_tot[2][base] = bases_tot[0][base]/anno_genome.count(base)

for line in bases_tot:
	print(line)
