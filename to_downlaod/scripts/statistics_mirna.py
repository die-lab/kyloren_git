from Bio import SeqIO
import csv
import sys
import re

mature_file = sys.argv[1]
length = int(re.findall(r'\d+',mature_file)[0])
output_file = mature_file.replace('.txt','.csv')
mature = []


for i in SeqIO.parse(mature_file, 'fasta'):
	 mature.append(i.seq)

bases_dict={}

for i in range(0,length):
	bases_dict[i] = 'X'

for curr_mature in mature:
	for index_base in range(0,len(curr_mature)):
		bases_dict[index_base] = bases_dict[index_base] + curr_mature[index_base]

bases_type = ['A','T','C','G']

base_summary = []
for i in range(0,len(bases_dict)):     
    base_index=[]
    for bt in bases_type:
        base_index.append(bases_dict[i].count(bt))
    base_summary.append(base_index)

for i in range(0,len(base_summary)):
    print(*base_summary[i])
transposed = zip(*base_summary)

with open(output_file, 'w', newline='') as file:
	writer = csv.writer(file)
	for row in transposed:
		writer.writerow(row)
