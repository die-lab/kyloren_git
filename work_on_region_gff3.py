##working on the previous version, that was doing fine, but it was not inserting the pipe to highlight the miRNA in the final dG file

import csv
import sys
from Bio import SeqIO

limit = 100 #set a limit of size of the region, below it the whole region is fold
part = 50 #set the two tail to take from both ends of mirna, when above the limit size

#genome_file = 'sequence.fasta'
#gff_file = 'tmp.gff3'
#smith_file = 'list.txt'

mode = sys.argv[1]
genome_file = sys.argv[2]
gff_file = sys.argv[3]
smith_file = sys.argv[4]

#FUNCTION TO INSERT PIPE
def ipipe(a,s1,s2):
	b = a[:s1] + '|' + a[s1:s2] + '|' + a[s2:]
	return b

#OPEN FILES
tab = []
with open(smith_file) as list:
	for line in csv.reader(list, delimiter="\t"):
	        tab.append(line)

gff = []
with open(gff_file) as list:
	for line in csv.reader(list, delimiter="\t"):
	        gff.append(line)

allowed_value = ['CDS', 'rRNA', 'tRNA', 'D_loop']
gff_pass = []
for line in gff:
	if line[2] in allowed_value:
		gff_pass.append(line)

for i in SeqIO.parse(genome_file, 'fasta'):
	genome = i.seq

#WORK ON
long_end = int(gff_pass[len(gff_pass) - 1][4])
long_mito = long_end - int(gff[0][4])

if ( long_mito > 0):
	tmp_first_line = gff_pass[len(gff_pass)-1][:]
	tmp_first_line[3] = '1'
	tmp_first_line[4] = str(long_mito)
	gff_pass.insert(0,tmp_first_line)	
	gff_pass[len(gff_pass)-1][4] = int(gff_pass[len(gff_pass)-1][4]) - long_mito

selected_regions = {}
for line in tab:
	start = int(line[2])
	end = int(line[3])
	min_s_diff = 100000
	min_e_diff = 100000	
	for a in gff_pass:
		diff_start = start-int(a[3])
		diff_end = int(a[4])-end
		if (min_s_diff > diff_start) and (diff_start >= 0):
			min_s_diff = diff_start
			s_line_to_find = a
		if (min_e_diff > diff_end) and (diff_end >= 0):
			min_e_diff = diff_end
			e_line_to_find = a
	if (int(s_line_to_find[4]) < start) and (int(e_line_to_find[3]) > end):
				replace_line = e_line_to_find[:]
				replace_line[2] = 'UR'
				replace_line[3] = str(int(s_line_to_find[4])+1)
				replace_line[4] = str(int(e_line_to_find[3])-1)
				replace_line[5]	= str(s_line_to_find[5]) + ',' + str(e_line_to_find[5])			
				s_line_to_find = replace_line[:]
				e_line_to_find = replace_line[:]
	selected_regions[line[0]] = [s_line_to_find, e_line_to_find]

#remove the second entry if the two are the same. i.e the miRNA is mapping on a single region, not in between.
for line in selected_regions:
	if selected_regions[line][0] == selected_regions[line][1]:
		selected_regions[line] = selected_regions[line][0][:]

##warning: looking into file at the current state, some selected regions are considered as mapping in between two genes, while mapping at the begining or the end of a single tRNA
#gotta tolerate some bases of mismatch? in this way, we fold only the tRNA while the candidate miRNA is at the boundaries of it.

coor_dict = {}

if mode == 'create_pre':
	for line in selected_regions:
		if len(selected_regions[line]) > 2 :
			region_seq = genome[int(selected_regions[line][3])-1:int(selected_regions[line][4])]
			region_type = selected_regions[line][2]
			if region_type == 'tRNA' or len(region_seq) < limit:
				for tab_line in tab:
					if line in tab_line:
						name = str('_'.join(tab_line).replace('>',''))
						coor_dict[name] = [int(tab_line[2]), int(tab_line[3]), int(selected_regions[line][3])-1, int(selected_regions[line][4])]
						with open(f'{name}.pre.fa', "w") as pre:
							pre.write('_'.join(tab_line) + '\n')
							pre.write(str(region_seq) + '\n')
							pre.close()
			else:
				for tab_line in tab:
					if line in tab_line:
						start = int(tab_line[2]) - part
						end = int(tab_line[3]) + part
						if (start - int(selected_regions[line][3])) < 0: 
							start = int(selected_regions[line][3])
						if (int(selected_regions[line][4]) - end) < 0:
							end = int(selected_regions[line][4])
						region_seq = genome[start:end]		
						name = str('_'.join(tab_line).replace('>',''))	
						coor_dict[name] = [int(tab_line[2]), int(tab_line[3]), start, end]			
						with open(f'{name}.pre.fa', "w") as pre:
							pre.write('_'.join(tab_line) + '\n')
							pre.write(str(region_seq) + '\n')
							pre.close()
		else:
			for tab_line in tab:
				if line in tab_line:
					start = int(tab_line[2]) - part
					end = int(tab_line[3]) + part
					region_seq = genome[start:end]
					name = str('_'.join(tab_line).replace('>',''))				
					coor_dict[name] = [int(tab_line[2]), int(tab_line[3]), start, end]								
					with open(f'{name}.pre.fa', "w") as pre:
						pre.write('_'.join(tab_line) + '\n')
						pre.write(str(region_seq) + '\n')
						pre.close()
	with open('coordinate.txt', 'w') as coor:
		for line in coor_dict:
			name = str(line)
			a = str(coor_dict[line][0])
			b = str(coor_dict[line][1])
			c = str(coor_dict[line][2])
			d = str(coor_dict[line][3])
			coor.write(name + '\t' + a + '\t' + b + '\t' + c + '\t' + d + '\n')
		coor.close()		

if mode == 'append_anno':
	conv_dict = {'tRNA':'T', 'CDS':'C', 'rRNA':'R', 'D_loop':'D'}
	anno_genome = []
	for nuc in range(0,len(genome)):
		anno_genome.append('U')
	for line in gff_pass:
		region_type = line[2]
		start = int(line[3])
		end = int(line[4])
		anno_genome[start:end+1] = conv_dict[region_type]*((end-start)+1)
	#nuc_genome = [a for a in genome]
	anno_genome = ''.join(anno_genome)
	#append annotation line at the end of every dG file created
	for line in selected_regions:
		if len(selected_regions[line]) > 2 :
			region_seq = anno_genome[int(selected_regions[line][3]):int(selected_regions[line][4])+1]
			region_type = selected_regions[line][2]
			if region_type == 'tRNA' or len(region_seq) < limit:
				for tab_line in tab:
					if line in tab_line:
						name = str('_'.join(tab_line).replace('>',''))				
						with open(f'{name}.pre.dG', "a") as pre:
							pre.write(str(region_seq) + '\n')
							pre.close()
			else:
				for tab_line in tab:
					if line in tab_line:
						start = int(tab_line[2]) - part
						end = int(tab_line[3]) + part
						if (start - int(selected_regions[line][3])) < 0: 
							start = int(selected_regions[line][3])
						if (int(selected_regions[line][4]) - end) < 0:
							end = int(selected_regions[line][4])
						region_seq = anno_genome[start:end+1]		
						name = str('_'.join(tab_line).replace('>',''))				
						with open(f'{name}.pre.dG', "a") as pre:
							pre.write(str(region_seq) + '\n')
							pre.close()
		else:
			for tab_line in tab:
				if line in tab_line:
					start = int(tab_line[2]) - part
					end = int(tab_line[3]) + part
					region_seq = anno_genome[start:end+1]
					name = str('_'.join(tab_line).replace('>',''))				
					with open(f'{name}.pre.dG', "a") as pre:
						pre.write(str(region_seq) + '\n')
						pre.close()

	coor = []
	with open('coordinate.txt') as list:
		for line in csv.reader(list, delimiter = '\t'):
			coor.append(line)
	for line in coor:
		curr = []	
		with open(f'{line[0]}.pre.dG','w') as cluster:
			for row in cluster:				
				curr.append(row.strip('\n'))
		skip1 = int(line[1]) - int(line[3])
		skip2 = int(line[2]) - int(line[3])
			cluster.close()
		with open(f'{line[0]}.pre.dG','w') as cluster:
			cluster.write(curr[0] + '\n' + ipipe(curr[1], skip1, skip2) + '\n' + ipipe(curr[2], skip1, skip2) + '\n' + ipipe(curr[3], skip1, skip2) + '\n')
			cluster.close()
		




