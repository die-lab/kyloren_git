from Bio import SeqIO
import csv
import sys
import re
import statistics

fasta_file = sys.argv[1]

sequence = []

def sum_characters(strings):
    character_count = 0
    for string in strings:
        character_count += len(string)
    return character_count

def meida(strings):
	len_each = []	
	for string in strings:
		len_each.append(len(string))
	return statistics.median(len_each)	

for i in SeqIO.parse(fasta_file, 'fasta'):
         sequence.append(i.seq)

bases_dict={}
bases_type = ['A','T','C','G']

for base in bases_type: 
        bases_dict[base] = 0

for base in bases_type: 
        for curr_sequence in sequence:
                bases_dict[base] = bases_dict[base] + curr_sequence.count(base)

num_UTR = len(sequence)
len_UTR = sum_characters(sequence)
avg_UTR = len_UTR/num_UTR
median_UTR = meida(sequence)

#print('the number of nucleotide is ' + str(len_UTR))
#print('the avarage length is' + str(avg_UTR))
#print('the median length is' + str(median_UTR))
gcskew = round((bases_dict['G']-bases_dict['C'])/(bases_dict['G']+bases_dict['C']),5)
atskew = round((bases_dict['A']-bases_dict['T'])/(bases_dict['A']+bases_dict['T']),5)
#print('\n')
#print('the composition is ')
#print(str(bases_type))
composition = []
for base in bases_type:
        composition.append(round(bases_dict[base]/len_UTR,5))
#print(str(composition))

#first line for 3UTR and mito, second one for smith 
#final_list1 = [ fasta_file.split('_')[0].split('.')[0].lower() + '_' + 'forward', fasta_file.split('_')[1].lower().split('.')[0], str(num_UTR), str(len_UTR), str(avg_UTR), str(median_UTR), str(composition), str(gcskew), str(atskew)]
#final_list2 = [ fasta_file.split('_')[0].split('.')[0].lower() + '_' + 'reverse', fasta_file.split('_')[1].lower().split('.')[0], str(num_UTR), str(len_UTR), str(avg_UTR), str(median_UTR), str(composition), str(gcskew), str(atskew)]
final_list = [ fasta_file.split('_')[1].split('.')[0] + '_' + fasta_file.split('_')[0], 'smith', str(num_UTR), str(len_UTR), str(avg_UTR), str(median_UTR), str(composition), str(gcskew), str(atskew)]

#output1 = ';'.join(final_list1)
#output2 = ';'.join(final_list2)
output = ';'.join(final_list)

#print(output1)
#print(output2)
print(output)
