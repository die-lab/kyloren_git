from Bio import SeqIO
import glob
import re
import pandas as pd

micromito_files= sorted(glob.glob('*.txt'))


#empty matrix
columns = rows = list(range(int(min(re.findall(r'\d+', str(micromito_files)))), int(max(re.findall(r'\d+', str(micromito_files))))+1))
df_duplication = pd.DataFrame(columns = columns, index = rows)
df_duplication_perc = pd.DataFrame(columns = columns, index = rows)
df_exclusion = pd.DataFrame(columns = columns, index = rows)
df_exclusion_perc = pd.DataFrame(columns = columns, index = rows)

for a in range(0,len(micromito_files)):	
	a_file = micromito_files[a]	
	duplication = [0] * len(micromito_files)	
	dup_percentage = [0] * len(micromito_files)
	exclusion = [0] * len(micromito_files)
	excl_percentage = [0] * len(micromito_files)
	for b in range(a,len(micromito_files)):
		b_file = micromito_files[b]
		sequence_a = []
		sequence_b = []
		for i in SeqIO.parse(a_file, 'fasta'):
			sequence_a.append(i.seq)
		for i in SeqIO.parse(b_file, 'fasta'):
			sequence_b.append(i.seq)
		count_string = []
		count_seq_a = []
		count_seq_a_dict = {}
		excluded_seq = []
		for seq_a in sequence_a:
			count = 0
			for seq_b in sequence_b:	
				if ( str(seq_b).find(str(seq_a)) >= 0 ):
					count += 1
			count_string.append(count)	
			if ( count <= 0 ):
				excluded_seq.append(str(seq_a))
			else:
				count_seq_a.append(count)
				count_seq_a_dict[seq_a] = count
		duplication[b] = len(count_seq_a) - len(count_seq_a_dict)
		dup_percentage[b] = round(duplication[b]/len(sequence_a),5)
		exclusion[b] = len(excluded_seq)	
		excl_percentage[b] = round(exclusion[b]/len(sequence_a),5)
		#print('excluded: ' + a_file + b_file + str(exclusion))
		#print('duplicated in ' + a_file + ': ' + str(duplication))
	df_duplication.loc[rows[a]] = duplication	
	df_duplication_perc.loc[rows[a]] = dup_percentage
	df_exclusion.loc[rows[a]] = exclusion
	df_exclusion_perc.loc[rows[a]] = excl_percentage	
dup_trend_perc = []
for i in rows:
	dup_trend_perc.append(df_duplication_perc[i])

print(dup_trend_perc)

print(df_exclusion)
print(df_exclusion_perc)
