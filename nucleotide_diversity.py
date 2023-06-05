from Bio import SeqIO
from Bio import Align
from datetime import datetime

startTime = datetime.now()

def compare_sequences(seq1, seq2):
    aligner = Align.PairwiseAligner()
    alignment = aligner.align(seq1, seq2)
    a_similairty = alignment.score / max(len(seq1),len(seq2))
    a_sequences = alignment.sequences
    return [a_similairty,a_sequences]

def rank_sequences(dataset1, dataset2):
    rankings = []
    for seq1 in dataset1:
        for seq2 in dataset2:
            similarity = compare_sequences(seq1, seq2)
            rankings.append((seq1, seq2, similarity))

    rankings.sort(key=lambda x: x[2], reverse=True)
    return rankings

# Example usage
dataset1 = ['ATCGA', 'AGTCA', 'CTGAT']
dataset2 = ['ATCGT', 'AGTCA', 'CTAAT']

hairpin_file = 'hairpin.fa'
hairpin_dict = {}
hairpin_sequences = SeqIO.parse(open(hairpin_file),'fasta')

for fasta in hairpin_sequences:
    hairpin_dict[str(fasta.seq)] = fasta.id

dataset1 = list(hairpin_dict.keys())[1:1000]

generated_file = 'reads_fasta.fa'
generated_dict = {}
generated_sequences = SeqIO.parse(open(generated_file),'fasta')

for fasta in generated_sequences:
    generated_dict[str(fasta.seq)] = fasta.id

dataset2 = list(generated_dict.keys())[1:1000]

#ranked_sequences = rank_sequences(dataset1, dataset2)

#for seq1, seq2, similarity in ranked_sequences:
#   print(f"Sequence 1: {seq1}\nSequence 2: {seq2}\nSimilarity: {similarity}\n")

max_scores = {}
for seq1 in dataset1:
        scores = []
        similarity_list = {}
        for seq2 in dataset2:
            similarity_curr = compare_sequences(seq1,seq2)
            #scores.append(similarity_curr[0])
            similarity_list[str(similarity_curr[1][1])] = (similarity_curr[0])

        max_scores[str(seq1)] = [max(similarity_list.values()), list(similarity_list.keys())[list(similarity_list.values()).index(max(similarity_list.values()))]]

with open('best_hit.txt','w') as best_file:
        for i in max_scores:
            best_file.write(str(i) + '\t' + str(max_scores[i][0]) + '\t' + max_scores[i][1] + '\n')

print(datetime.now() - startTime)
