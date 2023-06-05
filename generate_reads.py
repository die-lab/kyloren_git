from Bio import SeqIO
import os
import numpy as np
import random
import sys

# choose the method to use. Available are 'random' and 'slinding'
method = 'sliding'

genome_fasta_file = 'MuMu_mit.fasta'
#window = [int(sys.argv[1])]
window = [100]

# only if random method is choosen
reads_population = 10000
putative = 'putative.smithRNAs.fa'

reads_dict = []

for i in SeqIO.parse(genome_fasta_file, 'fasta'):
    genome = i.seq
    genome_rev = genome.reverse_complement()
genome_size = len(genome)


def delete_item(item):
    if os.path.exists(item):
        os.remove(item)
        print("The file " + item + " has been deleted successfully, a new file will replace it")
    else:
        print("The file " + item + " does not exist!, it will be created now")


if method == 'random':
    # create a vector from where randomly sample length of reads, based on putative smithRNAs (created by true (i.e. with clusters centroids) run of B.sh)
    if os.path.isfile(putative):
        window = []
        putative_sequences = SeqIO.parse(open(putative), 'fasta')
        for fasta in putative_sequences:
            window.append(len(str(fasta.seq)))


    def build_reads(pos, ref):
        read = ref[int(pos):(int(pos) + int(np.random.choice(window, 1)))]
        return read


    def reads_name(b, pop):
        first = random.randint(1, len(pop))
        second = random.randint(1, len(pop))
        third = random.randint(1, len(pop))
        name = str(b) + "_" + str(first) + "_" + str(second) + "_" + str(third)
        return name


    for i in range(int((reads_population) / 2)):
        start = int(np.random.choice(range(1, genome_size), size=1))
        read = build_reads(start, genome)
        reads_dict.append(read)

    for i in range(int((reads_population) / 2)):
        start = int(np.random.choice(range(1, genome_size), size=1))
        read = build_reads(start, genome_rev)
        reads_dict.append(read)

    np.random.shuffle(reads_dict)

    # write fasta file
    delete_item("reads_fasta.fa")
    with open("reads_fasta.fa", "a") as fasta:
        a = 0
        while a < reads_population:
            fasta.write(">smithRNA_C" + reads_name(a, np.array(range(1, reads_population))) + "\n")
            fasta.write(str(reads_dict[a]) + "\n")
            a += 1
        fasta.close()

else:
    def build_reads(pos, ref):
        read = ref[int(pos):(int(pos) + int(window[0]))].replace('T','U')
        return read


    def reads_name(b, pop):
        first = int(b)
        second = (int(b) + int(window[0]))
        third = random.randint(1, len(pop))
        name = str(b) + "_" + str(first) + "_" + str(second) + "_" + str(third)
        return name


    for i in range(int(genome_size - window[0])):
        read = build_reads(i, genome)
        reads_dict.append(read)

    # write fasta file
    delete_item("reads_fasta.fa")
    with open("reads_fasta.fa", "a") as fasta:
        a = 0
        while a < (genome_size - window[0]):
            fasta.write(">smithRNA_C" + reads_name(a, np.array(range(1, genome_size))) + "\n")
            fasta.write(str(reads_dict[a]) + "\n")
            a += 1
        fasta.close()

# write genome
delete_item("genome.fa")
with open("genome.fa", "a") as gen:
    gen.write(">test_genome" + "\n")
    gen.write(str(genome))
    gen.close()