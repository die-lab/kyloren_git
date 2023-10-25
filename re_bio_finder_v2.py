from Bio import Restriction
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter
from statistics import mean

### pattern to meet, depends on species order given to orthofinder ###
ro_excl = [1,1,1,2,2,2,2,2,2]
gr_excl = [1,1,1,2,2,2,2,1,1]
at_excl = [1,1,1,1,1,1,1,2,2]
all_excl = [1,1,1,2,2,2,2,3,3]
#wanted_pattern = [ro_excl,gr_excl,at_excl,all_excl]

#for rrna analyses
diciottos_alex = [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
ventottos_alex = [1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2]
wanted_pattern = [ventottos_alex]

#check if current match follow given code(s)
def match_comp(matches):
    matches_past = {}
    matches_check = []
    n = 1
    matches_past[str(matches[0])] = n
    for current_match in matches:
        if str(current_match) not in matches_past.keys():
            n = n + 1
            matches_past[str(current_match)] = n
        matches_check.append(matches_past[str(current_match)])
    return matches_check

#misure size of pieces, aftter cut by enzymes
def lengths_cut(match,len_record):
    match.append(len_record)
    len_pieces = []
    start = 0
    for cut in match:
        len_pieces.append(cut-start)
        start = cut
    return len_pieces

#sort a dataset based on number of cuts, length of enzymes (6 preferred)
def sort_dataset(line):
    line_to_sort = []
    favorite_re_len = 0
    if len(line[0].site) == 6:
        favorite_re_len = 0.1
    avoid_same_cut = 0
    if len(line[2]) < 2:
        avoid_same_cut = 100
    line_to_sort.append([line, mean(line[2]) - favorite_re_len - check_diff_pieces(line)+ avoid_same_cut])
    #print(line_to_sort)
    return line_to_sort

#check if every element in a list is identical, to see if it cuts in (significantly) different position.
#note that, if cut in the same position, the code is not respected, so it will be not consider by match_comp function
def all_equal(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)

def check_diff_pieces(line):
    diff_lengths = []
    diff_score = 0
    for element in line[3]:
        diff_lengths.append([abs(j - i) for i, j in zip(element, element[1:])])
    diff_check = []
    for a in diff_lengths:
        diff_check.append(all(i >= 30 for i in a))
    if True in diff_check:
        diff_score = 0.05
    return diff_score


def dealwithalign_v1(align_file):
    for i in align_file:
        i.seq = i.seq.replace('-','')

def dealwithalign_v2(align_file):
    for i in align_file:
        i.seq = i.seq.replace('-','n')
        bobo = list(i.seq)
        for element in range(0, len(bobo)):
            if bobo[element] == 'n':
                if bobo[element-1] != 'n' and bobo[element+1] != 'n':
                    bobo[element] = ''
        i.seq = Seq(''.join(bobo))

### read enzymes list ###
list_re = []
for enzyme_name in Restriction.CommOnly:
    list_re.append(enzyme_name)

### read orthogroups ###
orthogroup = []
with open('28s_clean.fa') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        orthogroup.append(record)

dealwithalign_v2(orthogroup)

#####
workin_enzymes = []
for enzyme_name in Restriction.CommOnly:
    matches = []
    all_cuts = []
    current_cuts = []
    for record in orthogroup:
        match = enzyme_name.search(record.seq)
        current_cuts.append(len(match))
        all_cuts.append(lengths_cut(match,len(record)))
        matches.append(match)
    check_on_current = match_comp(matches)
    if check_on_current in wanted_pattern:
        #workin_enzymes.append(enzyme_name)
        workin_enzymes.append([enzyme_name, check_on_current, set(current_cuts), set(tuple(row) for row in all_cuts)])
        #print('{0}\t{1}\t{2}\t{3}'.format(str(enzyme_name), str(check_on_current), str(set(current_cuts)),str(set(tuple(row) for row in all_cuts))))

dataset_to_sort = []
for row in workin_enzymes:
    dataset_to_sort.append(sort_dataset(row)[0])
enzymes_to_buy = [x for x in sorted(dataset_to_sort, key=itemgetter(1)) if x[1] < 90]

for enzyme in enzymes_to_buy:
    print(str(enzyme[0][:1]) + '\t' + str(enzyme[0][2:]) + '\t' + str(round(enzyme[1],2)))
