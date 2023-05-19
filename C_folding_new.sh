#!/bin/bash


########### OPTIONS #########################################
#Default
home=$PWD
RNAfoldT=25
organism="RuPh"
while getopts "T:O:" opt; do
case "$opt" in
	T)RNAfoldT=$OPTARG;;
        O)organism=$OPTARG;;
        \?) echo "Argument Error in command line"
        exit 1;;
esac
done
##############################################################

conda activate smithRNA_env2

out_mit=$home/$organism"_mit_doubled"/$organism"_mit.fasta"
gff_mit=$home/$organism"_mit.gff3"


#folder creation
if [ -d "11_"$organism"_folding" ]
        then rm -rf "11_"$organism"_folding"
        fi
mkdir "11_"$organism"_folding"

grep ">" "7_"$organism"_smithRNAs"/smithRNAs.fa | sed 's/_/\t/g' | sed 's/size//g' | sed 's/strand//g' | sed 's/pos//g' > "11_"$organism"_folding"/list.txt

cd "11_"$organism"_folding"

mkdir plots

grep -v '#' $gff_mit > tmp.gff3
python $home/scripts/work_on_region_gff3.py create_pre $out_mit tmp.gff3 list.txt

echo -e "Pre_smith_ID\tDG1\tGG2" > DG.txt

num_pre_smith=$(cat list.txt | wc -l)
pre_smith=1

while  [ $pre_smith -le $num_pre_smith ]; do
        cluster=$(head -n $pre_smith list.txt | cut -f 1 | tail -n 1 | sed 's/>//g')
        size=$(head -n $pre_smith list.txt | cut -f 2 | tail -n 1)
        start=$(head -n $pre_smith list.txt | cut -f 3 | tail -n 1)
        end=$(head -n $pre_smith list.txt | cut -f 4 | tail -n 1)
        strand=$(head -n $pre_smith list.txt | cut -f 5 | tail -n 1)
        ID=$(echo $cluster"_"$size"_"$start"_"$end"_"$strand)

	if [ $strand == "-" ];
		then	
			revseq $ID".pre.fa" -tag FALSE -outseq "RC_"$ID".pre.fa"
			RNAfold -T $RNAfoldT "RC_"$ID".pre.fa" > "RC_"$ID".pre.dG"

			sed -i "s/>.*/>$ID/g" $ID".pre.fa"
			sed -i "s/>.*/>$ID/g" "RC_"$ID".pre.fa"
			sed -i "s/>.*/>$ID/g" "RC_"$ID".pre.dG"

			RNAplot "RC_"$ID".pre.dG" -o svg

		else
			RNAfold -T $RNAfoldT $ID".pre.fa" > $ID".pre.dG"

			sed -i "s/>.*/>$ID/g" $ID".pre.fa"
			sed -i "s/>.*/>$ID/g" $ID".pre.dG"

			RNAplot $ID".pre.dG" -o svg
	fi

	pre_smith=$(( pre_smith + 1 ))
done

python $home/scripts/work_on_region_gff3.py append_anno $out_mit tmp.gff3 list.txt

mv *svg plots

conda deactivate 

exit
