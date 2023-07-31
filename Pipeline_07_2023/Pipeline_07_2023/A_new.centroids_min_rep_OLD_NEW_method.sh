#!/bin/bash

#need to copy manually in the current directory the nuclear and mitochondrial genomes of interest named <organism_nuc.fa> and <organism_mit.fasta> respectively
#need to copy manually in the current directory the folder with the scripts for graphs generation nemed <scripts>
#need to run with bash -i to fix conda activate bias


########### OPTIONS #########################################
#Default
home=$PWD
fastq_folder=$home/fastq
scripts=$home/scripts
organism="Unknown"
Trimming=SE
clusteringID=0.95
clustering_mode=NEW
stringency=0.50
threads=5
min_rep=1

while getopts "H:F:O:T:I:S:t:M:C:" opt; do
case "$opt" in
        H)home=$(echo $OPTARG | sed  "s/\/$//");;
        F)fastq_folder=$(echo $OPTARG | sed  "s/\/$//");;
        O)organism=$OPTARG;;
        T)Trimming=$OPTARG;;
        I)clusteringID=$OPTARG;;
        S)stringency=$OPTARG;;
        t)threads=$OPTARG;;
	M)min_rep=$OPTARG;;
	C)clustering_mode=$OPTARG;;
        \?) echo "Argument Error in command line"
        exit 1;;
esac
done

#controllo gli argomenti
echo ""
echo "Working directory=$home"
echo ""

if [ -d $fastq_folder ] ; then
        echo "Raw fastq files directory=$fastq_folder"
        echo ""
else
        echo "fastq files missing !"
        exit 1
fi

if [ -f $home/$organism"_mit.fasta" ] && [ -f $home/$organism"_nuc.fna" ]; then
        echo ""Organism=$organism""
        echo ""
else
        echo "fasta genomes missing!"
        exit 1
fi


if [ "$Trimming" == "PE" ] ; then
        echo "Trimming in PE mode"
        echo ""
elif
[ "$Trimming" == "SE" ] ; then
        echo "Trimming in SE mode"
        echo ""
else
        echo "Error in -T arguments"
        echo "Usage: $0 -H <working dir path> -F <fastq dir path> -O <organism ID> -T <trimming:PE-SE> -I <clustering identity:0-1> -S <stringency:0-1> -M <min replicates> -C <clustering mode:OLD-NEW> -t <threads:1-40>"
        exit 1
fi

if (( $(echo "$clusteringID >= 0 && $clusteringID <= 1" | bc -l) )); then
        echo "Clustering ID=$clusteringID"
        echo ""
else
        echo "Error in -I opion. Accepted values betweem 0 and 1"
        echo "Usage: $0 -H <working dir path> -F <fastq dir path> -O <organism ID> -T <trimming:PE-SE> -I <clustering identity:0-1> -S <stringency:0-1> -M <min replicates> -C <clustering mode:OLD-NEW> -t <threads:1-40>"
        exit 1
fi

if (( $(echo "$stringency >= 0 && $stringency <= 1" | bc -l) )); then
        echo "Stringency=$stringency"
        echo ""
else
        echo "Error in -S opion. Accepted values between 0 and 1"
        echo "Usage: $0 -H <working dir path> -F <fastq dir path> -O <organism ID> -T <trimming:PE-SE> -I <clustering identity:0-1> -S <stringency:0-1> -M <min replicates> -C <clustering mode:OLD-NEW> -t <threads:1-40>"
        exit 1
fi


if [[ $threads -gt 0 && $threads -le 40 ]]; then
        echo "threads=$threads"
        echo ""
else
        echo "Error in -t opion. Accepted values between 1 and 40"
        echo "Usage: $0 -H <working dir path> -F <fastq dir path> -O <organism ID> -T <trimming:PE-SE> -I <clustering identity:0-1> -S <stringency:0-1> -M <min replicates> -C <clustering mode:OLD-NEW> -t <threads:1-40>"
        exit 1
fi

num_samples=$(ls $fastq_folder/*fastq.gz | sed s"/\_[1|2].fastq.gz//g" | sort | uniq | wc -l);

if [[ $min_rep -le $num_samples ]]; then
	echo "Min rep=$min_rep"
	echo ""
else
	echo "Error in -M opion. Must be lower than $num_samples"
	echo "Usage: $0 -H <working dir path> -F <fastq dir path> -O <organism ID> -T <trimming:PE-SE> -I <clustering identity:0-1> -S <stringency:0-1> -M <min replicates> -C <clustering mode:OLD-NEW> -t <threads:1-40>"
	exit 1
fi

if [ "$clustering_mode" == "NEW" ] ; then
        echo "Clustering mode=NEW"
elif
[ "$clustering_mode" == "OLD" ] ; then
        echo "Clustering mode=OLD"
else
        echo "Error in -C arguments"
        echo "Usage: $0 -H <working dir path> -F <fastq dir path> -O <organism ID> -T <trimming:PE-SE> -I <clustering identity:0-1> -S <stringency:0-1> -M <min replicates> -C <clustering mode:OLD-NEW> -t <threads:1-40>"
        exit 1
fi


#Folder and files definition
#SRA_entries=""
#sex=male
genome_mit=$home/$organism"_mit_doubled"
#genome_mit_file=$home/$organism"_mit"_$sex".fasta"
genome_nuc=$home/$organism"_nuc"
trimming=$home/"1_"$organism"_trimmed"
fastqc=$home/"2_"$organism"_fastqc"
alignments=$home/"3_"$organism"_alignments"
BedFiles=$home/"4_"$organism"_COV_and_BedFiles"
Clustering=$home/"5_"$organism"_clustering"
Plots=$home/"6_"$organism"_plots"

########### BOWTIE INDEX BUILDING AND FOLDERS CREATION##################################

#conda enviroment activation
conda activate smithRNA_env2

echo "START BUILDING BOWTIE2 INDEXES AT $(date +%X)" > $home/smith.log
cat $home/smith.log
sleep 4

#create bowtie indexes of the two genomes provided
bowtie2-build -f $home/$organism"_mit.fasta" $genome_mit
bowtie2-build -f $home/$organism"_nuc.fna" $genome_nuc

#folder creation
if [ -d $genome_mit ]
        then rm -rf $genome_mit
        fi
mkdir $genome_mit

if [ -d $genome_nuc ]
        then rm -rf $genome_nuc
        fi
mkdir $genome_nuc

mv $genome_mit.* $genome_mit
mv $home/$organism"_mit.fasta" $genome_mit

if [ -f $genome_mit_file ]
        then mv $genome_mit_file $genome_mit
        fi
mv $genome_nuc.* $genome_nuc

if [ -f $genome_nuc_file ]
        then mv $genome_nuc_file $genome_nuc
        fi

if [ -d $trimming ]
        then rm -rf $trimming
        fi
mkdir $trimming

if [ -d $fastqc ]
        then rm -rf $fastqc
        fi
mkdir $fastqc

if [ -d $alignments ]
        then rm -rf $alignments
        fi
mkdir $alignments

if [ -d $BedFiles ]
        then rm -rf $BedFiles
        fi
mkdir $BedFiles

if [ -d $Clustering ]
        then rm -rf $Clustering
        fi
mkdir $Clustering

if [ -d $Plots ]
        then rm -rf $Plots
        fi
mkdir $Plots


#starting with smallRNA reads (download fastq)
#fastq-dump --defline-seq '@$sn[_cds     $rn]/$ri' --split-files $SRA_entries


cd $fastq_folder;

########### TRIMMING ###########################################
#In SE (OPTION S: cutadapt) o in PE (OPTION P: fastp). I fastq PE devono essere formattati: replicaX_1.fastq.gz, replicaX_2.fastq.gz. I fastq SE:replicaX_1.fastq.gz.

#trimming variables 
#for illumina
adapter_sequence='TGGAATTCTCGGGTGCCAAGG';
adapter_sequenceR2='TGGAATTCTCGGGTGCCAAGG';
#for NEBNext Small RNA library
#adapter_sequence='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC';
#adapter_sequenceR2='GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT';
samples=$(ls *fastq.gz | sed "s/\_[1|2].fastq.gz//g" | sort | uniq);
num_threads_trim=5

for sample in $samples; do 
	echo $sample >> $home/samples.txt;	
done

sed -i '/^$/d' $home/samples.txt

if [ "$Trimming" == "PE" ] ; then
        echo "" >> $home/smith.log
        echo "START TRIMMING IN PE MODE AT $(date +%X)" >> $home/smith.log
	sleep 4
	echo ""
        cat $home/smith.log

        for sample in $samples; do
                echo "fastp -c -x -g -a $adapter_sequence --adapter_sequence_r2 $adapter_sequenceR2 -w $num_threads_trim -i $fastq_folder/$sample"_1.fastq.gz" -I $fastq_folder/$sample"_2.fastq.gz" -o $trimming/$sample"_1.trimmed.fp.fastq" -O $trimming/$sample"_2.trimmed.fp.fastq"" >> $trimming/parallel_fp.sh
        done

        cd $trimming

        parallel -j $num_samples < parallel_fp.sh

        rm *2.trimmed.fp.fastq

elif [ "$Trimming" == "SE" ] ; then
        echo "" >> $home/smith.log
        echo "START TRIMMING IN SE MODE AT $(date +%X)" >> $home/smith.log
	sleep 4 
	echo ""
	cat $home/smith.log

        for sample in $samples; do
                echo "cutadapt  --match-read-wildcards --times 1 -e 0.1 -O 5 --cores $num_threads_trim --quality-cutoff 6 -m 15 --discard-untrimmed -a $adapter_sequence -o $trimming/$sample"_1.trimmed.ctdp.fastq" $fastq_folder/$sample"_1.fastq.gz"" >> $trimming/parallel_ctdp.sh
       done
	
	cd $trimming
	
	parallel -j $num_samples < parallel_ctdp.sh ;
fi

#Assegno il nome del campione a  ogni campione trimmato e lo stampo all'inizio della sequenza "dopo @" e creo il file parallel per i fastqc
for trimmed in $(ls *.trimmed*); do
        a=$(echo $trimmed | sed s"/.trimmed.*//g");
        sed -i '1~4 s/^@/@'$a'_/g' $trimmed;
        echo "fastqc -t $num_threads_trim -o $fastqc $trimmed" >> parallel_fastqc.sh
done

parallel -j $num_samples < parallel_fastqc.sh

gzip *.trimmed*


#############ALIGNMENTS#################################################################
echo "" >> $home/smith.log
echo "START ALIGNMENTS AT $(date +%X)" >> $home/smith.log
echo ""
cat $home/smith.log

for infile in $(ls *.trimmed* | sed "s/.trimmed.*//g"); do
        outfile=$infile\_on_MITO1.sam
        file=$infile\_MITO1_mapping.fastq
        outfile1=$infile\_on_MITO1.bam
        prefile1=$infile\_on_NUCL.sam
	prefile2=$infile\_on_NUCL.bam
        postfile1=$infile\_on_NUCL_unmapping.bam
        postfile2=$infile\_NUCL_unmapping.fastq
        postfile3=$infile\_on_MITO2.sam
        mito=$infile\_MitoUnique.bam

        #allineo reads su mitocondrio e faccio statistiche 
        bowtie2 -x $genome_mit/$organism"_mit_doubled" -U $infile*trimmed* -S $alignments/$outfile -N 1 -i C,1 -L 18 -p $threads --no-unal
        samtools sort -O BAM -o $alignments/$outfile1 $alignments/$outfile
        samtools flagstat $alignments/$outfile1 > $alignments/$infile.Log_MITO1.txt
        samtools depth -a $alignments/$outfile1 -o $BedFiles/$outfile1.cov

       #trasformo in fastq, le mappo sul nucleo, faccio le stats 
        samtools fastq  $alignments/$outfile1 > $alignments/$file
        bowtie2 -x $genome_nuc/$organism"_nuc" -q $alignments/$file -S $alignments/$prefile1 -i C,1 -L 22 -p $threads
        samtools sort -O BAM -o $alignments/$prefile2 $alignments/$prefile1
        samtools flagstat $alignments/$prefile2 > $alignments/$infile.Log_NUCL.txt

        #estraggo le reads che non mappano sul nucleo, le converto in fastq, le allineo nuovamente sul mitocondrio
	samtools view -b -f 4 $alignments/$prefile2 > $alignments/$postfile1
        samtools fastq $alignments/$postfile1 > $alignments/$postfile2
        bowtie2 -x $genome_mit/$organism"_mit_doubled" -q $alignments/$postfile2 -S $alignments/$postfile3 -N 1 -i C,1 -L 18 -p $threads
        samtools sort -O BAM -o  $alignments/$mito $alignments/$postfile3
        samtools flagstat $alignments/$mito > $alignments/$infile.Log_MITO2.txt
        samtools depth -a $alignments/$mito -o $BedFiles/$mito.cov

        #cancello allineamenti intermedi
	rm $alignments/*sam
	rm $alignments/*MITO1.bam $alignments/*NUCL.bam $alignments/*unmapping.bam
        rm $alignments/*fastq

done;

cd $alignments;

#genero files con statistiche con le reads mappanti sul mito, sul nucleo e uniche del mito (rispettivamente Mapped1, Mapped2 e Mapped3)
for i in $(ls *Log_MITO1.txt | sed 's/.Log_MITO1.txt//g'); do
        name=$i;
        mapped=$(grep "mapped" $i.Log_MITO1.txt | head -n 1);
        echo $name $mapped >> Log_MITO1.txt;
done;


for i in $(ls *Log_NUCL.txt | sed 's/.Log_NUCL.txt//g'); do
        name=$i;
        mapped=$(grep "mapped" $i.Log_NUCL.txt | head -n 1);
        echo $name $mapped >> Log_NUCL.txt
done;

for i in $(ls *Log_MITO2.txt | sed 's/.Log_MITO2.txt//g'); do
        name=$i;
        mapped=$(grep "mapped" $i.Log_MITO2.txt | head -n 1);
        echo $name $mapped >> Log_MITO2.txt;
done;

#creo files bed degli allineamenti sul MITO e li sposto in cartella BEdFiles
for infile in $(ls *Unique.bam | sed 's/.bam//g'); do
        bedtools bamtobed -i $infile.bam > $infile.bed
        bedtools sort -i $infile.bed > $infile'_sort.bed'
        rm *Unique.bed
done

mv $alignments/*bed $BedFiles

cd $Clustering

############CLUSTERING#################################
output=$organism"_results"
output_centroids=$output.centroids.fa
output_selected_centroids=$output.centroids.selected.fa
clusters_bed="5.2_"$output.clusters.bedfiles
clusters_fasta="5.1_"$output.clusters.fasta


#creo cartelle di OUT
if [ -f $output".fa" ]
        then rm $output".fa"
        fi

if [ -f $output_selected_centroids ]
        then rm $output_selected_centroids
        fi
if [ -d $clusters_bed ]
        then rm -rf $clusters_bed
        fi
mkdir $clusters_bed
if [ -d $clusters_fasta ]
        then rm -rf $clusters_fasta
        fi
mkdir $clusters_fasta


#estraggo le reads mappate sul MITO di tutti i campioni in un file fasta
for bam in $alignments/*Unique.bam; do 
        samtools fasta -n $bam >> $output".fa"
done


if [ "$clustering_mode" == "NEW" ] ; then
	echo "" >> $home/smith.log
	echo "START CLUSTERING NEW (%ID=$clusteringID) AT $(date +%X)" >> $home/smith.log
	echo ""
	cat $home/smith.log

	#sort manuale
	for i in $(grep ">" $output".fa"); do 
		ID=$i; seq=$(grep -A1 $i $output.fa | tail -n 1);
		echo -e "$(grep -w "$seq" $output.fa | wc -l)\t$ID" >> lista;
	done;

	sort -n -r -k 1 lista > lista_sorted

	for i in $(cut -f 2 lista_sorted); do
		grep -w -A 1 "$i" $output.fa >> $output.sorted.fa;
	done;

	sed -i 's/--//g' $output.sorted.fa
	sed -i '/^$/d'  $output.sorted.fa

	rm lista

	vsearch --cluster_smallmem $output".sorted.fa" --usersort --threads 12 --centroids $output"_centroids.fa" --sizeout --clusterout_id --clusterout_sort --consout $output"_consensus.fa" --id $clusteringID --profile clusters_genomecov.txt --msaout alignment-consensu.txt --uc vsearch_table.txt --relabel_keep --minseqlength 15 --clusters clusterino;

elif [ "$clustering_mode" == "OLD" ] ; then
	echo "" >> $home/smith.log
        echo "START CLUSTERING OLD (%ID=$clusteringID) AT $(date +%X)" >> $home/smith.log
        echo ""
        cat $home/smith.log
	
	vsearch --cluster_fast $output".fa" --threads 12 --centroids $output"_centroids.fa" --sizeout --clusterout_id --clusterout_sort --consout $output_"consensus.fa" --id $clusteringID --profile clusters_genomecov.txt --msaout alignment-consensu.txt --uc vsearch_table.txt --relabel_keep --minseqlength 15 --clusters clusterino
fi


########### CLUSTERING THRESHOLDS ################################################
echo "" >> $home/smith.log
echo "START CLUSTERS FILTERING (STRINGENCY=$stringency) AT $(date +%X)"  >> $home/smith.log
echo ""
cat $home/smith.log

mv clusterino* $clusters_fasta;

cd $clusters_fasta;

#definisco nome e numero di repliche
num_samples=$(cat $home/samples.txt | wc -l)

#creo file INFO.txt con: clusterID, reads totali nel cluster, reads di ogni replica nel cluster -> INFO.txt
echo -e "CLUSTER\tCLUSTER_SIZE\t$(paste -sd "\t" $home/samples.txt)" > INFO.txt

for cluster in clusterino*; do
	cluster_size=$(grep ">" $cluster | wc -l)
        echo -e "$cluster\t$cluster_size" >> File.tmp;
        n=1
        for sample in $(cat $home/samples.txt); do
                grep ">" $cluster | grep  "$sample" | wc -l >> $n.$sample.file.tmp;
                n=$(($n+1))
done
done

paste File.tmp  *file.tmp  >> INFO.txt

#definisco Threshold1: calcolo la dimensione ed il numero di clusters con dimensione unica ed estraggo la dimensione dell'Nesimo cluster (in base alla stringenza)
cut -f 2 INFO.txt | grep -v "CLUSTER" | sort -n -r| uniq > clusters_uniq.tmp
num_clusters_uniq=$(cut -f 2 INFO.txt | grep -v "CLUSTER" | sort | uniq | wc -l)
raw=$(echo "$num_clusters_uniq*(1-$stringency)" | bc | sed  "s/\..*//g")
sort -n -r clusters_uniq.tmp | head -n $raw | tail -n 1 > T1
echo ""
echo "Threshold1=$(cat T1)"

rm clusters_uniq.tmp


#definisco T2: come T1 ma per ogni replica
for sample in $(cat $home/samples.txt); do
	#prima devo trovare il numero della colonna corrispondente ad ogni replica in INFO.txt! -> column
        column=$(head -n 1 INFO.txt | sed -e 's/\t/\n/g' | grep -w -n "$sample" | sed s'/:/\t/g' | cut -f 1)
	cut -f $column INFO.txt | grep -v "$sample" | sort -n -r| uniq > $sample.clusters_uniq.tmp
        num_clusters_uniq_size=$(cut -f $column INFO.txt | grep -v "$sample" | sort | uniq | wc -l)
        raw=$(echo "$num_clusters_uniq_size*(1-$stringency)" | bc | sed  "s/\..*//g")
        sort -n -r $sample.clusters_uniq.tmp | head -n $raw | tail -n 1 > $sample.T2
        echo "Threshold2 for $sample =$(cat $sample.T2)"
done

rm *tmp

#Stampo i clusters che passano e non passano la T1 -> FILTER_T1.txt e la T2 -> FILTER_T2.txt
for cluster in clusterino*; do
        Cluster_size=$(grep -c ">" $cluster)
        T1=$(cat T1)
        if [ $Cluster_size -ge $T1 ]
        then
        echo -e "$cluster\t$Cluster_size\tPASS" >> FILTER_T1.txt
        else
        echo -e "$cluster\t$Cluster_size\tNOT_PASS" >> FILTER_T1.txt
        fi
        for sample in $(cat $home/samples.txt); do
                Sample_size=$(grep $sample $cluster | wc -l )
                T2=$(cat $sample.T2)
                if [ $Sample_size -ge $T2 ]
                then
                        echo -e "$cluster\t$sample\t$T2\t$Sample_size\tPASS" >> replicates_threshold2.txt
                else
                        echo -e "$cluster\t$sample\t$T2\t$Sample_size\tNOT_PASS" >> replicates_threshold2.txt
                fi
        done

        #Conto dal file replicates_threshold2.txt il numero di repliche che passano e non passano la threshold2
        N_rep_PASS=$(grep -w $cluster replicates_threshold2.txt | grep -w "PASS" | wc -l)
        N_rep_NOT_PASS=$(grep -w $cluster replicates_threshold2.txt | grep -w "NOT_PASS" | wc -l)

        #Se in un cluster, il numero di repliche che passa la threshold2 è >= N-1 il cluster è PASS -> Threshold2.txt
        echo -e "$cluster\t$N_rep_PASS\t$N_rep_NOT_PASS" | awk -v num=$min_rep '{ if ($2>=num) print $1"\t"$2"\t"$3"\t""PASS"; else print $1"\t"$2"\t"$3"\t""NOT_PASS"}' >> FILTER_T2.txt

done

rm replicates_threshold2.txt

#formatto le tabelle
sed -i '1s/^/CLUSTER\tSIZE\tFILTER\n/' FILTER_T1.txt
sed -i '1s/^/CLUSTER\tN_REP_PASS\tN_REP_NOT_PASS\tFILTER\n/' FILTER_T2.txt

#stampo i risultati della T1
echo ""
echo "CLUSTERING THRESHOLD 1:"

sleep 2

cat FILTER_T1.txt

#stampo i risultati della soglia 2 di clustering
echo ""
echo "CLUSTERING THRESHOLD 2:"

sleep 2

cat FILTER_T2.txt


##################ELIMINO CLUSTERS SOTTO SOGLIA#################

#Elimino clusters che non passano T1 e T2
for clusters_NOT_PASS in $(cat FILTER_T1.txt FILTER_T2.txt | grep -w "NOT_PASS" | cut -f 1 | sort | uniq); do
	rm $clusters_NOT_PASS;
done

#rinomino clusters rimasti
for clusters_PASS in clusterino*; do
	mv $clusters_PASS ${clusters_PASS#clusterino}.fa; 
done


mv *T2 ../
mv T1 ../
mv FILTER* ../
mv INFO.txt ../
################ CREO BED DEI CLUSTERS PASS ###################################

#creo files bed per ogni cluster "PASS" e li sposto nella cartella BedFiles
cp $BedFiles/*sort.bed .;
cat *bed > ALL.bed;
for cluster in $(ls *fa); do
        for seq in $(grep ">" $cluster | sed 's/>//g' | sed 's/;/\t/g' | cut -f 1); do 
		grep $seq ALL.bed >> $Clustering/$clusters_bed/${cluster%"fa"}"bed" ;
	done;
done

rm *bed


#########################CREO COVERAGE CLUSTERS PASS######################################

#creo file (*.genome) con nome e lunghezza della sequenza del mitocondrio da utilizzare nel passaggio successivo (creazione file *genomecov)
samtools faidx $genome_mit/$organism"_mit.fasta"
name=$(cut -f 1 $genome_mit/$organism"_mit.fasta.fai")
length=$(cut -f 2 $genome_mit/$organism"_mit.fasta.fai")
echo -e $name"\t"$length > $Clustering/$organism".genome"


cd $Clustering/$clusters_bed

#A partire da ogni bed creo files di copertura tot, al 3' e al 5' (genomecov) per ogni cluster
for bed in *bed
        do
        bedtools sort -i $bed > ${bed%.bed}.sorted.bed
        bedtools genomecov -d -i ${bed%.bed}.sorted.bed -g $Clustering/$organism.genome > ${bed%.bed}.genomecov
        bedtools genomecov -5 -d -i ${bed%.bed}.sorted.bed -g $Clustering/$organism.genome > ${bed%.bed}.genomecov.5
        bedtools genomecov -3 -d -i ${bed%.bed}.sorted.bed -g $Clustering/$organism.genome  > ${bed%.bed}.genomecov.3
done

rm $Clustering/*genome

############## CREO FASTA DEI MICRO PASS #############################################
cd $Clustering/$clusters_fasta

#creo fasta dei centroidi PASS
for i in $(ls *fa | sed 's/.fa//g'); do
grep -w -A 1 $(echo clusterid=$i) $Clustering/$output"_centroids.fa" > $Clustering/$i"_centroid.fa";
done

cd $Clustering

#per ogni sequenza definisco il centroide e lo cercaro nei file bed per definire la posizione. Aggiungo info posizione e strand nell'header
for i in $(ls *_centroid.fa | sed 's/_centroid.fa//g'); do
        centroid=$(grep ">" $i"_centroid.fa" | sed 's/>//g' | sed 's/;/\t/g' | cut -f 1)
        position=$(grep -w $centroid $Clustering/$clusters_bed/*sorted.bed | cut -f 2,3,6 | sed 's/\t/;/g')
        sed -i "1 s/$/;$position/g" $i"_centroid.fa"
done

#concateno i fasta micro in un multifasta
cat *_centroid.fa > $output"_centroids.PASS.fasta"

#formatto il multifasta
sed -i 's/;/_/g' $output"_centroids.PASS.fasta" 
sed -i  's/>.*clusterid/>clusterid/g' $output"_centroids.PASS.fasta"
sed -i 's/=//g' $output"_centroids.PASS.fasta"
sed -i 's/_/_pos/2' $output"_centroids.PASS.fasta"
sed -i 's/_/_strand/4' $output"_centroids.PASS.fasta"

#rm *_centroid.fa

####################### PLOTS R ##################################################
echo "" >> $home/smith.log
echo "START PLOTTING DATA AT $(date +%X)"  >> $home/smith.log
echo ""
cat $home/smith.log

cd $Plots


#copio e formatto i file di copertura in tabelle COV1 e COV2 (copertura delle reads sul MITO prima e dopo la rimozione delle reads nucleari)
cp $BedFiles/*Unique.bam.cov .
cp $BedFiles/*MITO1.bam.cov .

#creo file di Posizioni
for i in *MitoUnique.bam.cov; do 
	cut -f 2 $i > POS.temp; 
done

#per ogni campione creo file temporanei di copertura del primo mappaggio
for i in *MITO1.bam.cov; do
        cut -f 3 $i > $i.COV1.temp;
done

#unisco i file e sommo ogni riga (copertura totale per ogni posizione)
paste *COV1.temp | awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}'> SUM1.temp

#creo tabella COV1:copertura delle reads sul MITO PRIMA della rimozione delle reads nucleari
paste POS.temp *COV1.temp SUM1.temp > COV1.txt


#per ogni campione creo file temporanei di copertura del secondo mappaggio
for i in *MitoUnique.bam.cov; do
        cut -f 3 $i > $i.COV2.temp;
done

#unisco i file e sommo ogni riga (copertura totale per ogni posizione)
paste *COV2.temp | awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}'> SUM2.temp

#creo tabella COV2:copertura delle reads sul MITO DOPO la rimozione delle reads nucleari
paste POS.temp *COV2.temp SUM2.temp > COV2.txt

#copio file per i grafici 
cp $Clustering/T1 .
cp $Clustering/*T2 .
cp $Clustering/$output"_centroids.PASS.fasta" .

#eseguo script e sposto i risultati in cartella output.R
Rscript $scripts/esplora2.R COV1.txt COV2.txt $PWD $output"_centroids.PASS.fasta" 

#eseguo script PLOT Clusters
$scripts/Makeplots.py -d $Clustering/$clusters_bed

rm *temp
rm *cov
rm T1 *T2 *_centroids.PASS.fasta

cd $home

#disattivo ambiente conda
conda deactivate

exit 0
