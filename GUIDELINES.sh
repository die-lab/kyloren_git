#!/bin/bash

##HOW TO ASSEMBLY, SHORT GUIDE

echo 'do u compile manually the SRA entries?'
echo 'have the reads paired layout?'
echo '10 sec to discard, or let me work'
sleep 10

#download the raw reads, having SRA accession codes.
organism=MuMu
SRA_entries='SRR3520994 SRR3520995'
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $SRA_entries

#first thing first, asses the quality and other parameters of raw reads using fastq tool.
#the report can be be pulled together using multiqc.
#fastqc *.fastq

#find k-mer read error using Rcorrector (https://gigascience.biomedcentral.com/articles/10.1186/s13742-015-0089-y)
#it would delete many SNP as weel, but it does not matter for our purpose.
for fastq in *_1.fastq
do 
perl /home/PERSONALE/diego.carli2/app/rcorrector/run_rcorrector.pl -1 $fastq -2 ${fastq%_1.fastq}_2.fastq
#perl /home/PERSONALE/diego.carli2/app/rcorrector/run_rcorrector.pl -s $fastq
done
#not sure but at a first glance it looks like many quality paramters are worse then before. Need to know if this habit is common or not.

#now removing adapters, low quality base and so.
conda activate fastp
for fastq in *_1.cor.fq
do 
fastp --in1 $fastq --in2 ${fastq%_1.cor.fq}_2.cor.fq --out1 ${fastq%_1.cor.fq}_1.trimmed.fastq --out2 ${fastq%_1.cor.fq}_2.trimmed.fastq --detect_adapter_for_pe --average_qual 28
done

#filter out contaminations using kraken2 and kraken databases on the server
conda activate base
for trim in *1.trimmed.fastq
do 
kraken2 -db /home/dbs/kraken_dbs/ --paired --unclassified-out ${trim%.trimmed.fastq}.#.kraken.fastq $trim ${trim%1.trimmed.fastq}2.trimmed.fastq
#kraken2 -db /home/dbs/kraken_dbs/ --single --unclassified-out ${trim%.trimmed.fastq}.#.kraken.fastq $trim
done 

#maybe need to filter out rRNA using sortmerna, already installed in itsown conda environment (but a database is nedded).

#read normalization using khmer

#merge reads
cat *_1.kraken.fastq > left_kraken.fastq
cat *_2.kraken.fastq > right_kraken.fastq
#cat *kraken.fastq > single_kraken.fastq

#assembling using trinity, cat every fastq from where u want to start(after kraken probably) to a left file(1), an a right file(2). Or to a single file in case of single layout. 
Trinity --seqType fq --max_memory 10G --full_cleanup --CPU 12 --left left_kraken.fastq --right right_kraken.fastq 
#Trinity --seqType fq --max_memory 10G --full_cleanup --CPU 12 --single single_kraken.fastq 

transcriptome=$organism'_transcriptome'
mv trinity_out_dir.Trinity.fasta $transcriptome'.fasta'

#check the ortholog presence using busco command line
conda activate busco
busco -i $transcriptome'.fasta' --mode transcriptome -l metazoa -o assessment

#check the quality alignin the RNA-seq reads to the de novo assembled transcriptome. It should have above 80% of aligned reads for being considered of a good quality.
conda activate base
bowtie2-build -f $transcriptome'.fasta' $transcriptome
bowtie2 -x $transcriptome -1 left_kraken.fastq -2 right_kraken.fastq -S reads_against_transcriptome.sam --met-file reads_against_transcriptome.out
#bowtie2 -x $transcriptome -U single_kraken.fastq -S reads_against_transcriptome.sam --met-file reads_against_transcriptome.out

#extracting UTR
conda activate exutr
mkdir 3UTR
cd 3UTR
perl /home/PERSONALE/diego.carli2/app/ExUTR/3UTR_orf_20170816.pl -i ../$transcriptome'.fasta' -d /home/dbs/swissprot/swissprot -a 12 -o "$transcriptome"_3UTR -l un
perl /home/PERSONALE/diego.carli2/app/ExUTR/3UTR_ext_20170816.pl -i1 ../$transcriptome'.fasta' -i2 $transcriptome'_3UTR_orfs.fa' -a 12 -o "$transcriptome"_3UTR.fasta -x 2500
