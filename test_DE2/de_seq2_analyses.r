##DE aanlyses

#import libraries
library(tximportData)
library(tximport)
library(rhd5f)
library(DESeq2)

#get working directory
dir <- getwd()

#read sample description file and location of 'abundance.h5'
samples <- read.table(file.path(dir,'samples.txt'), header=TRUE)
files <- file.path(dir, samples$run, 'abundance.h5')

#read kallisto files
names(files) <-  samples$run
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)

#import as DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(txi.kallisto$counts), colData = samples, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
summary(res)

#vedi il manuale a https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#tximport
#vedi comìè fatto il samples.txt a /home/PERSONALE/diego.carli2/kallisto_analyses/test_DE




#get all pairwise condition


pairs_lists <- list()

for (column in colnames(samples)[-1]) {
pairs_vector <- character(0)
my_list <- unique(samples[[column]])
for (i in 1:(length(my_list) - 1)) {
    for (j in (i + 1):length(my_list)) {
      pair_string <- paste(my_list[i], my_list[j], sep = "_")
      pairs_vector <- c(pairs_vector, pair_string)
pairs_lists[[column]] <- pairs_vector
    }
  }
}  

# Iterate through column names, excluding the first column
for (column in colnames(samples)[-1]) {
	dds <- DESeqDataSetFromMatrix(countData = round(txi.kallisto$counts), colData = samples, design = formula(paste("~", column)))
	dds <- DESeq(dds)
	res <- results(dds)
  my_list <- unique(samples[[column]])
  
 for (pair in pairs_lists[[column]]) {
	results_tmp <- results(dds, contrast=c(column, strsplit(pair,'_')[[1]][1], strsplit(pair,'_')[[1]][2]), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)
	assign(pair, results_tmp)
	print(summary(pair))
	#pair <- lfcShrink(dds, contrast=c(column, strsplit(pair,'_')[[1]][1], strsplit(pair,'_')[[1]][2]), res=pair)
		}
	}



