#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
#if (length(args)!=2) {
#  stop("correct syntax is Rscript esplora.R /home/gm/smithRNA_pipeline/PROVE/PROVA4/C.capitata_BedFiles/R/COV1.txt /home/gm/smithRNA_pipeline/PROVE/PROVA4/C.capitata_BedFiles/R/COV2.txt")
#}
#fileprima <- args[1]
#filedopo <- args[2]

# input alternativo se non legge da commandline
fileprima <- 'COV1.txt'
filedopo <- 'COV2.txt'

# output
plot_nuc_mt <- gsub('.txt$' ,'.pdf', filedopo)
plot_replicates <- gsub('.txt$' ,'.replicates.pdf', filedopo)
thresholds <- gsub('.txt$' ,'.stats', filedopo)  
terminal_out <- c() # will print to screen after execution


# reads prima, coverage after remapping of all reads
# columns are: base, coverage in each replicate, total coverage
# reads dopo, coverage after remapping of non nuclear reads
prima <- read.csv(file=fileprima, header=F, sep='\t')
dopo <- read.csv(file=filedopo, header=F, sep='\t')
n_replicates <- ncol(prima)-2

# gives informative column names to both dataframes
names <- colnames(prima)
names[1] <- "base"
for (r in 1:n_replicates) {
    names[r+1] <- paste('R', r, sep='')
}
names[n_replicates+2] <- 'tot'
colnames(prima) <- names
colnames(dopo) <- names


# PLOT PICCHI NUC MT
pdf(plot_nuc_mt, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54)
par(mfrow=c(2,1), mar=c(2, 4.1, 1, 1))

# threshold globale
# soglia sopra la quale sono coperte il 5% delle basi nel totale del secondo remapping
tail5pc_dopo <- sort(dopo$tot)[round((length(dopo$tot)/100)*95)]
cat(paste('Global_threshold: ',tail5pc_dopo, '\n', sep=''),file=thresholds)
terminal_out <- c(terminal_out, tail5pc_dopo) # to screen after execution

# coverage prima e dopo, superimposti con threshold totale
# si vede quali e quanti picchi sono scomparsi col remapping nucleare
plot(prima$base, prima$tot, type='l', xlab=NA, ylab='coverage (all replicates), first and second remapping', col='red')
points(dopo$base, dopo$tot, type='l', col='darkgreen')
abline(h=tail5pc_dopo, col='blue')

# coverage dopo, con threshold totale, zoom sui picchi bassi
plot(dopo$base, dopo$tot, type='l', xlab=NA, ylab='coverage (all replicates, second remapping, detail)', col='darkgreen', ylim=c(0, 5*tail5pc_dopo))
abline(h=tail5pc_dopo, col='blue')

trash <- dev.off()

# print to Rstudio AND file, non funziona da terminale ma è milgiore in Rstudio
# cancellare i precedenti pdf() and dev.off()
# dev.print(pdf, plot_nuc_mt, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54) 

# PLOT REPLICHE
pdf(plot_replicates, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54)
par(mfrow=c(n_replicates,1), mar=c(2, 4.1, 1, 1))

# plot cycle
for (r in 1:n_replicates) {
    rname <- paste('R', r, sep='')
    tail5pc_replica <- sort(dopo[,r+1])[round(((length(dopo[,r+1]))/100)*95)]
    plot(dopo$base, dopo[,r+1], type='l', xlab=NA, ylab=paste('coverage (', rname, ')', sep=''), col='darkgreen', ylim=c(0, 5*tail5pc_replica))
    abline(h=tail5pc_replica, col='blue')
    cat(paste(rname, '_threshold: ',tail5pc_replica, '\n', sep=''),file=thresholds,append=TRUE)
    terminal_out <- c(terminal_out, tail5pc_replica) # to screen after execution
}
trash <- dev.off()

# print to Rstudio AND file, non funziona da terminale ma è migliore in Rstudio
# cancellare i precedenti pdf() and dev.off()
# dev.print(pdf, plot_replicates, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54) 




# INFO SU COVERAGE
# total
cat(paste('All_replicates_meancoverage_first_remapping: ', round(mean(prima$tot)), '\n', sep=''),file=thresholds,append=TRUE)
cat(paste('All_replicates_meancoverage_second_remapping: ', round(mean(dopo$tot)), '\n', sep=''),file=thresholds,append=TRUE)
cat(paste('Coverage_(mean)_decreased_to: ', round((mean(dopo$tot)/mean(prima$tot))*100), '%\n', sep=''),file=thresholds,append=TRUE)
# replicates
for (r in 1:n_replicates) {
    rname <- paste('R', r, sep='')
    cat(paste(rname, '_meancoverage_first_remapping: ', round(mean(prima[,r+1])), '\n', sep=''),file=thresholds,append=TRUE)
    cat(paste(rname, '_meancoverage_second_remapping: ', round(mean(dopo[,r+1])), '\n', sep=''),file=thresholds,append=TRUE)
}
# compare replicates
cat(paste('Relative_coverage_(means): '),file=thresholds,append=TRUE)
for (r in 1:n_replicates) {
    cat(paste(round((mean(dopo[,r+1])/mean(dopo$tot))*100), '% ', sep=''),file=thresholds,append=TRUE)
}
cat('\n',file=thresholds,append=TRUE)

# output thresholds to terminal
print(terminal_out)

