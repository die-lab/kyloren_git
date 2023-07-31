#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("correct syntax is: Rscript esplora.R /path/COV1.txt /path/COV2.txt /path/folderconsoglie /path/fastacandidati.fas")
}
fileprima <- args[1]
filedopo <- args[2]
foldersoglie <- args[3]
fastacandidati <- args[4]


#input alternativo per testing
#fileprima <- 'COV1.txt'
#filedopo <- 'COV2.txt'
#foldersoglie <- './files_per_aggiornamento_script'
#fastacandidati <- './files_per_aggiornamento_script/nuoviheaders'


# legge soglie da files T1 e T2
T1 <- as.numeric(readLines(paste(foldersoglie, '/T1', sep = ""))) # soglia globale come numero
T2files <- list.files(path = foldersoglie, pattern = 'T2', full.names = TRUE) # lista files threshold repliche, namesorted
library(gtools)
T2files <- mixedsort(T2files) # numerical sort
T2 <- vector() # will take T2 of replicates
for (T2file in T2files) {
    T2 <- append(T2, as.numeric(readLines(T2file)))
} # T2 of replicates as one vector


# output file names
plot_nuc_mt <- gsub('.txt$' ,'.pdf', filedopo) # plot rossoverde
plot_replicates <- gsub('.txt$' ,'.replicates.pdf', filedopo) # plot repliche e soglie
plot_candidates <- gsub('.txt$' ,'.candidates.pdf', filedopo) # plot clusters
thresholds <- gsub('.txt$' ,'.stats', filedopo) # testo riassuntivo


# READS COVERAGE DATA
# reads prima, coverage after remapping of all reads
# reads dopo, coverage after remapping of non nuclear reads
# columns are: base, coverage in each replicate, total coverage
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


# PLOT PICCHI PRIMA E DOPO FILTRO NUCLEARE
pdf(plot_nuc_mt, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54)
par(mfrow=c(2,1), mar=c(2, 4.1, 1, 1))

# coverage prima e dopo, superimposti con threshold totale
# si vede quali e quanti picchi sono scomparsi col remapping nucleare
plot(prima$base, prima$tot, type='l', xlab=NA, ylab='coverage (all replicates), first and second remapping', col='red')
points(dopo$base, dopo$tot, type='l', col='darkgreen')
abline(h=T1, col='blue')

# coverage dopo, con threshold totale, zoom sui picchi bassi
# nota che la threshold è qui visualizzata sul coverage, ma in realtà si applica ai clusters.
plot(dopo$base, dopo$tot, type='l', xlab=NA, ylab='coverage (all replicates, second remapping, detail)', col='darkgreen', ylim=c(0, 5*T1))
abline(h=T1, col='blue')

trash <- dev.off()



# PLOT REPLICHE
# nota che la threshold è qui visualizzata sul coverage, ma in realtà si applica ai clusters.
pdf(plot_replicates, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54)
par(mfrow=c(n_replicates,1), mar=c(2, 4.1, 1, 1))

# plot cycle
for (r in 1:n_replicates) {
    rname <- paste('R', r, sep='')
    T2_replica <- T2[r] # soglia di questa replica, da T2
    plot(dopo$base, dopo[,r+1], type='l', xlab=NA, ylab=paste('coverage (', rname, ')', sep=''), col='darkgreen', ylim=c(0, 5*T2_replica))
    abline(h=T2_replica, col='blue')
}
trash <- dev.off()



# INFO SU THRESHOLDS E COVERAGE
# thresholds
cat(paste('Global_threshold: ',T1, '\n', sep=''),file=thresholds)
cat(paste('Replicates_thresholds: ', paste(T2, collapse=" "), '\n', sep=''),file=thresholds,append=TRUE) 
# total coverage
cat(paste('All_replicates_meancoverage_first_remapping: ', round(mean(prima$tot)), '\n', sep=''),file=thresholds,append=TRUE)
cat(paste('All_replicates_meancoverage_second_remapping: ', round(mean(dopo$tot)), '\n', sep=''),file=thresholds,append=TRUE)
cat(paste('Coverage_decreased_to: ', round((mean(dopo$tot)/mean(prima$tot))*100), '%\n', sep=''),file=thresholds,append=TRUE)
# replicate coverage
for (r in 1:n_replicates) {
    rname <- paste('R', r, sep='')
    cat(paste(rname, '_meancoverage_first_remapping: ', round(mean(prima[,r+1])), '\n', sep=''),file=thresholds,append=TRUE)
    cat(paste(rname, '_meancoverage_second_remapping: ', round(mean(dopo[,r+1])), '\n', sep=''),file=thresholds,append=TRUE)
}
# compare replicates
cat(paste('Relative_coverage: '),file=thresholds,append=TRUE)
for (r in 1:n_replicates) {
    cat(paste(round((mean(dopo[,r+1])/mean(dopo$tot))*100), '% ', sep=''),file=thresholds,append=TRUE)
}
cat('\n',file=thresholds,append=TRUE)


#PLOT CANDIDATES
# read candidates from fasta
fasta <- readLines(fastacandidati)
headers <- fasta[grepl("^>", fasta)]

# headers to dataframe
candidates = data.frame(from=numeric(0),to=numeric(0),strand=character(0),depth=numeric(0),clusterid=numeric(0))
for (c in 1:length(headers)) { # for each candidate
        candidate <- unlist(strsplit(gsub('>clusterid','',headers[c]), "_size|_pos|_|_strand")) # from, to, strand, depth, clusterid
        candidater <- list(as.numeric(candidate[3]), as.numeric(candidate[4]), candidate[5] , as.numeric(candidate[2]), as.numeric(candidate[1])) # ordine, tipo
        candidates <- rbind(candidates, candidater) # aggiunge i candidati al dataframe
}   
colnames(candidates) <- c('from', 'to', 'strand', 'depth', 'clusterid') 
candidates <-candidates[order(candidates$from),] # ordina per start

# cerca spazio verticale giusto fra scritte in base alla scala y
vs <- ((max(candidates[,4])*1.1)/30)

# aggiunge posizione media ed altezza per posizionare le etichette come colonne 6 e 7 
candidates$center <- (candidates[,1]+candidates[,2])/2
candidates$ylabel <- candidates[,4]+vs # proprio sopra il cluster, se va bene

# sposta in alto le labels che sovrascrivono, 5 successive approssimazioni
for (rounds in 1:5){ # aggiustamenti successivi
    for (r in 2:nrow(candidates)){
        for (rp in 1:(r-1)){ # ho messo r-1, era r
            # se una delle precedenti è vicina sia su x che su y
            if (((candidates[r,6] - candidates[rp,6]) < 400) & (abs(candidates[r,7] - candidates[rp,7]) < vs)){
                candidates[r, 7] = candidates[r, 7]+vs 
            } # alza quelli vicini
        }
    } 
}

# ylim del plot
mp <- (max(candidates[,7])*1.2) # etichetta massima +20%

# divide il genoma in 4 campi per aumentare la leggibilità
genomelen <- dopo[nrow(dopo),1]
primo <- c(-50, floor(genomelen/4)+50)
secondo <- c(ceiling(genomelen/4)-50, floor(genomelen/2)+50)
terzo <- c(ceiling(genomelen/2)-50, floor(genomelen*3/4)+50)
quarto <- c(ceiling(genomelen*3/4)-50, genomelen+50)
campi <- list(primo, secondo, terzo, quarto)

# PLOT
pdf(plot_candidates, onefile=T, paper='A4', width = 21/2.54, height = 29.7/2.54)

par(mfrow=c(4,1), mar=c(2, 4.1, 0.5, 1), xaxs="i", yaxs="i") # xaxs='i' rende i margini precisi
for (campo in campi) { # for each of the four plotting fields
    # plot coverage and global threshold
    plot(dopo$base, dopo$tot, type='l', xlim=campo, ylim=c(0, mp), xlab=NA, ylab='depth (detail)', col='grey')
        abline(h=T1, col='blue')
    # plot rectangles
    for (c in 1:nrow(candidates)) { # for each candidate
        if (candidates[c,3] == '+'){ # se + plotta rosso
            rect(candidates[c,1], 0, candidates[c,2], candidates[c,4], border='red')
        }
        if (candidates[c,3] == '-'){ # se - plotta verde
            rect(candidates[c,1], 0, candidates[c,2], candidates[c,4], border='green')
        }
    }
    text(candidates[,6], candidates[,7], candidates[,5], cex=0.7) # plotta etichette
}

trash <- dev.off()    
    
  
