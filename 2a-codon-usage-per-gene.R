rm(list=ls())
source('./func-room.R')
codonPerGene <- read_fasta_and_compute_codon(
  fa_file = './data/Mycobacterium_tuberculosis_H37Rv.fasta'
  )
load('./data/codon.table.RData')
proline.codons <- codon.table['Pro'][[1]]
geneCodon <- plyr::dlply(
  .data = codonPerGene,
  .variables = 'fastaHeader',
  .fun = function(x){
    idx <- which(x$codons %in% proline.codons)
    x <- x[idx,]
    x$geneType <- sapply(as.character(x$fastaHeader),getType)
    x$geneLocus <- sapply(as.character(x$fastaHeader),getLocus)
    x$geneSymbol <- sapply(as.character(x$fastaHeader),getGene)
    x$percentage <- (x$Freq/sum(x$Freq)) * 100 
    return(x)
  }
)

geneCodon.per <- plyr::ldply(
  .data = geneCodon,
  .fun = function(x){
    df <- t(data.frame(x$percentage))
    colnames(df) <- x$codons
    return(df)
  }
)
geneCodon.per$geneSymbol <- sapply(
  as.character(geneCodon.per$fastaHeader),
  getGene
)
geneCodon.per$geneType <- sapply(
  as.character(geneCodon.per$fastaHeader),
  getType
)
geneCodon.per$geneLocus <- sapply(
  as.character(geneCodon.per$fastaHeader),
  getLocus
)
geneCodon.per<- geneCodon.per[,-1]
#--------------
geneCodon.freq <- plyr::ldply(
  .data = geneCodon,
  .fun = function(x){
    df <- t(data.frame(x$Freq))
    colnames(df) <- x$codons
    return(df)
  }
)
geneCodon.freq$geneLocus <- sapply(
  as.character(geneCodon.freq$fastaHeader),
  getLocus
)
geneCodon.freq$geneType <- sapply(
  as.character(geneCodon.freq$fastaHeader),
  getType
)
geneCodon.freq$geneSymbol <- sapply(
  as.character(geneCodon.freq$fastaHeader),
  getGene
)
geneCodon.freq<- geneCodon.freq[,-1]
#----------
proline.codon.name <- c(
  'CCT' = 'CCT(proX)',
  'CCC' = 'CCC(proT)',
  'CCA' = 'CCA(proU)',
  'CCG' = 'CCG(proY)'
)
colnames(geneCodon.per)[1:4] <- proline.codon.name[colnames(geneCodon.per)[1:4]]
colnames(geneCodon.freq)[1:4] <- proline.codon.name[colnames(geneCodon.freq)[1:4]]
idx.rbs <- which(geneCodon.per$geneType %in% 'RBS')
write.table(
  geneCodon.per[-idx.rbs,],
  './data/proline-codon-percentage-per-gene-and-locus-Mtb-H37Rv-5June2018.csv',
  sep = ',',
  row.names = FALSE
)
write.table(geneCodon.freq[-idx.rbs,],
            './data/proline-codon-frequency-per-gene-and-locus-Mtb-H37Rv-5June2018.csv',
            sep = ',',
            row.names = FALSE
)