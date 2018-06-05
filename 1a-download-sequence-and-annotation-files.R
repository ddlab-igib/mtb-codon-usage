rm(list=ls())
dir.create('./data',showWarnings = FALSE)
mtb.gene.seq.loc = 'https://mycobrowser.epfl.ch/releases/2/get_file?dir=fasta_genes&file=Mycobacterium_tuberculosis_H37Rv.fasta'
download.file(
  url = mtb.gene.seq.loc,
  destfile = './data/Mycobacterium_tuberculosis_H37Rv.fasta'
  
)
mtb.annotation.gff.loc = 'https://mycobrowser.epfl.ch/releases/2/get_file?dir=gff&file=Mycobacterium_tuberculosis_H37Rv.gff'
download.file(
  url = mtb.annotation.gff.loc,
  destfile = './data/Mycobacterium_tuberculosis_H37Rv.gff')
#----------
# CODON TABLE
codon.table <- list(
  'Ala' = c(
    'GCG',
    'GCC',
    'GCA',
    'GCT'
    ),
  'Gly' = c(
    'GGG',
    'GGC',
    'GGA',
    'GGT'
    ),
  'Pro' = c(
    'CCG',
    'CCC',
    'CCT',
    'CCA'
    ),
  'Thr' = c(
    'ACG',
    'ACC',
    'ACA',
    'ACT'
    ),
  'Val' = c(
    'GTG',
    'GTC',
    'GTA',
    'GTT'
  ),
  'Ser' = c(
    'TCG',
    'AGC',
    'TCC',
    'AGT',
    'TCA',
    'TCT'
  ),
  'Arg' = c(
    'CGG',
    'CGC',
    'CGA',
    'CGT',
    'AGG',
    'AGA'
  ),
  'Leu' = c(
    'CTC',
    'CTG',
    'CTT',
    'CTA',
    'TTG',
    'TTA'
  ),
  'Phe' = c(
    'TTC',
    'TTT'
  ),
  'Asn' = c(
    'AAC',
    'AAT'
  ),
  'Lys' = c(
    'AAG',
    'AAA'
  ),
  'Asp' = c(
    'GAC',
    'GAT'
  ),
  'Glu' = c(
    'GAG',
    'GAA'
  ),
  'His' = c(
    'CAC',
    'CAT'
  ),
  'Gln' = c(
    'CAG',
    'CAA'
  ),
  'Ile' = c(
    'ATC',
    'ATA',
    'ATT'
  ),
  'Met' = c(
    'ATG'
  ),
  'Tyr' = c(
    'TAC',
    'TAT'
  ),
  'Cys' = c(
    'TGC',
    'TGT'
  ),
  'trp' = c(
    'TGG'
  ),
  'Stop' = c(
    'TAG',
    'TGA',
    'TAA'
  )
)
save(
  codon.table,
  file = './data/codon.table.RData'
)
##################
############ ANNOTATION
####
#### Get the annotation
gff.file <- read.delim(
  './data/Mycobacterium_tuberculosis_H37Rv.gff',
  sep = '\t',
  header = FALSE,
  stringsAsFactors = FALSE)
getLocus <- function(x){
  strsplit(strsplit(x,split = 'Locus=')[[1]][2],split = '\\;')[[1]][1]
}
getFuncCat <- function(x){
  strsplit(x,split = 'Functional_Category=')[[1]][2]
}

gff.file.cds <- gff.file[which(gff.file$V3 %in% 'CDS'),]
cds.funcCategory <- as.character(sapply(gff.file.cds$V9, getFuncCat))
cds.locus <- as.character(sapply(gff.file.cds$V9,getLocus))

cds.locus.funcCat <- cds.funcCategory
names(cds.locus.funcCat) <- cds.locus
################
################
save(cds.locus.funcCat,
     file = './data/cds.locus.funcCat.RData')

names(cds.funcCategory) <- cds.locus
cds.funcCategory.locus <- cds.funcCategory
save(cds.funcCategory.locus,
     file = './data/cds.funcCategory.locus.RData')
