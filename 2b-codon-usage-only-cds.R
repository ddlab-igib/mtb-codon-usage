rm(list=ls())
####
geneCodon.perc <- read.table(
  './data/proline-codon-percentage-per-gene-and-locus-Mtb-H37Rv-5June2018.csv',
  sep=',',
  header = TRUE,
  stringsAsFactors = FALSE)
##### only CDS
idx.cds <- which(geneCodon.perc$geneType %in% 'CDS')
geneCodon.cds.perc <- geneCodon.perc[idx.cds,]
#####
idx.ord <- order(geneCodon.cds.perc$CCC.proT.)
locus.ord <- geneCodon.cds.perc$geneLocus[idx.ord]
perc.proT <- geneCodon.cds.perc$CCC.proT.[idx.ord]
perc.proY <- geneCodon.cds.perc$CCG.proY.[idx.ord]
cds.perc.proT <- perc.proT
names(cds.perc.proT) <- geneCodon.cds.perc$geneLocus[idx.ord]
save(cds.perc.proT,
     file = './data/cds.perc.proT.RData')
cds.perc.proY <- perc.proY
names(cds.perc.proY) <- geneCodon.cds.perc$geneLocus[idx.ord]
save(cds.perc.proY,
     file = './data/cds.perc.proY.RData')
###########
# Figures
dir.create(
  'figures',
  showWarnings = FALSE
)
old.mar <- par()$mar
png('./figures/1a-genome-wide-dist-proT-proY-cds.png',
    width = 800,
    height = 700,
    res = 150
  )
par(mar = c(3.1,4.1,4.1,2.1))
barplot(perc.proY,
        border = '#ffcc00ff',
        col = '#ffcc00ff',
        ylab = '% proline codon usage',
        xlab = '',
        main = 'Genome-wide distribution of proT and\nproY codon usage in CDS')
title(xlab = 'CDS index',line=0, cex.lab=1.2)
bar.proT <- barplot(
  perc.proT,
  border = '#2a7fff10',
  col = '#2a7fff10',
  add = TRUE)
points(x = bar.proT[,1],
       y = perc.proT,
       type = 'l',
       lwd = 2.5,
       col = '#2a7fffff'
)
legend(x=3000,
       y = 90,
       legend = c('proY','proT'),
       fill = c('#ffcc00ff',
                '#2a7fffff')
)
dev.off()
par(mar = old.mar)
################
#####
genes.proT <- geneCodon.cds.perc[,2]
names(genes.proT) <- geneCodon.cds.perc$geneLocus
genes.proT <- na.omit(genes.proT)
genes.proY <- geneCodon.cds.perc[,3]
names(genes.proY) <- geneCodon.cds.perc$geneLocus
genes.proY <- na.omit(genes.proY)
png('./figures/1b-genomewide-proT-and-proY-codon-usage.png',
    width = 350,
    height = 650,
    res = 150
)
boxplot(list(proT=genes.proT,
             proY=genes.proY),
        ylab='Percentage Proline Codon Usage',
        col = c('#2a7fffff','#ffcc00ff')
)
dev.off()
###############
#####
