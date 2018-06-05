rm(list=ls())
load('./data/cds.perc.proT.RData')
load('./data/cds.perc.proY.RData')
load('./data/cds.locus.funcCat.RData')
func.Cat <- unique(cds.locus.funcCat)
png('./figures/2a-genomewide-functional-category-proT-proY.png',
    width = 800,
    height = 1800,
    res = 150
)
par(mfrow=c(5,2), mar = c(3.1,4.1,4.1,2.1))
for(i in 1:length(func.Cat)){
  idx.cat <- which(cds.locus.funcCat %in% func.Cat[i])
  locus.cat <- names(cds.locus.funcCat[idx.cat])
  locus.perc.proT <- cds.perc.proT[locus.cat]
  locus.perc.proY <- cds.perc.proY[locus.cat]
  names(locus.perc.proT) <- NULL
  names(locus.perc.proY) <- NULL
  ###########
  idx.ord <- order(locus.perc.proT)
  ###########
  barplot(locus.perc.proY[idx.ord],
          border = '#ffcc00ff',
          col = '#ffcc00ff',
          ylab = '% proline codon usage',
          main = func.Cat[i],
          ylim = c(0,100))
  title(xlab = 'CDS index',line = 0)
  bar.proT <- barplot(
    locus.perc.proT[idx.ord],
    border = '#2a7fff10',
    col = '#2a7fff10',
    add = TRUE,
    ylim = c(0,100))
  title(xlab = 'CDS index',line = 0)
  points(x = bar.proT[,1],
         y = locus.perc.proT[idx.ord],
         type = 'l',
         lwd = 2.5,
         col = '#2a7fffff'
  )
  x.point <- round(length(bar.proT) * (70/100))
  legend(x=x.point,
         y = 95,
         legend = c('proY','proT'),
         fill = c('#ffcc00ff',
                  '#2a7fffff')
  )
}
dev.off()
#################
#################
### get PE/PPE category
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
getName <- function(x){
  strsplit(strsplit(x,split = 'Name=')[[1]][2],split = '\\;')[[1]][1]
}

func.df <- data.frame(
  locus = sapply(gff.file$V9, function(x)getLocus(x)),
  name = sapply(gff.file$V9, function(x)getName(x)),
  funcCat = sapply(gff.file$V9, function(x)getFuncCat(x)),
  stringsAsFactors = FALSE
)
func.df.pe_ppe <- func.df[func.df$funcCat %in% 'PE/PPE',]
subgroup <- sapply(func.df.pe_ppe$name,function(x){
  strsplit(x,split = '[0-9]')[[1]][1]
})
# https://mycobrowser.epfl.ch/quicksearch.php?gene+name=Rv1169c
subgroup[which(subgroup %in% 'lipX')] <- 'PE'
# https://mycobrowser.epfl.ch/quicksearch.php?gene+name=Rv3097c
subgroup[which(subgroup %in% 'lipY')] <- 'PE_PGRS'
#https://mycobrowser.epfl.ch/quicksearch.php?gene+name=Rv1759c
subgroup[which(subgroup %in% 'wag')] <- 'PE_PGRS'
func.df.pe_ppe <-  data.frame(
  func.df.pe_ppe,
  subgroup =subgroup,
  stringsAsFactors = FALSE
  )
idx.dup <- which(duplicated(func.df.pe_ppe$locus))
func.df.pe_ppe <- func.df.pe_ppe[-idx.dup,]
#####-------------
func.df.pe_ppe$proT <- cds.perc.proT[func.df.pe_ppe$locus]
func.df.pe_ppe$proY <- cds.perc.proY[func.df.pe_ppe$locus]
func.df.pe_ppe.proT.proY <- func.df.pe_ppe
save(func.df.pe_ppe.proT.proY,
     file = './data/func.df.pe_ppe.proT.proY.RData')

