rm(list=ls())
load('./data/func.df.pe_ppe.proT.proY.RData')
library(plyr)
pe.ppe.pepgrs <- dlply(func.df.pe_ppe.proT.proY,'subgroup')
##########
########## PE
png('./figures/2b-PE-func-category-PE-subcategory.png',
    width = 500,
    height = 500,
    res = 100
    )
perc.proT <- pe.ppe.pepgrs[[1]]$proT
perc.proY <- pe.ppe.pepgrs[[1]]$proY

idx.ord <- order(perc.proT)
###########
barplot(perc.proY[idx.ord],
        border = '#ffcc00ff',
        col = '#ffcc00ff',
        ylab = '% proline codon usage',
        xlab = 'CDS index',
        main = names(pe.ppe.pepgrs)[1],
        ylim = c(0,100)
)
bar.proT <- barplot(
  perc.proT[idx.ord],
  border = '#2a7fff40',
  col = '#2a7fff40',
  add = TRUE,
  ylim = c(0,100))
points(x = bar.proT[,1],
       y = perc.proT[idx.ord],
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
dev.off()
##########
##########
png('./figures/2c-PE-func-category-PE_PGRS-subcategory.png',
    width = 500,
    height = 500,
    res=100)
perc.proT <- pe.ppe.pepgrs[[2]]$proT
perc.proY <- pe.ppe.pepgrs[[2]]$proY

idx.ord <- order(perc.proT)
###########
barplot(perc.proY[idx.ord],
        border = '#ffcc00ff',
        col = '#ffcc00ff',
        ylab = '% proline codon usage',
        xlab = 'CDS index',
        main = names(pe.ppe.pepgrs)[2],
        ylim = c(0,100)
)
bar.proT <- barplot(
  perc.proT[idx.ord],
  border = '#2a7fff40',
  col = '#2a7fff40',
  add = TRUE,
  ylim = c(0,100))
points(x = bar.proT[,1],
       y = perc.proT[idx.ord],
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
dev.off()
#############
#############
png('./figures/2d-PE-func-category-PPE-subcategory.png',
    width = 500,
    height = 500,
    res = 100)
perc.proT <- pe.ppe.pepgrs[[3]]$proT
perc.proY <- pe.ppe.pepgrs[[3]]$proY

idx.ord <- order(perc.proT)
###########
barplot(perc.proY[idx.ord],
        border = '#ffcc00ff',
        col = '#ffcc00ff',
        ylab = '% proline codon usage',
        xlab = 'CDS index',
        main = names(pe.ppe.pepgrs)[3],
        ylim = c(0,100)
)
bar.proT <- barplot(
  perc.proT[idx.ord],
  border = '#2a7fff40',
  col = '#2a7fff40',
  add = TRUE,
  ylim = c(0,100))
points(x = bar.proT[,1],
       y = perc.proT[idx.ord],
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
dev.off()
