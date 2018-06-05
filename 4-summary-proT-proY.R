rm(list=ls())
load('./data/cds.perc.proY.RData')
load('./data/cds.perc.proT.RData')
load('./data/cds.locus.funcCat.RData')
funcCat.df <- data.frame(
  funCat = cds.locus.funcCat,
  proT = cds.perc.proT[names(cds.locus.funcCat)],
  proY = cds.perc.proY[names(cds.locus.funcCat)],
  stringsAsFactors = FALSE
)
g.mean.proT <- mean(na.omit(funcCat.df$proT))
g.mean.proY <- mean(na.omit(funcCat.df$proY))
prop.g <- round((c(g.mean.proT,g.mean.proY)/(sum(c(g.mean.proT,g.mean.proY))))*100,2)

funcCat.l <- plyr::ddply(funcCat.df,'funCat',.fun = function(x){
  m.pro <- c(proT = mean(na.omit(x$proT)),
             proY = mean(na.omit(x$proY)))
  per.pro <- round(m.pro/(sum(m.pro))*100,2)
  names(per.pro) <- paste('p',names(m.pro),sep='.')
  c(m.pro,per.pro)
})
##########
load('./data/func.df.pe_ppe.proT.proY.RData')
funcSubcat <- plyr::ddply(
  func.df.pe_ppe.proT.proY,
  'subgroup',
  .fun = function(x){
    m.pro <- c(proT = mean(na.omit(x$proT)),
               proY = mean(na.omit(x$proY)))
    per.pro <- round((m.pro/(sum(m.pro)))*100,2)
    names(per.pro) <- paste('p.',names(m.pro),sep='')
    c(m.pro,per.pro)
  }
)
source('./func-room.R')
library(shapes2)
pdf('./figures/3b-summary-map-genome-func-category.pdf',10,10)
plot(x = seq(from = -40,
             to = 40,length.out = 10),
     y = seq(from = -40,
             to = 40,length.out = 10),
     type = 'n',
     ann = FALSE,
     axes = FALSE
)

drawDonut(x = 0,
          y = 0,
          prop = rep(2,10),
          label = 'Genome',
          stroke.size = 0.65)
idx <- seq(1, 360, length.out = 11)[-1]

coord.xy <- shapes2::get.xy(x = 0,
                            y = 0,
                            r.x = 15,
                            r.y = 15,
                            start.angle = 0,
                            end.angle = 360
)

for(i in 1:length(idx)){
  drawDonut(x = coord.xy$x[idx[i]],
            y = coord.xy$y[idx[i]],
            prop = funcCat.l$p.proT[i],
            label = i,r = 0.5)
}


y.sub <- rep(-25,3)
x.sub <- seq(from = -10, 10,length.out = 3)
for(i in 1:nrow(funcSubcat)){
  drawDonut(x = x.sub[i],y = y.sub[i],
            prop = funcSubcat$p.proT[i],
            label = funcSubcat$subgroup[i],r = 0.5)
}
dev.off()