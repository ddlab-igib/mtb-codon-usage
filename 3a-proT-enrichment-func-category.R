rm(list=ls())
load('./data/cds.perc.proT.RData')
####
enriched.proT <- names(which(cds.perc.proT >= 60))
load('./data/cds.locus.funcCat.RData')
df.funcCat.proT60 <- data.frame(processes = cds.locus.funcCat[enriched.proT])
df.funcCat.proT60 <- cbind(rvID = rownames(df.funcCat.proT60),
                           df.funcCat.proT60)
write.csv(df.funcCat.proT60,
          file = './data/proT-60%-above-funcCategories.csv',
          row.names = FALSE)
png('./figures/3a-functional-annotation-of-proT-enriched-cds.png',
    width = 700,
    height = 700,
    res=100)
par(mar=c(5.1,16.1,4.1,2.1))
barplot(table(cds.locus.funcCat[enriched.proT]),las=2,horiz = TRUE,
        xlim = c(0,40),
        col = c(
          '#c83737ff',
          '#de8787ff',
          '#536c53ff',
          '#93ac93ff',
          '#37abc8ff',
          '#afdde9ff',
          '#d4aa00ff',
          '#e9ddafff',
          '#918a6fff',
          '#c8c4b7ff'),
        xlab = '# CDS',
        main = 'Functional annotation of proT\nenriched(>=60%) CDS'
)
dev.off()
#############
