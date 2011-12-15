#!/usr/bin/Rscript

source("../rfDistancesWithLikelihood.R")

## png("allTrees.png", width=1200, height=1200)
pdf("allTrees.pdf")
rfDistFile  = "RAxML_RF-Distances.tmp"
lnlFile = "lnl.tab"
lnlCol = 4
catCol = 2

rfDistancesWithLikelihood(rfDistFile, lnlFile, lnlCol, catCol) # , findML = TRUE
bla =dev.off()
