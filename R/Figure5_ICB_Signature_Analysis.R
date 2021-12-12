rm(list = ls())

date = format(Sys.time(), "%b.%d.%Y")
source('R/Tigger.ICB.R')
#### Read the signatures ####
degres = readRDS("data/Figure5/DESeqRes_koTRAF3_IFNg.rds")
degres = degres[order(-abs(degres$stat)), ]
degres = degres[!duplicated(degres$Human), ]
degres = degres[order(degres$stat), ]
weights = degres$stat; names(weights) = degres$Human
weights = weights / max(abs(weights))
gset = list(UG = degres$Human[(nrow(degres)-199):nrow(degres)],
			            DG = degres$Human[1:200])
genesig = weights[unlist(gset)]
# All ICB cohorts are collected from published studies, for more details, please visit http://tide.dfci.harvard.edu/download/
# For each dataset, we pre-processed expression data and clinical information and saved them into together into an ExpressionSet object. Then for each dataset, we have a rds file including the ExpressionSet object, which should be included in the "path/to/ICB" for automatic loading by the Tigger.ICB function.
res = Tigger.ICB(genesig, ICBdir = "path/to/ICB", outdir = "./")


