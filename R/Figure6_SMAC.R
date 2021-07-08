library(MAGeCKFlute)
require(ggplot2)
require(ggpubr)
library("ggfortify")
library(ggcorrplot)
options(stringsAsFactors = FALSE)
source("R/DEAnalyze.R")
rm(list = ls())
dir.create("Figure6", showWarnings = FALSE)

normdata = readRDS("data/Figure6/SMAC_ICB_RNASeq_normalized.rds")
scaledata = t(scale(t(normdata)))
genes = c("TAP1", "B2M", "CD8A", "CD8B", "PRF1", "GZMA", "GZMB")
genes = c("Tap1", "Tap2", "Tapbp", "B2m", "H2-K1", "H2-D1",
          "Traf3", "Rela", "Relb", "Nfkb2", "Nfkb1", "Rel",
          "Cd3d", "Cd3e", "Cd3g","Cd8a", "Cd8b1", "Gzma", "Gzmb", "Stat1",
          "Jak1", "Jak2", "Ifngr1", "Ifngr2")
HeatmapView(scaledata[genes, ], limit = c(-3, 3), cluster_cols = FALSE, display_number = TRUE)

#### Enrichment analysis ####
enrichRes1 = readRDS("data/Figure6/EnrichGSEA_SMAC_vs_Vehicle.rds")
enrichRes2 = readRDS("data/Figure6/EnrichGSEA_SMAC+ICB_vs_ICB.rds")
terms = c("GOBP:0045087", "GOBP:0006954", "GOBP:0006955", "GOBP:0050776",
          "GOBP:0042060", "GOBP:0019221", "GOBP:0002250", "GOBP:0033209",
          "GOBP:0002479", "GOBP:0002224", "GOBP:0070098", "GOBP:0038061")
gg = enrichRes1@result[rownames(enrichRes1@result)%in%terms, ]
gg["GOBP:0002479",2] = "antigen processing and presentation via MHC-I"
p = EnrichedView(gg, x = "NES", rank_by = "NES", top = 15, charLength = 100)
p = p + labs(title = "SMAC VS Vehicle")
p = p + scale_color_manual(values = "#66c2a5")
p
ggsave("Figure6/Fig6_EnrichedView_SMAC_vs_Vehicle.pdf", p, width = 5.5, height = 2.8)

gg = enrichRes2@result[rownames(enrichRes2@result)%in%terms, ]
gg["GOBP:0002479",2] = "antigen processing and presentation via MHC-I"
p = EnrichedView(gg, x = "NES", rank_by = "NES", top = 15, charLength = 100)
p = p + labs(title = "SMAC+ICB VS ICB")
p = p + scale_color_manual(values = "#fc8d62")
p
ggsave("Figure6/Fig6_EnrichedView_SMAC+ICB_vs_ICB.pdf", p, width = 6, height = 3)

# Immune cell infiltration
timer = read.csv("data/Figure6/TIMER_estimation.csv")
gg = reshape2::melt(timer, id.vars = "sampleID")
gg$Condition = gsub("_.*", "", gg$sampleID)
for(ct in unique(gg$variable)){
  p = ggplot(gg[gg$variable==ct, ], aes(Condition, y=value, color = Condition), alpha=0.5)
  p = p + geom_boxplot()
  p = p + geom_jitter()
  p = p + scale_color_manual(values = c("gray70", "#66c2a5", "#8da0cb", "#fc8d62"))
  p = p + labs(x = NULL, y = "TIMER estimation", fill = NULL, title = ct)
  p = p + theme_pubr()
  p = p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  p
  ggsave(paste0("Figure6/Fig6_Boxview_TIMER_", ct, ".pdf"), p, width = 4, height = 3.5)
}

# TRAF3 signature analysis
HumanGene = TransGeneID(rownames(normdata), "Symbol", "Symbol", fromOrg = "mmu")
idx = duplicated(HumanGene)|is.na(HumanGene)
normdata = normdata[!idx, ]
rownames(normdata) = HumanGene[!idx]

degres = readRDS("data/Figure4/DESeqRes_koTRAF3_IFNg.rds")
degres = degres[order(-abs(degres$stat)), ]
degres = degres[!duplicated(degres$Human), ]
degres = degres[order(degres$stat), ]
weights = degres$stat; names(weights) = degres$Human
weights = weights / max(abs(weights))
gset = list(UG = degres$Human[(nrow(degres)-199):nrow(degres)], DG = degres$Human[1:200])
genes = intersect(unlist(gset), rownames(normdata))
TRAF3_sig = as.vector(t(normdata[genes, ])%*%weights[genes] / (length(genes)^0.5))
gg = data.frame(Sample = colnames(normdata),
                Condition = gsub("_.*", "", colnames(normdata)),
                TRAF3_sig = TRAF3_sig)
p = ggplot(gg[gg$Condition%in%c("SMAC", "Vehicle"), ],
           aes(Condition, TRAF3_sig, color = Condition))
p = p + geom_boxplot()
p = p + geom_jitter()
p = p + scale_color_manual(values = c("gray70", "#66c2a5", "#8da0cb", "#fc8d62"))
p = p + theme_pubr()
p = p + labs(x = NULL, y = "TRAF3 knockout signature", color = NULL)
p = p + theme(legend.position = "none")
p
ggsave("Figure6/Fig6_Boxview_TRAF3_signature_SMAC.pdf", p, width = 2.5, height = 3)

p = ggplot(gg[gg$Condition%in%c("ICB", "SMAC.ICB"), ],
           aes(Condition, TRAF3_sig, color = Condition))
p = p + geom_boxplot()
p = p + geom_jitter()
p = p + scale_color_manual(values = c("#8da0cb", "#fc8d62"))
p = p + theme_pubr()
p = p + labs(x = NULL, y = "TRAF3 knockout signature", color = NULL)
p = p + theme(legend.position = "none")
p
ggsave("Figure6/Fig6_Boxview_TRAF3_signature_SMAC+ICB.pdf", p, width = 2.5, height = 3)
