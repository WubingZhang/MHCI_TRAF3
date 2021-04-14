## Set environment
rm(list = ls())
library(MAGeCKFlute);
library(data.table)
source("R/rScatterView.R")
source("R/RankView2.R")
############### Load CRISPR screen and prepare data for plotting ##################
date = format(Sys.time(), "%b.%d.%Y")
geneSum = readRDS("data/Figure1/MAGeCK_merged.GeneSummary.rds")
sgrnaSum = readRDS("data/Figure1/MAGeCK_merged.sgRNASummary.rds")
countSum = readRDS("data/Figure1/MAGeCK_normcounts.rds")

groups = unique(gsub("_1|_2$", "", names(geneSum)))
labels = c(expression("MHC-I"^"lo"*" PD-L1"^"lo"), expression("MHC-I"^"lo"),
           expression("PD-L1"^"lo"), expression("PD-L1"^"hi"),
           expression("MHC-I"^"hi"*" PD-L1"^"hi"), expression("MHC-I"^"hi"))
names(labels) = groups
labels = c("MHC-I & PD-L1\npositive regulator", "MHC-I specific\npositive regulator",
           "PD-L1 specific\npositive regulator", "PD-L1 specific\nnegative regulator",
           "MHC-I & PD-L1\nnegative regulator", "MHC-I specific\nnegative regulator")
names(labels) = groups
colors = c("#ff7f00", "#FC6665", "#1f78b4", "#1f78b4", "#ff7f00", "#FC6665")
names(colors) = groups
colors = c("#f768a1", "#FC6665", "#1f78b4", "#02818a", "#FF7F00", "#BF9000")
names(colors) = c("MHCI.low", "MHCI.high", "PDL1.low", "PDL1.high", "MHCI.low_PDL1.low", "MHCI.high_PDL1.high")
topgenes = list(c("Stat1", "Ifngr2", "Jak1", "Ifngr1", "Jak2"),
                c("B2m", "Tap1", "Tapbp", "Nlrc5", "H2-K1"),
                c("Cd274", "Jak2", "Jak1", "Ifngr2"),
                c("Urod", "Gtf2h3", "Krr1", "Dicer1"),
                c("Stub1", "Ube2n"),
                c("Traf3", "Tada3", "Med13", "Ezh2")) # "Traf3", "Tada3", "Cers3", "Emc1", "Syndig1"
names(topgenes) = groups
Mouse2Human = TransGeneID(unique(countSum$Gene), fromType = "Symbol", "Symbol",
                          fromOrg = "mmu", toOrg = "hsa")
Mouse2Human = Mouse2Human[!(is.na(Mouse2Human)|Mouse2Human=="")]

############### Visualize CRISPR screen results ##################
tophits = list()
merged = c()
for(gr in groups){
  tmp1 = ReadRRA(geneSum[[paste0(gr, "_1")]]);
  tmp2 = ReadRRA(geneSum[[paste0(gr, "_2")]]);
  tmp1 = tmp1[!is.na(tmp1$id), ]
  tmp2 = tmp2[!is.na(tmp2$id), ]
  rownames(tmp1) = tmp1$id
  rownames(tmp2) = tmp2$id
  tmp1 = tmp1[order(-tmp1$Score), ]
  tmp2 = tmp2[order(-tmp2$Score), ]
  tmp1$Score[tmp1$Score<0] = 0
  tmp2$Score[tmp2$Score<0] = 0
  gg = tmp1; gg$LFC = (tmp1$Score+tmp2[rownames(tmp1), "Score"])/2
  gg = gg[order(gg$LFC, decreasing = TRUE), ]
  gg$Human = Mouse2Human[gg$id]
  tmp = Mouse2Human[gg$id[gg$LFC>2]]
  tophits[[gr]] = tmp[!is.na(tmp)]

  color_pal = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                colors[gr], "gray70")[c(1:length(topgenes[[gr]]), 6:7)];
  names(color_pal) = c(topgenes[[gr]], "Z1", "Z2")

  gg = gg[1:500, ]
  p1 = RankView2(gg, text = "id", cutoff = 10, color_pal = color_pal,
                 top = 0, genelist = topgenes[[gr]],
                 legend.pos = c(0.7, 0.58), title = labels[gr])
  ggsave(paste0("Fig1_RankView_", gr, "_", date, ".pdf"), p1,
         width = 2.2, height = 2.8, dpi = 200, useDingbats=FALSE)
  genescore = gg$LFC; names(genescore) = gg$id

  color_pal = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                colors[gr], "gray70")[c(1:length(topgenes[[gr]]), 6:7)];
  names(color_pal) = c(topgenes[[gr]], "Z1", "Z2")
  p2 = rScatterView(genescore, top = 0, genelist = topgenes[[gr]], cutoff = 10,
                    color_pal = color_pal, title = labels[gr])
  p2 = p2 + labs(x = "RandomIndex", y = "Log2FC")
  ggsave(paste0("Fig1_rScatterView_", gr, "_", date, ".pdf"), p2,
         width = 2.8, height = 2.8, dpi = 200, useDingbats=FALSE)
}
