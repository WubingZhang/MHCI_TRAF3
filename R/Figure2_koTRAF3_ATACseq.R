#### Cistrome-GO ####
## Quantification
library(ggplot2)
library(ggrepel)
cistromego1 = read.table("data/Figure2/Traf3_Veh_1k_CistromeGO_RP_result.txt",
                        sep = "\t", row.names = 1, quote = "")
colnames(cistromego1) = c("Coordinate", "PeakNum", "RP", "adjRP", "Rank")
cistromego1$Human = TransGeneID(rownames(cistromego1), "Symbol", "Symbol",
                                fromOrg = "mmu")
cistromego1 = na.omit(cistromego1)
cistromego1 = cistromego1[order(cistromego1$Rank), ]
gg = cistromego1[cistromego1$RP>0,]
gg$Rank = rank(-gg$RP)
gg$Label = rownames(gg)
toplabels = c("H2-K1", "H2-K2", "H2-Q4", "Nfkb2", "Stat4",
              "Il34", "Ccl9", "Nfkbie", "Ccl6", "Stum", "Lgals3")
p = ScatterView(gg, "Rank", "RP", label = "Label", main = "sgTraf3 vs sgRosa26")
p = p + scale_color_manual(values = "#377eb8", labels = "Vehicle")
p = p + geom_label_repel(color = "#377eb8", data = gg[rownames(gg)%in%toplabels, ], force = 2)
p = p + theme_bw(base_line_size = NA)
p = p + theme(legend.position = c(0.7,0.8), plot.title = element_text(hjust = 0.5))
p
ggsave("RankView_Traf3_veh_RP.pdf", p, width = 2.8, height = 3.2)


cistromego2 = read.table("data/Figure2/Traf3_IFNg_1k_CistromeGO_RP_result.txt",
                         sep = "\t", row.names = 1, quote = "")
colnames(cistromego2) = c("Coordinate", "PeakNum", "RP", "adjRP", "Rank")
cistromego2$Human = TransGeneID(rownames(cistromego2), "Symbol", "Symbol",
                                fromOrg = "mmu")
cistromego2 = na.omit(cistromego2)
cistromego2 = cistromego2[order(cistromego2$Rank), ]
gg = cistromego2[cistromego2$adjRP>0,]
gg$Rank = rank(-gg$adjRP)
gg$Label = rownames(gg)
toplabels = c("H2-K1", "H2-K2", "H2-Q4", "Nfkb2", "Stat4",
              "Il34", "Ccl9", "Nfkbie", "Ccl6", "Ccl5", "Stum", "Lgals3")
p = ScatterView(gg, "Rank", "RP", label = "Label", main = "sgTraf3 vs sgRosa26")
p = p + scale_color_manual(values = "#377eb8", labels = "IFNg")
p = p + geom_label_repel(color = "#377eb8", data = gg[rownames(gg)%in%toplabels, ], force = 2)
p = p + theme_bw(base_line_size = NA)
p = p + theme(legend.position = c(0.7,0.8), plot.title = element_text(hjust = 0.5))
p
ggsave("RankView_Traf3_IFNg_RP.pdf", p, width = 2.8, height = 3.2)

## Quantification
library(ggplot2)
bp1 = read.table("data/Figure2/Traf3_IFNg_1k_CistromeGO_go_kegg_result.txt",
                 sep = "\t", row.names = 1, quote = "")
bp2 = read.table("data/Figure2/Traf3_Veh_1k_CistromeGO_go_kegg_result.txt",
                 sep = "\t", row.names = 1, quote = "")
bp1$Description = gsub("\\(.*", "", rownames(bp1))
bp1$ID = gsub(".*\\(|\\).*", "", rownames(bp1))
colnames(bp1)[1:5] = c("statistics", "NES", "pvalue", "p.adjust", "Count")
bp1 = bp1[order(bp1$p.adjust), ]
bp2$Description = gsub("\\(.*", "", rownames(bp2))
bp2$ID = gsub(".*\\(|\\).*", "", rownames(bp2))
colnames(bp2)[1:5] = c("statistics", "NES", "pvalue", "p.adjust", "Count")
bp2 = bp2[order(bp2$p.adjust), ]
terms = c("mmu04612", "mmu04060", "mmu04145", "mmu04668", "mmu04064",
          "mmu04620", "mmu04062", "mmu04218")
p1 = EnrichedView(bp1, mode = 2, subset = terms)
p1 = p1 + scale_color_manual(values = "#3182bd", guide = FALSE)
p1 = p1 + labs(title = "Cistrome-GO (IFNg)")
ggsave("EnrichView_CistromeGO_TRAF3_IFNg.pdf", p1, width = 5.8, height = 2.5)
p2 = EnrichedView(bp2, mode = 2, subset = terms)
p2 = p2 + scale_color_manual(values = "#6baed6", guide = FALSE)
p2 = p2 + labs(title = "Cistrome-GO (Vehicle)")
ggsave("EnrichView_CistromeGO_TRAF3_Veh.pdf", p2, width = 5.8, height = 2.5)

#### Cistrome-toolkit ####
Toolkit = read.csv("data/Figure2/koTRAF3_IFNg_diff_peaks_top1k_Toolkit_result.csv",
                   header = TRUE, stringsAsFactors = FALSE)
Toolkit = Toolkit[Toolkit$Factor%in%unique(Toolkit$Factor)[1:8], ]
Toolkit$Factor = factor(Toolkit$Factor, levels = rev(unique(Toolkit$Factor)))
p <- ggplot(Toolkit, aes(x=Factor, y=GIGGLE_score, color=Factor))
p = p + geom_point()
p = p + scale_color_brewer(type = "qual", palette=3)
p = p + labs(x = NULL, y = "GIGGLE score", title = "sgTraf3 vs sgRosa26 (IFNg)")
p = p + theme_bw(base_line_size = NA)
p = p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
p = p + coord_flip()
p
ggsave("CistromeToolkit_koTRAF3_IFNg.pdf", p, width = 3.5, height = 2.3)

Toolkit = read.csv("data/Figure2/koTRAF3_Veh_diff_peaks_top1k_Toolkit_result.csv",
                   header = TRUE, stringsAsFactors = FALSE)
Toolkit = Toolkit[Toolkit$Factor%in%unique(Toolkit$Factor)[1:8], ]
Toolkit$Factor = factor(Toolkit$Factor, levels = rev(unique(Toolkit$Factor)))
p <- ggplot(Toolkit, aes(x=Factor, y=GIGGLE_score, color=Factor))
p = p + geom_point()
p = p + scale_color_brewer(type = "qual", palette=3)
p = p + labs(x = NULL, y = "GIGGLE score", title = "sgTraf3 vs sgRosa26 (Vehicle)")
p = p + theme_bw(base_line_size = NA)
p = p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
p = p + coord_flip()
p
ggsave("CistromeToolkit_koTRAF3_Veh.pdf", p, width = 3.5, height = 2.3)
