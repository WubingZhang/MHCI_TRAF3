library(MAGeCKFlute)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)

rm(list = ls())
options(stringsAsFactors = FALSE)
dir.create("Figure4", showWarnings = FALSE)
normdata = readRDS("data/Figure4/MHCI_RNASeq_vst_normalized.rds")
normdata = as.matrix(normdata)
scaleexpr = t(scale(t(normdata)))
metadata = read.csv("data/Figure4/MHCI_RNASeq_metasheet.csv", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
metadata = metadata[order(metadata$MHC.I), ]
metadata$MHC.I = factor(tolower(metadata$MHC.I), levels = c("low", "high"))
colnames(metadata)[2] = "MHCI"

genes = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
          "TAPBP", "NLRC5", "PSMB8", "PSMB9")
gg = data.frame(Gene = factor(rep(genes,2), levels = genes),
                Group = rep(c("MHC-I.high", "MHC-I.low"), each=length(genes)))
gg$Expression = c(rowMeans(normdata[genes,1:5]), rowMeans(normdata[genes,6:10]))
gg$Sd = c(matrixStats::rowSds(normdata[genes,1:5]),
          matrixStats::rowSds(normdata[genes,6:10]))
tmp = apply(scaleexpr[genes,], 1, function(x) t.test(x[1:5], x[6:10])$p.value)
gg$Significance = tmp[gg$Gene]
gg$Symbol = "*"
gg$Symbol[gg$Significance<0.01] = "**"
gg$Symbol[gg$Significance<0.001] = "***"

p = ggplot(gg, aes(Gene, Expression, color = Group))
p = p + geom_bar(stat="identity", width = 0.7, position=position_dodge(0.8), fill = "white")
p = p + geom_errorbar(aes(ymin=Expression-Sd, ymax=Expression+Sd),
                      width=.2, position=position_dodge(0.8))
p = p + geom_text(aes(label = Symbol), hjust=-1, vjust = 0.5,
                  data = gg[!duplicated(gg$Gene),], color = "black")
p = p + scale_color_manual(values = c("#fb8072", "#1f7fb4"))
p = p + theme_pubr(legend = "bottom") + ylim(0, 17)
p = p + labs(x = NULL, y = "Normalized expression", color = NULL)
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + coord_flip()
p
ggsave("Figure4/Fig4_Barview_MHCI_expression.pdf", p, width = 3.8, height = 5)


gg = as.data.frame(t(scaleexpr[genes, ]))
gg$MHCI = paste0("MHC-I.", metadata[rownames(gg), 2])
gg = reshape2::melt(gg, id.vars = "MHCI")
tmp = apply(scaleexpr[genes,], 1, function(x) t.test(x[1:5], x[6:10])$p.value)
gg$Significance = tmp[gg$variable]
gg$Symbol = "*"
gg$Symbol[gg$Significance<0.01] = "**"
gg$Symbol[gg$Significance<0.001] = "***"
gg = gg[order(-gg$value), ]
p = ggplot(gg, aes(variable, value, fill = MHCI))
# p = p + geom_boxplot(position=position_dodge(0.8))
p = p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8))
p = p + stat_summary(aes(color = MHCI), fun.data=mean_sdl, fun.args = list(mult=1),
                     geom="pointrange", position=position_dodge(0.8))
p = p + geom_text(aes(y = min(gg$value), label = Symbol), hjust=0, vjust = 1,
                  data = gg[!duplicated(gg$variable),],
                  color = "black", fontface = "bold")
p = p + scale_color_manual(values = c("#fb8072", "#1f7fb4"))
p = p + scale_fill_manual(values = c("#fb8072", "#1f7fb4"))
p = p + theme_pubr(legend = "top")
p = p + labs(x = NULL, y = "Normalized expression", color = NULL, fill = NULL)
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + coord_flip()
p
ggsave("Figure4/Fig4_Dotview_MHCI_expression.pdf", p, width = 3, height = 4.5)

#### Read the signatures ####
degres = readRDS("data/Figure4/DESeqRes_koTRAF3_IFNg.rds")
degres = degres[order(-abs(degres$stat)), ]
degres = degres[!duplicated(degres$Human), ]
degres = degres[order(degres$stat), ]
weights = degres$stat; names(weights) = degres$Human
weights = weights / max(abs(weights))
gset = list(UG = degres$Human[(nrow(degres)-199):nrow(degres)], DG = degres$Human[1:200])
genes = intersect(unlist(gset), rownames(normdata))
metadata$TRAF3_sig = as.vector(t(normdata[genes, ])%*%weights[genes] / (length(genes)^0.5))
library(ggplot2)
library(ggpubr)
p = ggplot(metadata, aes(MHCI, TRAF3_sig, color = MHCI))
p = p + geom_boxplot(outlier.shape = NA)
p = p + geom_jitter()
p = p + stat_compare_means(method = "t.test")
p = p + scale_x_discrete(labels = c("MHC-I.low", "MHC-I.high"))
p = p + scale_color_manual(values = c("#1f7fb4", "#fb8072"), guide = FALSE)
p = p + theme_pubr()
p = p + labs(x = NULL, y = "Traf3-knockout signature")
p
ggsave("Figure4/Fig4_Boxplot_TRAF3_sig.pdf", p, width = 2.6, height = 3)

enrich2 = readRDS("data/Figure4/EnrichGSERes_Stevehodi_RNAseq.rds")
enrich2@result$pvalue = format(enrich2@result$p.adjust, scientific = TRUE, digits = 3)
enrich2@result$Description[enrich2@result$ID=="GOBP:0043123"] = "NF-kB signaling"
enrich2@result$Description[enrich2@result$ID=="GOBP:0050852"] = "T cell receptor signaling"
enrich2@result$Description[enrich2@result$ID=="GOBP:0019886"] = "Antigen presentation via MHC-I"
enrich2@result$Description[enrich2@result$ID=="GOBP:0033209"] = "TNF signaling"
enrich2@result$Description[enrich2@result$ID=="GOBP:0060333"] = "IFNg signaling"

cbp1 <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
p = ggView::gseaView(enrich2, c("GOBP:0043123", "GOBP:0050852", "GOBP:0019886",
                                "GOBP:0033209"),
             base_size = 12)
p = p + scale_color_manual(values = cbp1)
p = p + scale_y_continuous(limits = c(NA, 0.65), breaks = seq(0, 0.8, 0.2))
p = p + theme(legend.justification = c(0.5, 0.5), legend.position = c(0.65, 0.92))
p = p + theme(text = element_text(colour="black",size = 12, family = "Helvetica"),
              plot.title = element_text(hjust = 0.5, size=14),
              axis.text = element_text(colour="gray10"))
p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_blank(), panel.background = element_blank())
p = p + theme(legend.key.height = unit(4.5, "mm"))
p = p + labs(x = "Rank in ordered gene list")
p
ggsave("Figure4/Fig4_GSEAplot_MHCI_high.pdf", p, width = 5, height = 5)

idx = enrich2@result$ID %in% c("GOBP:0043123", "GOBP:0050852", "GOBP:0033209", "GOBP:0019886")
tmpgenes = strsplit(enrich2@result$geneName[idx], "\\/")
names(tmpgenes) = c("NF-kB \nsignaling", "TCR \nsignaling",
                    "TNF \nsignaling", "Antigen \npresentation")
cbp1 <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
names(cbp1) = names(tmpgenes)
tmpgenes = data.frame(Gene = unlist(tmpgenes), GOTERM = rep(names(tmpgenes), lengths(tmpgenes)))
topann = columnAnnotation("MHC-I"= anno_block(labels=c("MHC-I.low", "MHC-I.high"),
                                              gp = gpar(fill = c("#b3cde3", "#fbb4ae"))),
                          show_legend = FALSE, show_annotation_name = FALSE)
f1 = colorRamp2(seq(min(scaleexpr), max(scaleexpr), length = 3),
                c("blue", "#EEEEEE", "red"))
leftann = rowAnnotation(Pathway=tmpgenes$GOTERM, show_legend = FALSE,
                        col = list(Pathway = cbp1), show_annotation_name = FALSE)
p = Heatmap(scaleexpr[tmpgenes$Gene,],
            cluster_rows = FALSE, cluster_columns = FALSE,
            col = f1, show_row_names = FALSE, row_names_side = "left",
            show_column_names = FALSE, column_names_side = "bottom",
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            column_names_rot = 90,
            left_annotation = leftann,
            top_annotation = topann,
            row_split = tmpgenes$GOTERM,
            column_split = metadata$MHCI,
            show_heatmap_legend = TRUE,
            heatmap_legend_param = list(
              title = "Expression",
              title_gp = gpar(fontsize = 10),
              labels_gp = gpar(fontsize = 10),
              legend_height = unit(0.5, "cm"),
              legend_width = unit(2, "in"),
              title_position = "lefttop",
              legend_direction = "horizontal"))
pdf(paste0("Figure4/Fig4_Heatmap_Pathways_expression.pdf"), width = 3.5, height = 5)
draw(p, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "bottom")
dev.off()


############################################################################
#### SteveHodi Chip-seq datasets ####
############################################################################


#### Cistrome-GO ####
## Quantification
library(ggplot2)
library(ggrepel)
source('R/RankView2.R')
cistromego1 = read.table("data/Figure4/MHCI_chipseq_diffpeak_top1k_CistromeGO_RP_result.txt",
                         sep = "\t", row.names = 1, quote = "")
metadata = read.csv("data/Figure4/MHCI_chipseq_metasheet.csv", header = TRUE, row.names = 1)
colnames(cistromego1) = c("Coordinate", "PeakNum", "RP", "adjRP", "Rank")
cistromego1 = na.omit(cistromego1)
cistromego1 = cistromego1[order(cistromego1$Rank), ]
gg = cistromego1[cistromego1$adjRP>0,]
gg$Rank = rank(-gg$adjRP)
gg$Label = rownames(gg)
toplabels = c("B2M", "HLA-B", "HLA-C", "PSMB8", "PSMB9", "TAP1", "NLRC5", "HLA-A")
p = ScatterView(gg, "Rank", "RP", label = "Label", ylab = "Regulatory potential",
                main = "H3K27ac ChIP-Seq")
p = p + scale_color_manual(values = "#377eb8")
p = p + geom_label_repel(color = "#377eb8",
                         data = gg[rownames(gg)%in%toplabels, ], force = 3)
p = p + theme_bw(base_line_size = NA)
p = p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
p
ggsave("Figure4/Fig4_RankView_MHCI_chipseq_RP.pdf", p, width = 2.5, height = 3.2)

## Quantification
library(ggplot2)
bp1 = read.table("data/Figure4/MHCI_chipseq_diffpeak_top1k_CistromeGO_go_kegg_result.txt",
                 sep = "\t", row.names = 1, quote = "")
bp1$Description = gsub("\\(.*|pathway", "", rownames(bp1))
bp1$ID = gsub(".*\\(|\\).*", "", rownames(bp1))
colnames(bp1)[1:5] = c("statistics", "NES", "pvalue", "p.adjust", "Count")
bp1 = bp1[order(bp1$p.adjust), ]
terms = c("hsa04612", "hsa04621", "hsa04145", "hsa04658", "hsa04620")
p1 = EnrichedView(bp1, mode = 2, subset = terms)
p1 = p1 + scale_color_manual(values = "#3182bd", guide = FALSE)
p1 = p1 + labs(title = "Cistrome-GO (H3K27ac)")
p1 = p1 + xlim(1, NA)
p1
ggsave("Figure4/Fig4_EnrichView_CistromeGO_KEGG_H3K27ac.pdf", p1, width = 5, height = 2.5)

terms = c("GO:0019882", "GO:0001916", "GO:0050707", "GO:0032652", "GO:0051092")
bp2 = read.table("data/Figure4/MHCI_chipseq_diffpeak_top1k_CistromeGO_go_bp_result.txt",
                 sep = "\t", row.names = 1, quote = "")
bp2$Description = gsub("\\(.*|pathway", "", rownames(bp2))
bp2$ID = gsub(".*\\(|\\).*", "", rownames(bp2))
colnames(bp2)[1:5] = c("statistics", "NES", "pvalue", "p.adjust", "Count")
bp2 = bp2[order(bp2$p.adjust), ]
bp2 = bp2[bp2$ID%in%terms, ]
bp2$Description = gsub(".*regulation of ", "", bp2$Description)
p2 = EnrichedView(bp2, mode = 2, subset = terms)
p2 = p2 + scale_color_manual(values = "#3182bd", guide = FALSE)
p2 = p2 + labs(title = "Cistrome-GO (H3K27ac)")
p2
ggsave("Figure4/Fig4_EnrichView_CistromeGO_BP_H3K27ac.pdf", p2, width = 5.5, height = 2.5)

#### Cistrome-toolkit ####
Toolkit = read.csv("data/Figure4/MHCI_chipseq_diff_peaks_top1k_toolkit_result.csv",
                   header = TRUE, stringsAsFactors = FALSE)
Toolkit = Toolkit[Toolkit$Factor%in%unique(Toolkit$Factor)[1:11], ]
Toolkit = Toolkit[!Toolkit$Factor%in%c("PRDM1", "BRD4", "CIITA"), ]
Toolkit$Factor = factor(Toolkit$Factor, levels = rev(unique(Toolkit$Factor)))
p <- ggplot(Toolkit, aes(x=Factor, y=GIGGLE_score, color=Factor))
p = p + geom_point()
p = p + scale_color_brewer(type = "qual", palette=3)
p = p + labs(x = NULL, y = "GIGGLE score", title = "H3K27ac")
p = p + theme_bw(base_line_size = NA)
p = p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
p = p + coord_flip()
p
ggsave("Figure4/Fig4_CistromeToolkit_MHCI_chipseq.pdf", p, width = 3.5, height = 2.3)

RPmat = read.csv("data/Figure4/merged_RPScore.csv", header = TRUE, row.names = 1, check.names = FALSE)
RPmat = as.matrix(RPmat)
genes = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9")
idx1 = metadata[colnames(RPmat),1] == "MHCI_high"
idx2 = metadata[colnames(RPmat),1] == "MHCI_low"
gg = data.frame(Gene = factor(rep(genes,2), levels = genes),
                Group = rep(c("MHC-I.high", "MHC-I.low"), each=length(genes)))
gg$Expression = c(rowMeans(RPmat[genes,idx1]), rowMeans(RPmat[genes,idx2]))
gg$Sd = c(matrixStats::rowSds(RPmat[genes,idx1]),
          matrixStats::rowSds(RPmat[genes,idx2]))
tmp = apply(RPmat[genes,], 1, function(x) t.test(x[idx1], x[idx2])$p.value)
gg$Significance = tmp[gg$Gene]
gg$Symbol = ""
gg$Symbol[gg$Significance<0.1] = "*"
gg$Symbol[gg$Significance<0.01] = "**"
gg$Symbol[gg$Significance<0.001] = "***"

p = ggplot(gg, aes(Gene, Expression, color = Group))
p = p + geom_bar(stat="identity", width = 0.7, position=position_dodge(0.8), fill = "white")
p = p + geom_errorbar(aes(ymin=Expression-Sd, ymax=Expression+Sd),
                      width=.2, position=position_dodge(0.8))
p = p + geom_text(aes(y = 5.5, label = Symbol), hjust=-1, vjust = 0.5,
                  data = gg[!duplicated(gg$Gene),], color = "black")
p = p + scale_color_manual(values = c("#fb8072", "#1f7fb4"))
p = p + theme_pubr(legend = "bottom")
p = p + labs(x = NULL, y = "Normalized RP", color = NULL)
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + coord_flip()
p
ggsave("Figure4/Fig4_Barview_H3K27ac_MHCI_RP.pdf", p, width = 3.5, height = 5)


#### TRAF3 signature score ####
metadata = metadata[order(metadata$Type), ]
colnames(metadata)[1] = "MHCI"
metadata$MHCI = factor(metadata$MHCI, levels = c("MHCI_low", "MHCI_high"))
degres = readRDS("data/Figure4/DESeqRes_koTRAF3_IFNg.rds")
degres = degres[order(-abs(degres$stat)), ]
degres = degres[!duplicated(degres$Human), ]
degres = degres[order(degres$stat), ]
weights = degres$stat; names(weights) = degres$Human
weights = weights / max(abs(weights))
gset = list(UG = degres$Human[(nrow(degres)-199):nrow(degres)], DG = degres$Human[1:200])
RPmat[RPmat<3] = 0
genes = intersect(unlist(gset), rownames(RPmat))
metadata$TRAF3_sig = as.vector(t(RPmat[genes, rownames(metadata)])%*%weights[genes] / (length(genes)^0.5))

library(ggplot2)
library(ggpubr)
p = ggplot(metadata, aes(MHCI, TRAF3_sig, color = MHCI))
p = p + geom_boxplot(outlier.shape = NA)
p = p + geom_jitter()
p = p + scale_x_discrete(labels = c("MHC-I.low", "MHC-I.high"))
p = p + scale_color_manual(values = c("#1f7fb4", "#fb8072"), guide = FALSE)
p = p + theme_pubr()
p = p + labs(x = NULL, y = "Traf3-knockout signature", title = "Regulatory potential")
p = p + theme(plot.title = element_text(hjust = 0.5))
p
ggsave("Figure4/Fig4_Boxplot_RP_TRAF3_sig.pdf", p, width = 2.6, height = 3)

