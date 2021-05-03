## Public datasets and SMAC mimetics
rm(list=ls())
library(ComplexHeatmap); library(robustbase);
library(circlize);
library(quickAnalyze)
library(MAGeCKFlute)
library(ggplot2)
library(RColorBrewer)
date = format(Sys.time(), "%b.%d.%Y")
colorbar = c(Grey = "#BABABA", Salmon ="#f77d64", Blue = "#1f78b4")
options(stringsAsFactors = FALSE)

#### Environment preparation ####
#' To repeat this part, please process each drug treatment data downloaded from GEO and 
#' merge the LFCs from each dataset. The full annotation of drug treatment data is deposit in
#' data/Figure6/Drug_Perturbation_Annotation.txt.
mergedLFC = readRDS("path/to/Merged LFC.rds")
mergedLFC = limma::normalizeQuantiles(mergedLFC)
DrugAnn = read.table("data/Figure6/Drug_Perturbation_Annotation.txt", header = TRUE, sep = "\t", fill = TRUE, quote = "")
drug_type = DrugAnn$Category[!duplicated(DrugAnn$Treat)]
names(drug_type) = DrugAnn$Treat[!duplicated(DrugAnn$Treat)]
CellType = DrugAnn$CancerType[!duplicated(DrugAnn$CellType)]
names(CellType) = DrugAnn$CellType[!duplicated(DrugAnn$CellType)]
SampleAnn = t(as.data.frame(strsplit(colnames(mergedLFC), "_")))
rownames(SampleAnn) = colnames(mergedLFC)
colnames(SampleAnn) = c("GSE", "CellLine", "Treatment")
SampleAnn = as.data.frame(SampleAnn)
SampleAnn$CellLine = gsub("\\.T.*", "", SampleAnn$CellLine)
SampleAnn$Treatment = gsub("\\..*", "", SampleAnn$Treatment)
SampleAnn$CancerType = CellType[SampleAnn$CellLine]
SampleAnn$Category = drug_type[SampleAnn$Treatment]
cts = names(table(SampleAnn$CancerType))[table(SampleAnn$CancerType)<10]
SampleAnn$CancerType[SampleAnn$CancerType%in%cts] = "Others"
trts = names(table(SampleAnn$Category))[table(SampleAnn$Category)<10]
SampleAnn$Category[SampleAnn$Category%in%trts] = "Others"
SampleAnn$CancerType = factor(SampleAnn$CancerType)
SampleAnn$Category = factor(SampleAnn$Category,
                            levels = c("Targeted therapy","Chemotherapy",
                                       "Combination", "Others"))

#### Figure6A heatmap to show the log2FC of MHC-I and PD-L1 ####
gg = t(mergedLFC[c("3105","3106","3107","567", "6890", "29126", "80380"), ])
colnames(gg) = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "CD274", "PDCD1LG2")
gg = as.data.frame(gg)

gg = gg[order(rowMeans(gg[, c("HLA-A", "HLA-B", "HLA-C", "B2M")])), ]
bot_ann = columnAnnotation(foo = anno_lines(
  data.frame(MHCI = rowMeans(gg[, c("HLA-A", "HLA-B", "HLA-C", "B2M")]),
             PDL1 = rowMeans(gg[, c("CD274", "PDCD1LG2")])),
  gp = gpar(col = c(2,4)), pt_gp = gpar(col = c(2,4)),
  height = unit(1.5, "cm"), add_points = FALSE), show_annotation_name = FALSE)

top_ann = columnAnnotation(Cancer = SampleAnn[rownames(gg), "CancerType"],
                           Treatment = SampleAnn[rownames(gg), "Category"],
                           col = list(Cancer = getCols(unique(SampleAnn$CancerType),6),
                                      Treatment = getCols(unique(SampleAnn$Category),1)),
                           annotation_name_side = "right",
                           annotation_name_gp = gpar(fontsize = 11))

heatcolors = colorRamp2(seq(1.5,-1.5,length.out = 200),
                        colorRampPalette(brewer.pal(10, "RdYlBu"))(200))
p1 = Heatmap(t(gg), name = "logFC", color_space = "RGB",
             show_row_names = TRUE, show_column_names = FALSE,
             top_annotation = top_ann,
             bottom_annotation = bot_ann,
             cluster_rows = FALSE, cluster_columns = FALSE,
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 5),
             col = heatcolors,
             show_heatmap_legend = TRUE,
             heatmap_legend_param=list(
               labels_gp = gpar(fontsize = 10),
               title_gp = gpar(fontsize = 10),
               legend_width = unit(2, "mm"))
)
draw(p1, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "right")

pdf(paste0("Heatmap_DrugTreatment_MHCI_PDL1.pdf"),
    width = 7.5, height = 3)
draw(p1, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "right")
decorate_annotation("foo", {
  grid.lines(c(0, 1), c(0.42, 0.42), gp = gpar(lty = 1))
})
dev.off()

p2 = Heatmap(cor(gg), name = "Pearson correlation", color_space = "RGB",
             show_row_names = TRUE, show_column_names = FALSE,
             cluster_rows = TRUE, cluster_columns = FALSE,
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 10),
             col = heatcolors, show_heatmap_legend = TRUE,
             heatmap_legend_param=list(
               labels_gp = gpar(fontsize = 10),
               title_position = "leftcenter-rot",
               title_gp = gpar(fontsize = 10, fontface = "bold"),
               legend_width = unit(2, "mm"),
               legend_height = unit(40, "mm"))
)
draw(p2, padding = unit(c(2, 1, 1, 1), "mm"))
pdf("Heatmap_MHCI_PDL1_pearson.pdf", width = 3.5, height = 2)
draw(p2, padding = unit(c(2, 1, 1, 1), "mm"))
dev.off()


## Plot the TRAF3 signature scores of drugs
# TRAF3 signature
degres = readRDS("data/Figure5/DESeqRes_koTRAF3_IFNg.rds")
degres = degres[order(-abs(degres$stat)), ]
degres = degres[!duplicated(degres$Human), ]
degres = degres[order(degres$stat), ]
degres$Entrez = TransGeneID(degres$Human, "Symbol", "Entrez")
weights = degres$stat; names(weights) = degres$Entrez
weights = weights / max(abs(weights))
gset = list(UG = degres$Entrez[(nrow(degres)-199):nrow(degres)],
            DG = degres$Entrez[1:200])
genes = intersect(unlist(gset), rownames(mergedLFC))
TRAF3_sig = as.vector(t(mergedLFC[genes, ])%*%weights[genes] / (length(genes)^0.5))
gg = data.frame(Sample = colnames(mergedLFC), TRAF3_sig = TRAF3_sig,
                MHCI = colMeans(mergedLFC[c("3105","3106","3107","567", "6890",
                                            "6892", "29126", "80380"), ]),
                PDL1 = unlist(mergedLFC["29126",]),
                Label = gsub("\\..*", "", gsub(".*_", "", colnames(mergedLFC))))
p = ScatterView(gg, "TRAF3_sig", "MHCI")
p = p + geom_smooth(method = "lm")
p = p + theme(legend.position = "none")
p = p + labs(x = "TRAF3 knockout signature", y = "MHC-I expression change")
p = p + scale_color_manual(values = "gray30")
test = cor.test(gg$MHCI, gg$TRAF3_sig)
p = p + annotate(geom="text", label = paste0("r = ", round(test$estimate,3),
                                             "\np-value = ", round(test$p.value,3)),
                 parse = FALSE, x = max(gg$TRAF3_sig),
                 y = min(gg$MHCI), hjust=1, vjust=-1)
p
ggsave("Correlation_MHCI_Signature.pdf", p, width = 4, height = 3.2)

gg = data.frame(Sample = colnames(mergedLFC), TRAF3_sig = TRAF3_sig,
                MHCI = colMeans(mergedLFC[c("3105","3106","3107","567", "6890",
                                            "6892", "29126", "80380"), ]),
                PDL1 = unlist(mergedLFC["29126",]),
                Label = gsub("\\..*", "", gsub(".*_", "", colnames(mergedLFC))))

gg$Label[gg$Label%in%c("BV6", "SM164", "SM83")] = "SMAC mimetic"
gg$Label[gg$Label!="SMAC mimetic"] = "Others"
gg$CancerType = SampleAnn[rownames(gg), "CancerType"]
gg$CancerType[gg$Label!="SMAC mimetic"] = "Others"
p = ggplot(gg, aes(TRAF3_sig, MHCI, color=Label, size = Label))
p = p + geom_point(alpha = 0.7)
p = p + scale_size_manual(values = c(0.8,2.5), guide=FALSE)
p = p + scale_color_manual(values = c("#999999", "#ef8a62"), breaks = c("SMAC mimetic"))
p = p + theme_classic()
p = p + theme(legend.position = c(0.8,0.08),
              legend.margin = margin(-1, 2, 1, 1),
              legend.spacing.y = unit(0.5, "mm"),
              panel.border = element_rect(colour = "black", fill=NA),
              legend.background = element_blank(),
              legend.box.background = element_rect(colour = "black"))
p = p + labs(x = "TRAF3 knockout signature", y = "MHC-I expression change", color=NULL)
p
ggsave("Scatter_DrugTreatment_SMAC.pdf", p, width = 4, height = 3.2)

gg$Treat = gg$Label
gg$Label[gg$Label%in%c("BV6", "SM164")] = "SMAC mimetic"
gg$Label[gg$TRAF3_sig<0.5|gg$MHCI<0.5] = ""
gg$Size = "none"
gg$Size[gg$PDL1<0.5&gg$MHCI>0.5&gg$TRAF3_sig>0.5] = "selected"
gg$Label[gg$Size=="none"] = ""
library(ggpubr)
library(ggrepel)

p = ggplot(gg, aes(TRAF3_sig, MHCI, label = Label,
                   color = Color, size = Size))
p = p + geom_point(alpha = 0.6)
p = p + scale_color_manual(values = c("#8da0cb", "#e41a1c"))
p = p + scale_size_manual(values = c(1, 3))
p = p + theme_pubr()
p = p + theme(legend.position = "none")
p = p + geom_text_repel(size = 3.5, force = 5, color = "black")
p = p + labs(x = "TRAF3 knockout signature", y = "MHC-I expression change")
p
ggsave("Scatter_signature_SMAC.pdf", p, width = 6, height = 5.5)

