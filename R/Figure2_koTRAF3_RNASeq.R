rm(list = ls())
library(MAGeCKFlute)
library(rMAUPS)
library(ggpubr)
library(ComplexHeatmap)
source('R/VolcanoView.R')
date = format(Sys.time(), "%b.%d.%Y")

#### Prepare the koTRAF3 RNA-Seq data ####
rawcount = readRDS("data/Figure2/koTRAF3_rawcount.rds")
colnames(rawcount) = paste0(rep(c("sgCtrl_veh_", "sgTraf3_1_veh_", "sgTraf3_2_veh_",
                                  "sgCtrl_IFNg_", "sgTraf3_1_IFNg_", "sgTraf3_2_IFNg_"),
                                each = 2), c("r1", "r2"))
SampleAnn = data.frame(row.names = colnames(rawcount),
                       Condition = gsub("_1|_2|_3|_r.$", "", colnames(rawcount)),
                       stringsAsFactors = FALSE)
Mouse2Human = TransGeneID(rownames(rawcount), fromType = "Symbol", "Symbol",
                          fromOrg = "mmu", toOrg = "hsa")
Mouse2Human = Mouse2Human[!(is.na(Mouse2Human)|Mouse2Human=="")]
rawcount = rawcount[rownames(rawcount)%in%names(Mouse2Human), ]
normcount = TransformCount(rawcount, method = "vst")
rownames(normcount) = rownames(rawcount)

#### Heatmap shows the top differential genes ####
scaledata = t(scale(t(normcount[,1:6])))
gg = scaledata[rownames(degres1)[degres1$padj<0.05&degres1$log2FC>0], ]
topann = columnAnnotation("design"= anno_block(labels=c("sgCtrl_veh", "sgTraf3_veh_1", "sgTraf3_veh_2"),
                                               gp = gpar(fill = c("#bababa", "#fbb4ae", "#fbb4ae"))),
                          show_legend = FALSE, show_annotation_name = FALSE)
at = which(rownames(gg) %in% c("Nfkb2","Map3k14","Nfkbie","H2-K1","Tapbp","Relb","Nfkbia",
                         "Cxcl10","Ifngr2","H2-Q4","Mvp", "Ddx58", "H2-D1",
                         "Ccl9", "Irf1", "B2m", "Irf8", "Tnfrsf9", "Psmb9", "Psmb8",
                         "Psmb10", "Nlrc5", "Tap1", "Tap2"))
rightann = rowAnnotation(foo = anno_mark(at = at, labels = rownames(gg)[at]),
                         show_annotation_name = FALSE)
col1 = circlize::colorRamp2(c(min(gg), 0, max(gg)), c("#4A2FFF", "#F6F3F2", "#FF331E"))
col2 = circlize::colorRamp2(c(min(gg), 0, max(gg)), c("#67a9d8", "#ffffff", "#fc8d70"))
col3 = circlize::colorRamp2(c(min(gg), 0, max(gg)), c("#3c91d8", "#f7f7f7", "#e73739"))
p = Heatmap(gg[1:100,], cluster_rows = FALSE, cluster_columns = FALSE,
            col = col2, show_row_names = FALSE, row_names_side = "left",
            show_column_names = FALSE, column_names_side = "bottom",
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            column_names_rot = 90,
            right_annotation = rightann,
            top_annotation = topann,
            # row_split = tmpgenes$GOTERM,
            column_split = gsub("_r1$|_r2$","",rownames(SampleAnn)[1:6]),
            show_heatmap_legend = TRUE,
            heatmap_legend_param = list(
              title = "Expression",
              title_gp = gpar(fontsize = 10),
              labels_gp = gpar(fontsize = 10),
              legend_height = unit(0.5, "cm"),
              legend_width = unit(2, "in"),
              title_position = "lefttop",
              legend_direction = "horizontal"))
pdf(paste0("Heatmap_koTRAF3_RNASeq_Vehicle.pdf"), width = 5, height = 4.5)
draw(p, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "bottom")
dev.off()


scaledata = t(scale(t(normcount[,7:12])))
gg = scaledata[rownames(degres2)[degres2$padj<0.05&degres2$log2FC>0.5], ]
topann = columnAnnotation("design"= anno_block(labels=c("sgCtrl_IFNg", "sgTraf3_IFNg_1", "sgTraf3_IFNg_2"),
                                               gp = gpar(fill = c("#bababa", "#fbb4ae", "#fbb4ae"))),
                          show_legend = FALSE, show_annotation_name = FALSE)
at = which(rownames(gg) %in% c("Nfkb2","Map3k14","Nfkbie","H2-K1","Tapbp","Relb","Nfkbia",
                               "Cxcl10","Ifngr2","H2-Q4","Mvp", "Ddx58", "H2-D1",
                               "Ccl9", "Irf1", "B2m", "Irf8", "Tnfrsf9", "Psmb9", "Psmb8",
                               "Psmb10", "Nlrc5", "Tap1", "Tap2"))
rightann = rowAnnotation(foo = anno_mark(at = at, labels = rownames(gg)[at]),
                         show_annotation_name = FALSE)
col1 = circlize::colorRamp2(c(min(gg), 0, max(gg)), c("#4A2FFF", "#F6F3F2", "#FF331E"))
col2 = circlize::colorRamp2(c(min(gg), 0, max(gg)), c("#67a9d8", "#ffffff", "#fc8d70"))
col3 = circlize::colorRamp2(c(min(gg), 0, max(gg)), c("#3c91d8", "#f7f7f7", "#e73739"))
p = Heatmap(gg[1:100,], cluster_rows = FALSE, cluster_columns = FALSE,
            col = col2, show_row_names = FALSE, row_names_side = "left",
            show_column_names = FALSE, column_names_side = "bottom",
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            column_names_rot = 90,
            right_annotation = rightann,
            top_annotation = topann,
            # row_split = tmpgenes$GOTERM,
            column_split = gsub("_r1$|_r2$","",rownames(SampleAnn)[7:12]),
            show_heatmap_legend = TRUE,
            heatmap_legend_param = list(
              title = "Expression",
              title_gp = gpar(fontsize = 10),
              labels_gp = gpar(fontsize = 10),
              legend_height = unit(0.5, "cm"),
              legend_width = unit(1, "in"),
              title_position = "lefttop",
              legend_direction = "horizontal"))
pdf(paste0("Heatmap_koTRAF3_RNASeq_IFNg.pdf"), width = 5, height = 4.5)
draw(p, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "bottom")
dev.off()

#### Barplot shows the top enriched pathways ####
degres1 = readRDS("data/Figure2/DESeqRes_koTRAF3_Veh.rds")
degres1 = degres1[order(-degres1$stat), ]
genelist = degres1$stat; names(genelist) = rownames(degres1)
ortRes1 = EnrichAnalyzer(genelist[1:200], universe = rownames(degres1),
                         keytype = "Symbol", type = "GOBP", method = "HGT",
                         organism = "mmu")
degres2 = readRDS("data/Figure2/DESeqRes_koTRAF3_IFNg.rds")
degres2 = degres2[order(-degres1$stat), ]
genelist = degres2$stat; names(genelist) = rownames(degres2)
ortRes2 = EnrichAnalyzer(genelist[1:200], universe = rownames(degres2),
                         keytype = "Symbol", type = "GOBP", method = "HGT",
                         organism = "mmu")
df = ortRes1@result
df$logFDR = -log10(df$p.adjust)
df = df[df$ID%in%c("GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN",
                   "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
                   "GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY",
                   "GO_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE",
                   "GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
                   "GO_T_CELL_MEDIATED_CYTOTOXICITY",
                   "GO_RESPONSE_TO_TYPE_I_INTERFERON", "GO_NIK_NF_KAPPAB_SIGNALING",
                   "GO_RESPONSE_TO_INTERLEUKIN_1",
                   # "GO_MYD88_INDEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                   "GO_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE"), ]
df$Description = c("Antigen presentation of endogenous peptide",
                   "Interferon gamma signaling", "Adaptive immune response",
                   "Antigen presentation via MHC-I",
                   "T cell mediated cytotoxicity", "TNF signaling pathway",
                   "Response to type I interferon", "NIK NF kappab signaling",
                   "Response to interleukin 1", "Cytokine production")
p = BarView(df, "Description", 'logFDR', fill = "#8da0cb")
p = p + labs(x = NULL, y = expression(-log[10]*FDR),
             title = "sgTraf3 vs sgControl (Veh)") + coord_flip()
p
ggsave("BarView_enrichment_koTRAF3_Veh.pdf", p, width = 6.5, height = 3.5)

df = ortRes2@result
df$logFDR = -log10(df$p.adjust)
df = df[df$ID%in%c("GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN",
                   "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
                   "GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY", "GO_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE",
                   "GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY", "GO_T_CELL_MEDIATED_CYTOTOXICITY",
                   "GO_RESPONSE_TO_TYPE_I_INTERFERON", "GO_NIK_NF_KAPPAB_SIGNALING", "GO_RESPONSE_TO_INTERLEUKIN_1",
                   "GO_MYD88_INDEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY"), ]
df$Description = paste0(substr(df$Description,0,1), tolower(substr(df$Description,2,nchar(df$Description))))
df$Description = gsub("mhc class i", "MHC-I", df$Description)
df$Description = gsub("Nik nf", "NIK NF", df$Description)
df$Description = gsub("type i", "type I", df$Description)
df$Description = c("Antigen presentation via MHC-I",
                   "T cell mediated cytotoxicity", "TNF signaling pathway",
                   "NIK NF kappab signaling", "Antigen presentation of endogenous peptide",
                   "Interferon gamma signaling", "Response to type I interferon",
                   "Adaptive immune response", "Response to interleukin 1",
                   "Toll like receptor signaling")

p = BarView(df, "Description", 'logFDR', fill = "#8da0cb")
p = p + labs(x = NULL, y = expression(-log[10]*FDR),
             title = "sgTraf3 vs sgControl (IFNg)") + coord_flip()
p
ggsave("BarView_enrichment_koTRAF3_IFNg.pdf", p, width = 6.5, height = 3.5)
