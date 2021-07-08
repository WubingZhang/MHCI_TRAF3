rm(list = ls())
require("GSVA")
require("ggpubr")
require("survival")
library("survminer")
library("Biobase")
date = format(Sys.time(), "%b.%d.%Y")
### Too large, so the data are not provided
exprdir = "TCGA/ProcessedDat/"
source('R/TCGA_Correlation.R')
source('R/Tigger.TCGA.R')
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

#### Association between koTRAF3 signature and TIL #####
res = Tigger.TCGA(genesig)

for(i in setdiff(names(res), c("summary","summary.p"))){
  p = res[[i]]$MHCI + scale_color_manual(values = "#e7298a") + ylab("MHC-I")
  ggsave(paste0("CorrView_koTRAF3_Signature_MHCI_", i, ".pdf"), p,
         width = 4, height = 3, dpi = 200, useDingbats=FALSE)
  p = res[[i]]$CTL + scale_color_manual(values = "#d95f02")
  ggsave(paste0("CorrView_koTRAF3_Signature_CTL_", i, ".pdf"), p,
         width = 4, height = 3, dpi = 200, useDingbats=FALSE)
  p = res[[i]]$PDL1 + scale_color_manual(values = "#1b9e77") + ylab("PD-L1")
  ggsave(paste0("CorrView_koTRAF3_Signature_MHCI-PDL1_", i, ".pdf"), p,
         width = 4, height = 3, dpi = 200, useDingbats=FALSE)
}
summary_res = res$summary
summary_res = summary_res[order(-summary_res$Signature.Zscore), ]
idx = summary_res$Cancer %in% c("ACC", "MESO", "THCA", "SARC", "SKCM" ,"COAD",
                                "READ", "BRCA.LumA", "KIRP", "BLCA", "LIHC",
                                "BRCA.Her2", "OV", "BRCA.Basal", "KICH", "KIRC")
summary_res = summary_res[idx, ]

## TRAF3 signature correlation with MHC-I and PD-L1
gg = summary_res[,c(3,5)]; colnames(gg) = gsub(".PCC","",colnames(gg))
p = DensityView(gg, xlab = "Correlation with Traf3 signature") +
  theme(legend.position = "top", legend.justification = c(0.5, 0.5)) +
  geom_vline(xintercept = robustbase::colMedians(as.matrix(gg)),
             linetype = "dashed", color = c("#fb8072","#80b1d3"))
tmp = t.test(gg$MHCI, gg$PDL1, paired = TRUE, alternative = "greater")
p = p + annotate(geom="text", x=0, y=2.4, label=paste0("pval: ",
                 format(tmp$p.value, digits = 3, scientific = TRUE)))
p
ggsave("DensityView_TCGA_Correlation_TRAF3sig_MHCI_PDL1.pdf", p,
       width = 3.5, height = 3.2)

summary_res$Cancer = factor(summary_res$Cancer,
                            levels = summary_res$Cancer[order(-summary_res$Signature.Zscore)])
gg = reshape2::melt(summary_res[, seq(1,ncol(summary_res),2)], id.vars = "Cancer")
pvalue = reshape2::melt(summary_res[, c(seq(2,ncol(summary_res),2),
                                        ncol(summary_res))], id.vars = "Cancer")
gg$Pvalue = pvalue$value
gg$Symbol = ""
gg$Symbol[gg$Pvalue<0.05] = "*"
gg$Symbol[gg$Pvalue<0.01] = "**"
gg$Symbol[gg$Pvalue<0.001] = "***"
gg = gg[!gg$variable%in%c("CD4.PCC", "PDL1.PCC", "CD8.PCC"), ]
gg$variable = as.character(gg$variable)
gg$variable[gg$variable=="Signature.Zscore"] = "TRAF3-KO.Zscore"
gg$variable = factor(gg$variable, levels = c("MHCI.PCC", "CTL.PCC", "TRAF3-KO.Zscore"))
p5 = ggplot(gg, aes(Cancer, value, fill=variable, label = Symbol))
p5 = p5 + geom_bar(stat="identity", width=0.7, alpha = 0.5)
p5 = p5 + geom_text(hjust = 1)
p5 = p5 + facet_grid(~variable, switch = "y", scales="free_x")
p5 = p5 + scale_fill_manual(values = c("#e7298a", "#d95f02", "#1b9e77", "#7570b3", "#e6ab02"))
p5 = p5 + coord_flip()
p5 = p5 + labs(x = NULL, y = NULL)
p5 = p5 + theme_pubr()
p5 = p5 + theme(legend.position="none")
p5
ggsave("BarView_Summary_TCGA_koTRAF3_Signature.pdf", p5, width = 7, height = 4)

#### SKCM: Traf3 signature score and TIDE prediction ####
SKCM_TIDE = read.csv("data/Figure5/SKCM_TIDE_prediction.csv", row.names = 1,
                     header = TRUE, stringsAsFactors = FALSE)
tmp = log2(readRDS("data/Figure5/SKCM.RNAseq.TPM.rds")+1)
genes = intersect(unlist(gset), rownames(tmp))
SampleAnn = data.frame(Sample=colnames(tmp), Patient = gsub("-..$", "", colnames(tmp)),
                       UG_DG = as.vector(t(tmp[genes, ])%*%weights[genes] / (length(genes)^0.5)),
                       Response1 = SKCM_TIDE[colnames(tmp), 1],
                       Response2 = SKCM_TIDE[colnames(tmp), 2], stringsAsFactors = FALSE)
SampleAnn$Response = NA
SampleAnn$Response[SampleAnn$Response1=="True"] = "Non-responder"
SampleAnn$Response[SampleAnn$Response2=="True"] = "Responder"
SampleAnn = na.omit(SampleAnn)
p = ggboxplot(SampleAnn, "Response", "UG_DG", notch=TRUE,
              color = "Response", width = 0.5,
              palette = c("#1f78b4", "#FC6665"))
p = p + stat_compare_means(comparisons = list(c("Responder", "Non-responder")))
p = p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
p = p + labs(x = NULL, color = NULL, y = "Traf3-knockout signature", title = "TCGA - SKCM")
p = p + ylim(NA, max(SampleAnn$UG_DG)+0.5)
p
ggsave(paste0("BoxView_SKCM.TIDE_koTRAF3_Signature.pdf"), p,
       width = 3, height = 3, dpi = 200, useDingbats=FALSE)

