#' Partial correlation between two gene sets in TCGA RNA-seq data
#' @author Wubing Zhang
#'
TCGA_Correlation <- function(gene1, gene2, summary = NULL,
                             method = "pearson", adjust = "Timer",
                             TCGAdir = "~/Jobs/Data_Archive/TCGA/ProcessedDat/",
                             scale = TRUE, center = TRUE){
  tumor_purity = readRDS(file.path(TCGAdir, "TCGA_TumorPurity.rds"))
  if(is.null(summary[1])){
    summary = TCGA_subset(genes = list(gene1 = gene1, gene2 = gene2),
                          TCGAdir = TCGAdir)
    gene1 = "gene1"; gene2 = "gene2"
  }

  source('~/Jobs/Project/ImmResponse/_Code/APC.Response/pcor.R')
  ## Timer tumor purity
  gg = summary[, c("Sample", "Cancer", gene1, gene2, adjust)]
  # if(!is.null(adjust)) gg = gg[!is.na(gg[, adjust]), ]
  CorrRes <- t(sapply(unique(gg$Cancer), function(x){
    tmp = gg[which(gg$Cancer==x), ]
    rownames(tmp) = tmp$Sample; tmp = tmp[,-(1:2)]
    # PCC = sapply(1:30, function(i){
      # idx = sample(1:nrow(tmp), 50, replace = TRUE)
      # tmp0 = tmp[idx, ]
    tmp0 = tmp
      if(!is.null(adjust) & sum(!is.na(tmp0[, adjust]))>3)
        pc = pcor.test(tmp0[!is.na(tmp0[, adjust]), gene1], tmp0[!is.na(tmp0[, adjust]), gene2],
                       tmp0[!is.na(tmp0[, adjust]), adjust], method = method)
      else
        pc = cor.test(tmp0[, gene1], tmp0[, gene2], method = method)
      c(pc$estimate, pc$p.value)
    # })
    # rowMeans(PCC)
  }))
  colnames(CorrRes) <- c("PC", "Pvalue")
  CorrRes = CorrRes[!is.na(CorrRes[,1]), ]

  return(CorrRes)
}




