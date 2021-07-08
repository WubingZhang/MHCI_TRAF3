#' Differential expression analysis
#'
#' @docType methods
#' @name DEAnalyze
#' @rdname DEAnalyze
#'
#' @param obj Matrix like object or an ExpressionSet instance.
#' @param SampleAnn Matrix like object (only when obj is a matrix),
#' the rownames should match colnames in obj, and the first column should be Condition.
#' @param type "Array", "RNASeq" or "msms", only needed when obj is matrix like object.
#' @param method Differential expression analysis method, e.g. limma, DESeq2, GFOLD,
#' glm.pois, glm.qlll, and glm.nb.
#' @param paired Boolean, specifying whether perform paired comparison.
#' @param GeneAnn Matrix like object (only when obj is a matrix), the rownames should match rownames of obj.
#' @param return Character, either data.frame or ExpressionSet, specifying the return object type.
#' @param app.dir The path to application (e.g. GFOLD).
#'
#' @return An ExpressionSet instance.
#' @seealso \code{\link{ExpressionSet-class}}
#'
#' @author Wubing Zhang
#'
#' @importFrom Biobase exprs AnnotatedDataFrame pData sampleNames ExpressionSet
#' @importFrom limma lmFit eBayes topTable voom
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq lfcShrink
#' @importFrom msmsTests msms.glm.pois msms.glm.qlll msms.edgeR
#' @importFrom edgeR estimateDisp glmFit glmLRT topTags
#' @importFrom MSnbase MSnSet
#' @importFrom msmsEDA pp.msms.data
#' @importFrom stats na.omit model.matrix p.adjust
#' @export

DEAnalyze <- function(obj, SampleAnn = NULL, type = "Array",
                      method = "limma", paired = FALSE, GeneAnn = NULL,
                      return = c("data.frame", "ExpressionSet")[1],
                      app.dir = "/Users/Wubing/Applications/gfold/gfold"){
  #### Create a new object ####
  if(is.matrix(obj) | is.data.frame(obj)){
    obj = stats::na.omit(obj)
    colnames(SampleAnn)[1] = "Condition"
    if(paired) colnames(SampleAnn)[2] = "Sibs"
    samples = intersect(rownames(SampleAnn), colnames(obj))
    if(length(samples)<2) stop("Too small sample size !!!")
    expr <- as.matrix(obj[, samples])
    SampleAnn = SampleAnn[samples, , drop = FALSE]
    # obj = new("ExprDataSet", rawdata = expr, SampleAnn = SampleAnn, type = type)
    obj = Biobase::ExpressionSet(assayData = expr, phenoData = Biobase::AnnotatedDataFrame(SampleAnn))
    if(!is.null(GeneAnn)){
      tmp = GeneAnn[rownames(obj), drop = FALSE]
      rownames(tmp) = rownames(obj)
      slot(obj, "featureData") = Biobase::AnnotatedDataFrame(tmp)
    }
  }

  #### Build design matrix ####
  if(paired){
    Sibs = factor(Biobase::pData(obj)$Sibs)
    Condition = factor(Biobase::pData(obj)$Condition)
    design = stats::model.matrix(~Sibs+Condition)
  }else{
    design = stats::model.matrix(~1+Condition, Biobase::pData(obj))
    rownames(design) = Biobase::sampleNames(obj)
  }
  if(tolower(type) == "array"){
    # Biobase::exprs(obj) = normalizeQuantiles(Biobase::exprs(obj))
    #"ls" for least squares or "robust" for robust regression
    fit = limma::eBayes(limma::lmFit(Biobase::exprs(obj), design, na.rm=TRUE))
    res = limma::topTable(fit, adjust.method="BH", coef=ncol(design), number = Inf)
    res = res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
    colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
  }else if(tolower(type) == "rnaseq"){
    if(tolower(method) == "deseq2"){
      # Biobase::exprs(obj) = TransformCount(Biobase::exprs(obj), method = "vst")
      # DESeq2
      dds = DESeq2::DESeqDataSetFromMatrix(Biobase::exprs(obj), colData = Biobase::pData(obj), design = design)
      dds <- DESeq2::DESeq(dds)
      res0 = DESeq2::results(dds)
      res <- DESeq2::lfcShrink(dds, coef = ncol(design), quiet = TRUE)
      res$padj[is.na(res$padj)] = 1
      res = res[, c("log2FoldChange", "baseMean", "pvalue", "padj")]
      res$stat = res0$stat
      res = res[, c("log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(tolower(method) == "limma"){
      # Biobase::exprs(obj) = TransformCount(Biobase::exprs(obj), method = "voom")
      # limma:voom
      dge <- edgeR::DGEList(counts=Biobase::exprs(obj))
      dge <- edgeR::calcNormFactors(dge)
      dge <- limma::voom(dge, design, plot=FALSE)
      fit <- limma::eBayes(limma::lmFit(dge, design))
      res = limma::topTable(fit, adjust.method="BH", coef=ncol(design), number = nrow(Biobase::exprs(obj)))
      res = res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(tolower(method) == "edger"){
      # Biobase::exprs(obj) = TransformCount(Biobase::exprs(obj), method = "voom")
      dge <- edgeR::DGEList(counts=Biobase::exprs(obj))
      dge <- edgeR::calcNormFactors(dge)
      dge <- edgeR::estimateDisp(dge, design, robust=TRUE)
      fit <- edgeR::glmFit(dge, design)
      lrt <- edgeR::glmLRT(fit)
      res <- edgeR::topTags(lrt, n = nrow(Biobase::exprs(obj)))
      res = res$table[, c("logFC", "logCPM", "logFC", "PValue", "FDR")]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(tolower(method) == "gfold"){
      Biobase::exprs(obj) = TransformCount(Biobase::exprs(obj), method = "voom")
      # GFOLD
      tmp = mapply(function(x){
        write.table(cbind(NA, Biobase::exprs(obj)[,x], NA, NA),
                    file=paste0(colnames(Biobase::exprs(obj))[x], ".txt"),
                    sep="\t", col.names=FALSE)}, x=1:ncol(Biobase::exprs(obj)))
      lev = levels(Biobase::pData(obj)$Condition)
      ctrlname = rownames(Biobase::pData(obj))[Biobase::pData(obj)$Condition==lev[1]]
      treatname = rownames(Biobase::pData(obj))[Biobase::pData(obj)$Condition==lev[2]]
      system(paste0(app.dir, " diff -s1 ", paste0(ctrlname, collapse=","),
                    " -s2 ", paste0(treatname, collapse=","), " -suf .txt -o gfold_tmp")
      )
      res = read.table("gfold_tmp", row.names=1, stringsAsFactors = FALSE)
      res = res[, c(4, 5, 2, 3, 3)]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
      tmp = file.remove(paste0(ctrlname, ".txt"), paste0(treatname, ".txt"),
                        "gfold_tmp", "gfold_tmp.ext")
    }else{
      stop("Method not available for RNA-seq data !!!")
    }
  }else if(tolower(type) == "msms"){
    if (tolower(method) == "limma"){
      #"ls" for least squares or "robust" for robust regression
      fit = limma::eBayes(limma::lmFit(Biobase::exprs(obj), design))
      res = limma::topTable(fit, adjust.method="BH", coef=ncol(design), number = Inf)
      res = res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }else if(grepl("^glm\\.", tolower(method))){
      fd <- data.frame(gene = rownames(Biobase::exprs(obj)),
                       row.names = rownames(Biobase::exprs(obj)),
                       stringsAsFactors = FALSE)
      MSnSet_obj <- MSnbase::MSnSet(exprs=Biobase::exprs(obj), fData=fd,
                                    pData=Biobase::pData(obj))
      MSnSet_obj <- msmsEDA::pp.msms.data(MSnSet_obj)  # pp.msms.data function used to deleted genes which all expression is 0.

      null.f <- "y~1"
      alt.f <- "y~Condition"
      div <- colSums(Biobase::exprs(obj), na.rm = TRUE)
      ### msmsTests method
      if(tolower(method)=="glm.pois"){
        res <- msmsTests::msms.glm.pois(MSnSet_obj, alt.f, null.f, div=div)
      }else if(tolower(method)=="glm.qlll"){
        res <- msmsTests::msms.glm.qlll(MSnSet_obj, alt.f, null.f, div=div)
      }else if(tolower(method)=="glm.nb"){
        res <- msmsTests::msms.edgeR(MSnSet_obj, alt.f, null.f, div=div)
      }
      res$baseMean = rowMeans(Biobase::exprs(obj))[rownames(res)]
      res$padj = stats::p.adjust(res$p.value, method = "BH")
      res = res[, c(1,4,2:3,5)]
      colnames(res) = c("log2FC", "baseMean", "stat", "pvalue", "padj")
    }
  }else{
    stop("Data type error! ")
  }
  res = as.data.frame(res)
  if(return == "data.frame") return(res)
  tmp = cbind(obj@featureData@data, res[rownames(obj@featureData@data), ])
  rownames(tmp) = rownames(obj@featureData@data)
  slot(obj, "featureData") = Biobase::AnnotatedDataFrame(tmp)
  return(obj)
}
