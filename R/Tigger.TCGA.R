require(ggplot2)
require(quickAnalyze)
require(survival)
require(survminer)
#' Select subset features from TCGA RNA-seq data
#' @author Wubing Zhang
#'
TCGA_subset <- function(features = NULL,
                        cancers = NULL,
                        TIL.method = "CIBERSORT-ABS",
                        TCGAdir = "/TCGA/ProcessedDat/"){
  tumor_purity = readRDS(file.path(TCGAdir, "TCGA_TumorPurity.rds"))
  TCGA_clinical = readRDS(file.path(TCGAdir, "TCGA_Clinical.rds"))
  TPM_files = list.files(TCGAdir, pattern = "RNASeq.log2TPM.rds",
                         full.names = TRUE, recursive = TRUE)
  CancerType = gsub(".*\\/|_RNASeq.log2TPM.rds", "", TPM_files)
  if(!is.null(cancers)) CancerType = intersect(CancerType, cancers)
  summary = data.frame(stringsAsFactors = FALSE)
  for(CT in CancerType){
    path = grep(paste0("/", CT, "_RNASeq.log2TPM.rds$"), TPM_files, value = TRUE)
    tmp = readRDS(path)
    tmp = t(scale(t(tmp)))
    idx_norm = as.integer(gsub(".*-","",colnames(tmp)))>9
    if(sum(idx_norm,na.rm = TRUE)>1){
      tmp = tmp - rowMeans(tmp[,idx_norm],na.rm = TRUE)
    }else{
      tmp = t(scale(t(tmp), scale = FALSE))
    }
    if(class(features)=="list"){
      genescore = lapply(features, function(x){
        if(is.character(x)){
          colMeans(tmp[rownames(tmp)%in%x, , drop = FALSE], na.rm = TRUE)
        }else if(is.numeric(x)){
          ogene = intersect(rownames(tmp), names(x))
          as.vector(t(tmp[ogene,,drop=FALSE])%*%x[ogene])
        }
      })
    }else{
      genescore = t(tmp[rownames(tmp)%in%features, , drop = FALSE])
    }
    genescore = as.data.frame(genescore, stringsAsFactors = FALSE)
    genescore = cbind.data.frame(genescore, tumor_purity[rownames(genescore), -(1:2)])
    genescore$Cancer = CT
    genescore$Patient = gsub("...$", "", rownames(genescore))
    genescore$Sample = rownames(genescore)

    ## Immune cell infiltration
    TIL = readRDS(file.path(TCGAdir, "TCGA_TIL.rds"))
    TIL = TIL[, grepl(paste0(TIL.method,collapse = "|"),colnames(TIL), ignore.case = TRUE)]
    colnames(TIL) = gsub("_.*", "", colnames(TIL))
    TIL = as.data.frame(TIL, stringsAsFactors = FALSE)
    genescore = cbind.data.frame(genescore, TIL[genescore$Sample, ])

    ## Survival information
    clinical = TCGA_clinical[[gsub("\\W.*|_.*", "", CT)]]
    genescore$time = clinical[genescore$Patient, "OverallDays"]
    genescore$status = clinical[genescore$Patient, "Status"]
    genescore$race = NA
    if("race" %in% colnames(clinical))
      genescore$race = clinical[genescore$Patient, "race"]
    genescore$gender = NA
    if("gender" %in% colnames(clinical))
      genescore$gender = clinical[genescore$Patient, "gender"]
    genescore$pathologic_stage = NA
    if("pathologic_stage" %in% colnames(clinical))
      genescore$pathologic_stage = clinical[genescore$Patient, "pathologic_stage"]
    genescore$age = NA
    if("years_to_birth" %in% colnames(clinical))
      genescore$age = as.integer(clinical[genescore$Patient, "years_to_birth"])

    summary = rbind.data.frame(summary, genescore)
  }
  summary$stage = NA
  idx0 = grepl("0", summary$pathologic_stage)
  idx1 = grepl("i", summary$pathologic_stage)
  idx2 = grepl("ii", summary$pathologic_stage)
  idx3 = grepl("iii", summary$pathologic_stage)
  idx4 = grepl("iv", summary$pathologic_stage)
  summary$stage[idx0] = 0
  summary$stage[idx1] = 1
  summary$stage[idx2] = 2
  summary$stage[idx3] = 3
  summary$stage[idx4] = 4
  indicator = as.integer(gsub(".*-|.*\\.", "", summary$Sample))
  summary$Condition = NULL
  summary$Condition[indicator%in%(10:14)] = "Normal"
  summary$Condition[indicator%in%c(1,3,5,8,9)] = "Primary"
  summary$Condition[indicator%in%c(2,4)] = "Recurrent"
  summary$Condition[indicator%in%(6:7)] = "Metastatic"
  summary = summary[summary$Condition!="Normal", ]
  summary$time = summary$time/30
  colnames(summary) = gsub(" |\\+|\\(|\\)", "", colnames(summary))

  return(summary)
}



#' Tigger signature analysis
#' @author Wubing Zhang
#'
Tigger.TCGA <- function(dat,
                        sigName,
                        corFeatures = c(),
                        adjust = c("age", "gender", "race", "stage"),
                        interact = NULL,
                        partial = "IHC",
                        outdir = NULL,
                        xlab = "Signature",
                        time = "time",
                        status = "status"){
  SampleAnn = dat
  SampleAnn = SampleAnn[!(is.na(SampleAnn[,time])|is.na(SampleAnn[,status])), ]

  summary_res = list()
  res = list()
  for(f in unique(SampleAnn$Cancer)){
    message("Run TCGA: ", f)
    subSampleAnn = SampleAnn[SampleAnn$Cancer==f, ]
    if(nrow(subSampleAnn)<30) next
    # Remove variables with too many NAs
    subSampleAnn = subSampleAnn[, colSums(is.na(subSampleAnn))<nrow(subSampleAnn)/2]
    # Remove samples with NAs
    idx = !is.na(subSampleAnn[,sigName])
    subSampleAnn = subSampleAnn[idx, ]
    # Remove variables with unique value
    idx = apply(subSampleAnn, 2, function(x) sum(!is.na(unique(x))))
    subSampleAnn = subSampleAnn[, idx>1]
    ## Cox-regression
    tmpDat = subSampleAnn[, colnames(subSampleAnn)%in%c("time", "status", sigName, adjust)]
    if(length(interact)!=0 && interact%in%colnames(subSampleAnn)){
      tmpDat$Interaction = subSampleAnn[,sigName] * subSampleAnn[,interact]
      tmpDat$Partner = subSampleAnn[,interact]
    }
    cox <- coxph(Surv(time=time, event=status) ~ ., data=tmpDat)
    interestTerm = ifelse("Interaction"%in%colnames(tmpDat), "Interaction", sigName)
    summary_res[[f]] = c(summary(cox)$coefficients[interestTerm, 4:5])
    names(summary_res[[f]]) = c("Signature.Zscore", "Signature.pval")

    ## KM plot
    p1 = KMView(tmpDat, bio = sigName, os = "time", event = "status",
                adjust = adjust, interact = "Partner",
                labels = c("Signature.low", "Signature.high"),
                pval.pos = c(0, 0.2), optimalCut = TRUE, size = 0.6)
    p1 = suppressMessages(p1 + scale_color_manual(values = c("#1f78b4", "#FC6665")))
    p1 = p1 + theme(legend.position = "bottom")
    p1 = p1 + labs(title = f, x = "Time (month)")
    p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
    res[[f]][["km.p"]] = p1
    if(!is.null(outdir))
      ggsave(paste0(outdir, "/SurvView_", f, "_Signature.pdf"), p1,
             width = 4.5, height = 4, dpi = 200, useDingbats=FALSE)

    ## Association between Signature and immune features
    if(!is.null(partial) && partial%in%colnames(subSampleAnn)){
      subSampleAnn = subSampleAnn[!is.na(subSampleAnn[,partial]), ]
    }
    interstFeature = sigName
    if(length(interact)!=0 && interact%in%colnames(subSampleAnn)){
      subSampleAnn$Interaction = subSampleAnn[,sigName] * subSampleAnn[,interact]
      interstFeature = "Interaction"
    }
    for(i in corFeatures){
      if(!i%in%colnames(subSampleAnn)){
        summary_res[[f]] = c(summary_res[[f]], NA, NA)
        names(summary_res[[f]])[(length(summary_res[[f]])-1):length(summary_res[[f]])] =
          paste0(i, c(".PCC", ".pval"))
        next
      }
      # ## Correlation between features and signatures
      tmpDat = data.frame(x = subSampleAnn[,interstFeature], y = subSampleAnn[,i])
      if(!is.null(partial) && partial%in%colnames(subSampleAnn)){
        tmpDat$z = subSampleAnn[,partial]
        tmpDat = na.omit(tmpDat)
        pc = ppcor::pcor.test(tmpDat$x, tmpDat$y, tmpDat$z, )
      }else{
        tmpDat = na.omit(tmpDat)
        pc = cor.test(tmpDat$x, tmpDat$y)
      }
      summary_res[[f]] = c(summary_res[[f]], pc$estimate, pc$p.value)
      names(summary_res[[f]])[(length(summary_res[[f]])-1):length(summary_res[[f]])] =
        paste0(i, c(".PCC", ".pval"))
    }
  }
  summary_res = as.data.frame(t(as.data.frame(summary_res)))
  summary_res$Cancer = rownames(summary_res)
  res[["summary"]] = summary_res
  # summary_res = na.omit(summary_res)
  summary_res = summary_res[order(summary_res[,1]), ]
  summary_res$Cancer = factor(summary_res$Cancer, levels = summary_res$Cancer)

  gg = reshape2::melt(summary_res[, seq(1,ncol(summary_res),2)], id.vars = "Cancer")
  pvalue = reshape2::melt(summary_res[, c(seq(2,ncol(summary_res),2),
                                          ncol(summary_res))], id.vars = "Cancer")
  gg$Pvalue = pvalue$value
  gg$Symbol = ""
  gg$Symbol[gg$Pvalue<0.05] = "*"
  gg$Symbol[gg$Pvalue<0.01] = "**"
  gg$Symbol[gg$Pvalue<0.001] = "***"
  gg = na.omit(gg)
  p5 = ggplot(gg, aes(Cancer, value, fill=variable, label = Symbol))
  p5 = p5 + geom_bar(stat="identity", width=0.7, alpha = 0.5)
  p5 = p5 + geom_text()
  p5 = p5 + facet_grid(~variable, switch = "y", scales="free")
  # p5 = p5 + scale_fill_manual(values = c("#e7298a", "#7570b3", "#d95f02", "#1b9e77", "#e6ab02"))
  p5 = p5 + coord_flip()
  p5 = p5 + labs(x = NULL, y = NULL)
  p5 = p5 + theme_pubr()
  p5 = p5 + theme(legend.position="none")
  res[["summary.p"]] = p5
  if(!is.null(outdir))
    ggsave(paste0(outdir, "/BarView_summary_TCGA_Signature.pdf"), p5,
           width = 10, height = 8, dpi = 200, useDingbats=FALSE)
  return(res)
}
