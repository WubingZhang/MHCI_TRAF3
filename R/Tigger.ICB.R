#' Tigger signature analysis on ICB cohorts
#' @param signature A gene signature, either weighted or unweighted
#' @param datasets ICB cohorts for analysis
#' @param features A list, including features for correlation analysis with the signature.
#' Default features include MHCI, PDL1, IFNg, CTL, CD4.
#' @param adjust Confounding factor (one of the features) to adjust in the correlation analysis.
#' @param outdir Output directory, if provided, output all the figures and tables.
#' @param xlab The label of x-axis in the figures.
#' @param ICBdir The path to the ICB datasets.
#'
#' @return A list.
#'
#' @import plotROC ggplot2 survminer survival Biobase
#' @author Wubing Zhang
#'
Tigger.ICB <- function(signature,
					                          datasets = NULL,
											                         features = NULL,
											                         adjust = NULL,
																	                        outdir = NULL,
																	                        xlab = "Traf3 knockout signature",
																							                       ICBdir = "~/Jobs/Data_Archive/ClinicalTrial/RDS/"){
	  requireNamespace("ggplot2")
  requireNamespace("survminer")
    requireNamespace("survival")
    requireNamespace("Biobase")
	  requireNamespace("plotROC")
	  if(is.character(signature)){
		      tmp = rep(1, length(signature))
	      names(tmp) = signature
		      signature = tmp
		    }
	    defined = list(MHCI = c("HLA-A", "HLA-B", "HLA-C", "B2M"),
					                    PDL1 = c("CD274"),
										                 IFNg = c("STAT1", "IFNGR1", "IFNGR2", "JAK1", "JAK2"),
										                 CTL = c("CD8A", "CD8B", "PRF1", "GZMA", "GZMB"),
														                  CD4 = "CD4")
	    if(!is.null(features)) defined = c(defined, features)

		  exprsets = list.files(ICBdir, "ExprSet.rds", full.names = TRUE)
		  metadata = strsplit(gsub(".*\\/|_ExprSet.rds", "", exprsets), "_")
		    metadata = as.data.frame(metadata)
		    metadata = as.data.frame(t(metadata))
			  colnames(metadata) = c("Dataset", "ICB", "Cancer")
			  rownames(metadata) = paste0(metadata$Dataset, "_", metadata$ICB)
			    names(exprsets) = rownames(metadata)
			    if(!is.null(datasets))
					    metadata = metadata[rownames(metadata)%in%datasets, ]
				  summary = data.frame(stringsAsFactors = FALSE)
				    summary_res = list()
				    res = list()
					  for(ds in rownames(metadata)){
						      message("Run ICB dataset: ", ds)
					    tmp = readRDS(exprsets[ds])
						    expr = as.data.frame(exprs(tmp))
						    expr = t(scale(t(expr), scale = FALSE))
							    clinical = pData(tmp)
							    if(class(defined)=="list"){
									      genescore = lapply(defined, function(x)
															         colMeans(expr[rownames(expr)%in%x, , drop = FALSE], na.rm = TRUE))
								    }else{
										      genescore = t(expr[defined, , drop = FALSE])
									    }
								    genescore = as.data.frame(genescore, stringsAsFactors = FALSE)
									    # genescore$MHCI_PDL1 = genescore$MHCI - genescore$PDL1
									    genescore$Dataset = ds
									    idx = colnames(clinical)%in%c("OS","OS.Event","PFS","PFS.Event","RECIST","binaryResponse")
										    genescore = cbind.data.frame(genescore, clinical[rownames(genescore), idx])

										    if(nrow(genescore)<20) next
											    genes = intersect(names(signature), rownames(expr))
											    if(length(genes)<1) stop("Incorrect gene names ...")
												    tiggerScore = as.vector(t(expr[genes,,drop=FALSE])%*%signature[genes] / (length(genes)^0.5))
												    names(tiggerScore) = colnames(expr)
													    genescore$tiggerScore = tiggerScore[rownames(genescore)]
													    ## Correlation between survival and signature
													    ## Overall survival analysis
													    summary_res[[ds]] = c(NA,NA,NA,NA)
														    names(summary_res[[ds]]) = c("Coxph.Z", "Coxph.pval", "Response.stat", "Response.pval")
														    if(all(c("OS","OS.Event")%in%colnames(genescore)) &&
															          min(table(genescore$OS.Event))>2){
																      p1 = KMView(genescore, bio = "tiggerScore", os = "OS", event = "OS.Event",
																				                    labels = c("Signature.low", "Signature.high"),
																									                  optimalCut = TRUE, size = 0.6)
															      p1 = p1 + scale_color_manual(values = c("#1f78b4", "#FC6665"))
																        p1 = p1 + theme(legend.position = "bottom")
																        p1 = p1 + labs(title = ds, x = "Time (month)")
																		      p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
																		      res[[ds]][["OS.p"]] = p1
																			        if(!is.null(outdir))
																						        ggsave(paste0(outdir, "/SurvView_", ds, ".OS_Signature.pdf"), p1,
																									                  width = 4.5, height = 4, dpi = 200, useDingbats=FALSE)

																			        tmp <- surv_cutpoint(genescore, time = "OS", event = "OS.Event",
																										                            variables = "tiggerScore")
																					      tmpAnn<- surv_categorize(tmp, labels = c("low", "high"))
																					      tmpAnn$tiggerScore = factor(tmpAnn$tiggerScore, levels = c("low", "high"))
																						        cox <- coxph(Surv(time=OS, event=OS.Event) ~ tiggerScore, data=tmpAnn)
																						        if(!is.null(adjust)&&adjust%in%colnames(genescore)){
																									        cox <- coxph(Surv(time=OS, event=OS.Event) ~ tiggerScore+genescore[,adjust],
																														                      data=genescore)
																								      }
																								      summary_res[[ds]][1:2] = summary(cox)$coefficients[1, 4:5]
																								    }
															    if(all(c("PFS","PFS.Event")%in%colnames(genescore)) &&
																          min(table(genescore$PFS.Event))>2){
																	      p1 = KMView(genescore, bio = "tiggerScore", os = "PFS", event = "PFS.Event",
																					                    labels = c("Signature.low", "Signature.high"),
																										                  optimalCut = TRUE, size = 0.6)
																      p1 = p1 + scale_color_manual(values = c("#1f78b4", "#FC6665"))
																	        p1 = p1 + theme(legend.position = "bottom")
																	        p1 = p1 + labs(title = ds, x = "Time (month)", y = "Progression-Free Survival")
																			      p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
																			      res[[ds]][["PFS.p"]] = p1
																				        if(!is.null(outdir))
																							        ggsave(paste0(outdir, "/SurvView_", ds, ".PFS_Signature.pdf"), p1,
																										                  width = 4.5, height = 4, dpi = 200, useDingbats=FALSE)

																				        tmp <- surv_cutpoint(genescore, time = "PFS", event = "PFS.Event",
																											                            variables = "tiggerScore")
																						      tmpAnn<- surv_categorize(tmp, labels = c("low", "high"))
																						      tmpAnn$tiggerScore = factor(tmpAnn$tiggerScore, levels = c("low", "high"))
																							        cox <- coxph(Surv(time=PFS, event=PFS.Event) ~ tiggerScore, data=tmpAnn)
																							        if(!is.null(adjust)&&adjust%in%colnames(genescore)){
																										        cox <- coxph(Surv(time=PFS, event=PFS.Event) ~ tiggerScore+genescore[,adjust],
																															                      data=genescore)
																									      }
																									      summary_res[[ds]][1:2] = summary(cox)$coefficients[1, 4:5]
																									    }
																    if("binaryResponse" %in% colnames(genescore) &&
																	          min(table(genescore$binaryResponse))>2){
																		      genescore$binaryResponse = factor(genescore$binaryResponse, levels = c("SD/PD", "CR/PR"))
																	      p = ggboxplot(genescore[!is.na(genescore$binaryResponse), ], "binaryResponse",
																						                    "tiggerScore", notch=FALSE, color = "binaryResponse", width = 0.5,
																											                    palette = c("#1f78b4", "#FC6665"))
																		        p = p + stat_compare_means(comparisons = list(c("SD/PD", "CR/PR")),
																										                                    method = "wilcox.test", label = "p.signif",
																																			                                 method.args = list(alternative = "less"))
																		        # p = p + scale_x_discrete(labels = c("CR     PR", "SD     PD"))
																		        p = p + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
																				      p = p + labs(x = NULL, color = NULL, y = xlab, title = ds)
																				      p = p + ylim(NA, max(genescore$tiggerScore)+0.6)
																					        p
																					        res[[ds]][["biResponse.p"]] = p
																							      if(!is.null(outdir))
																									          ggsave(paste0(outdir, "/BoxView_", ds, ".biResponse_Signature.pdf"), p,
																													                width = 2.8, height = 3.5, dpi = 200, useDingbats=FALSE)
																							      tiggerscore_resp = genescore$tiggerScore[as.character(genescore$binaryResponse)=="CR/PR"]
																								        tiggerscore_nonresp = genescore$tiggerScore[as.character(genescore$binaryResponse)=="SD/PD"]
																								        test = wilcox.test(tiggerscore_resp, tiggerscore_nonresp, alternative = "greater")
																										      summary_res[[ds]][3] = median(tiggerscore_resp,na.rm = TRUE) - median(tiggerscore_nonresp,na.rm = TRUE)
																										      summary_res[[ds]][4] = test$p.value

																											        ## ROC curve
																											        gg = genescore
																											        gg$Response = 0
																													      gg$Response[gg$binaryResponse=="CR/PR"] = 1
																													      p <- ggplot(gg, aes_string(d = "Response", m = "tiggerScore")) +
																															          plotROC::geom_roc(color = "#fc8d59", n.cuts = 0) +
																																	          geom_abline(slope = 1, linetype = "dashed", color = "#91bfdb")
																																		        p <- p + theme_bw()
																														        p <- p + labs(x = "False positive rate", y = "True positive rate", title = ds)
																																      p <- p + annotate("text", x = 0.75, y = 0.25,
																																						                        label = paste0("AUC = ", round(calc_auc(p)$AUC,3)),
																																												                        color = "#f46d43")
																																      p
																																	        res[[ds]][["roc.p"]] = p
																																	        if(!is.null(outdir))
																																				        ggsave(paste0(outdir, "/ROC_", ds, ".biResponse_Signature.pdf"), p,
																																							                  width = 2.8, height = 3.5, dpi = 200, useDingbats=FALSE)
																																			    }
																	    res[[ds]][["summary"]] = genescore
																		    features = intersect(names(defined), colnames(genescore))
																		    # features = c(features, "MHCI_PDL1")
																		    for(i in features){
																				      # ## Correlation between features and signatures
																				      pc = cor.test(genescore$tiggerScore, genescore[,i])
																			      # if(!is.null(adjust)&&adjust%in%colnames(subSampleAnn)){
																			      #   pc = ppcor::pcor.test(genescore$tiggerScore, genescore[,i], genescore[,adjust])
																			      # }
																			      summary_res[[ds]] = c(summary_res[[ds]], pc$estimate, pc$p.value)
																				        names(summary_res[[ds]])[(length(summary_res[[ds]])-1):length(summary_res[[ds]])] =
																							        paste0(i, c(".PCC", ".pval"))
																								      p = ggplot(genescore, aes_string("tiggerScore", i))
																				        p = p + geom_point(aes(color = "#e7298a"), alpha = 0.3)
																						      p = p + geom_smooth(method='lm')
																						      p = p + theme_pubr(legend = "none")
																							        p = p + labs(title = ds, x = xlab)
																							        p = p + theme(plot.title = element_text(hjust = 0.5))
																									      label = paste("PCC: ", round(pc$estimate,3),
																														                    "\npval: ", format(pc$p.value, digits=2))
																									      p = p + annotate("text", x = min(genescore$tiggerScore),
																														                          y = max(genescore[,i]),
																																				                         label = label, hjust=0, vjust=1)
																										        res[[ds]][[i]] = p
																										        if(!is.null(outdir))
																													        ggsave(paste0(outdir, "/CorrView_", ds, "_Signature_", i, ".pdf"), p,
																																                  width = 4, height = 3, dpi = 200, useDingbats=FALSE)
																												    }
																			  }
					  summary_res = as.data.frame(t(as.data.frame(summary_res)))
					    summary_res$Cancer = rownames(summary_res)
					    summary_res = summary_res[order(summary_res$Coxph.Z), ]
						  res[["summary"]] = summary_res
						  summary_res = na.omit(summary_res)
						    summary_res$Cancer = factor(summary_res$Cancer, levels = summary_res$Cancer)

						    gg = reshape2::melt(summary_res[, seq(1,ncol(summary_res),2)], id.vars = "Cancer")
							  pvalue = reshape2::melt(summary_res[, c(seq(2,ncol(summary_res),2), ncol(summary_res))], id.vars = "Cancer")
							  gg$Pvalue = pvalue$value
							    gg$Symbol = ""
							    gg$Symbol[gg$Pvalue<0.05] = "*"
								  gg$Symbol[gg$Pvalue<0.01] = "**"
								  gg$Symbol[gg$Pvalue<0.001] = "***"
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
												    ggsave(paste0(outdir, "/BarView_summary_ICB_Signature.pdf"), p5,
														              width = 10, height = 8, dpi = 200, useDingbats=FALSE)
											  return(res)
}

