# source('~/Jobs/Project/ImmResponse/_Code/MHCI.TRAF3/R/theme_pubr.R')
#' Volcano View
#'
#' Volcano plot
#'
#' @docType methods
#' @name VolcanoView
#' @rdname VolcanoView
#'
#' @param df Data frame
#' @param x Colname of df specifying x-axis in Volcanno figure, 'logFC' (default).
#' @param y Colname of df specifying y-axis in Volcanno figure, 'adj.P.Val' (default).
#' @param Label Colname of df specifying labeled terms in Volcanno figure.
#' @param top Interger, the number of top significant terms to be labeled.
#' @param topnames Character vector, indicating interested terms to be labeled.
#' @param alpha Parameter in ggplot.
#' @param main Title of volcano figure.
#' @param xlab Label of x-axis in figure.
#' @param ylab Label of y-axis in figure.
#' @param filename Figure file name to create on disk. Default filename="NULL",
#' which means don't save the figure on disk.
#' @param width Width of figure.
#' @param height Height of figure.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(rra.gene_summary)
#' rra = ReadRRA(rra.gene_summary)
#' VolcanoView(rra, x = "LFC", y = "FDR", Label = "Official")
#'
#' @import ggrepel
#' @export

VolcanoView <- function(df, x = "logFC", y = "adj.P.Val",
                        Label = NA, top = 5, topnames = NULL,
                        mycolour=c("gray80", "#1f78b4", "#e31a1c", 
                                   "#33a02c", "#ff7f00", "#6a3d9a", 
                                   "#b15928", "#a6cee3", "#b2df8a", 
                                   "#fb9a99", "#fdbf6f", "#cab2d6", 
                                   "#ffff99"),
                        alpha=0.6, label.size = 3, force = 0.2, main = NULL, 
                        xlab = "Log2 Fold Change", ylab = "-Log10(Adjust.P)",
                        filename = NULL, width = 4, height = 2.5,
                        ...){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  gg = df[, c(x, y)]
  gg[, y] = -log10(gg[, y])
  
  if(!(top==0 & is.null(topnames))){
    if(!is.na(Label)) gg$Label = df[, Label] else gg$Label = rownames(gg)
    gg = gg[order(gg[,y], abs(gg[,x]), decreasing = TRUE), ]
    idx1 = idx2 = c()
    if(top>0){
      idx1 = which(gg[,x]>0)[1:min(top, sum(gg[,x]>0))]
      idx2 = which(gg[,x]<0)[1:min(top, sum(gg[,x]<0))]
    }
    idx = unique(c(idx1, idx2, which(gg$Label %in% unlist(topnames))))
    gg$Label = as.character(gg$Label)
    gg$Label[setdiff(1:nrow(gg), idx)] = ""
    # gg$Label = factor(gg$Label, levels = setdiff(unique(gg$Label), ""))
  }
  groups = rep(names(topnames), lengths(topnames))
  names(groups) = unlist(topnames)
  gg$group = "no"
  gg$group[gg$Label!=""] = "Others"
  idx = gg$Label %in% unlist(topnames)
  gg$group[idx] = groups[gg$Label[idx]]
  gg$group = factor(gg$group, levels = c("no", "Others", names(topnames)))
  gg$size = 0.6
  gg$size[gg$group!="no"] = 1.2
  gg = gg[order(gg$group, gg[,y]), ]
  #=========
  p = ggplot(gg, aes(x=gg[,x], y=gg[,y], label=Label, colour=group))
  p = p + geom_jitter(alpha=alpha, size = gg$size)
  p = p + scale_color_manual(values=mycolour, breaks=setdiff(levels(gg$group), c("no", "Others")))
  p = p + theme_pubr(14)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  # p = p + scale_fill_manual(values=mycolour)
  # if(y_cutoff>0 & y_cutoff<1)
  #   p = p + geom_hline(yintercept = -log10(y_cutoff), linetype = "dotted")
  # if(x_cutoff!=0)
  #   p = p + geom_vline(xintercept = c(-x_cutoff, x_cutoff), linetype = "dotted")
  p = p + labs(x=xlab, y=ylab, title=main, color=NULL)
  # p = p + annotate("text", color="#e41a1c", x = x_cutoff, y = max(gg[,y]), hjust = 0, vjust = 1,
  #                  label = paste("Up: ",dim(gg[gg$group=="up",])[1], sep=""))
  # p = p + annotate("text", color = "#377eb8", x = (-x_cutoff), y = max(gg[,y]), hjust = 1, vjust = 1,
  #                  label = paste("Down: ", dim(gg[gg$group=="down",])[1], sep=""))
  if(!(top==0 & is.null(topnames)))
    p = p + ggrepel::geom_text_repel(force = force, size = label.size,
                                     segment.color = 'grey50', segment.size = 0.2)
  # p = p + theme(legend.position = "none")
  
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, width = width, height = height, units = "in", ...)
  }
  return(p)
}



#' Volcano View
#'
#' Volcano plot
#'
#' @docType methods
#' @name VolcanoView2
#' @rdname VolcanoView2
#'
#' @param df Data frame
#' @param x Colname of df specifying x-axis in Volcanno figure, 'logFC' (default).
#' @param y Colname of df specifying y-axis in Volcanno figure, 'adj.P.Val' (default).
#' @param Label Colname of df specifying labeled terms in Volcanno figure.
#' @param top Interger, the number of top significant terms to be labeled.
#' @param topnames Character vector, indicating interested terms to be labeled.
#' @param x_cutoff Cutoff of x-axis.
#' @param y_cutoff Cutoff of y-axis.
#' @param alpha Parameter in ggplot.
#' @param main Title of volcano figure.
#' @param xlab Label of x-axis in figure.
#' @param ylab Label of y-axis in figure.
#' @param filename Figure file name to create on disk. Default filename="NULL",
#' which means don't save the figure on disk.
#' @param width Width of figure.
#' @param height Height of figure.
#' @param ... Other available parameters in ggsave.
#'
#' @return An object created by \code{ggplot}, which can be assigned and further customized.
#'
#' @author Wubing Zhang
#'
#' @examples
#' data(rra.gene_summary)
#' rra = ReadRRA(rra.gene_summary)
#' VolcanoView(rra, x = "LFC", y = "FDR", Label = "Official")
#'
#' @import ggrepel
#' @export

VolcanoView2 <- function(df, x = "logFC", y = "logP",
                         Label = NA, top = 5, bottom=0, topnames = NULL,
                         mycolour=c("gray80", "#ca0020", "#0571b0"),
                         x_cutoff = NULL, y_cutoff = NULL,
                         alpha=0.6, force = 0.2, labelsize = 2, main = NULL,
                         xlab = "Log2 Fold Change", ylab = "-Log10(Adjust.P)",
                         filename = NULL, width = 4, height = 2.5, ...){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  gg = df[, c(x, y)]
  # gg[, y] = -log10(gg[, y])
  
  # cutoff issues
  if(length(x_cutoff)==0) x_cutoff = c(0,0)
  if(length(y_cutoff)==0) y_cutoff = 0
  if(length(x_cutoff)==1) x_cutoff = sort(c(-x_cutoff, x_cutoff))
  if(length(x_cutoff)>2) x_cutoff = x_cutoff[1:2]
  if(length(y_cutoff)>1) y_cutoff = y_cutoff[1]
  gg$group = "no"
  gg$group[gg[,x]>x_cutoff[2] & gg[,y]>y_cutoff] = "up"
  gg$group[gg[,x]<x_cutoff[1] & gg[,y]>y_cutoff] = "down"
  
  # label issues
  if(!(top==0 & bottom==0 & is.null(topnames))){
    if(!is.na(Label)) gg$Label = df[, Label] else gg$Label = rownames(gg)
    gg = gg[order(gg[,y], abs(gg[,x]), decreasing = TRUE), ]
    idx1 <- idx2 <- c()
    if(top>0) idx1 = which(gg[,x]>0)[1:min(top, sum(gg[,x]>0))]
    if(bottom>0) idx2 = which(gg[,x]<0)[1:min(top, sum(gg[,x]<0))]
    idx = unique(c(idx1, idx2, which(gg$Label %in% topnames)))
    gg$Label = as.character(gg$Label)
    gg$Label[setdiff(1:nrow(gg), idx)] = ""
    # gg$Label = factor(gg$Label, levels = setdiff(unique(gg$Label), ""))
  }
  
  # color issues
  names(mycolour) = c("no", "up", "down")
  gg$color = mycolour[gg$group]
  gg$color[gg$Label!=""] = "black"
    
  gg = gg[order(-gg[,y]), ]
  p = ggplot(gg, aes(x=gg[,x], y=gg[,y], label=Label, fill=group))
  p = p + geom_point(shape = 21, alpha=alpha, color = gg$color)
  p = p + scale_fill_manual(values=mycolour)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  p = p + labs(x=xlab, y=ylab, title=main, color=NULL)
  p = p + theme_pubr(14)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  
  if(!(top==0 & bottom==0 & is.null(topnames)))
    p = p + ggrepel::geom_text_repel(force = force, 
                                     family = "Helvetica",
                                     size = labelsize,
                                     segment.color = 'grey50', 
                                     segment.size = 0.2)
  p = p + theme(legend.position = "none")
  
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, width = width, height = height, units = "in", ...)
  }
  return(p)
}





VolcanoView3 <- function(df, x = "log2FC", y = "padj", Label = NA, 
                         topnames = NULL, alpha=0.6, force = 0.2, 
                         labelsize = 2, dotsize = labelsize, textsize = labelsize, 
                         labelcolor = "#ca0020", main = NULL, 
                         xlab = expression(Log[2]*FC), ylab = expression(-Log[10]*FDR),
                         filename = NULL, width = 4, height = 2.5,
                         ...){
  requireNamespace("ggrepel", quietly=TRUE) || stop("need ggrepel package")
  gg = df[, c(x, y)]
  gg[, y] = -log10(gg[, y])
  colnames(gg) = c("X", "Y")
  gg = as.data.frame(gg, stringsAsFactors = FALSE)
  
  # label issues
  if(!is.na(Label)) gg$Label = df[, Label] else gg$Label = rownames(gg)
  gg$Label = as.character(gg$Label)
  idx = gg$Label %in% topnames
  gg$Label[!idx] = "None"
  mycolour = rep(labelcolor, length(topnames))
  mysize = rep(dotsize, length(topnames))
  names(mycolour) = names(mysize) = topnames
  
  gg = gg[order(-gg$Y), ]
  gg1 = gg[gg$Label=="None", ]
  gg2 = gg[gg$Label!="None", ]
  p = ggplot()
  p = p + geom_jitter(aes(x=X, y=Y), shape = ".", color = "gray50", size = 1, 
                      alpha=alpha, data = gg1)
  p = p + geom_vline(xintercept = 0, linetype = "dashed", color = "gray60")
  p = p + labs(x=xlab, y=ylab, title=main, color=NULL, size = NULL)
  p = p + geom_jitter(aes(x=X, y=Y, color = Label, size = Label), data = gg2)
  p = p + scale_color_manual(values=mycolour)
  p = p + scale_size_manual(values=mysize)
  p = p + ggrepel::geom_text_repel(aes(x=X, y=Y, label=Label),
                                   force = force,
                                   size = textsize,
                                   family = "Helvetica",
                                   segment.color = 'grey50', 
                                   segment.size = 0.5,
                                   data = gg2)
  p = p + theme_bw(base_line_size = NA)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  p = p + theme(legend.position = "none")
  
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, width = width, height = height, units = "in", ...)
  }
  return(p)
}