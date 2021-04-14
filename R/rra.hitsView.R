rra.hitsView <- function(rra, neg = TRUE, top = 5, genelist=c(),
                         color = "red", palette = 1, legend.pos = c(0.7,0.7)){
  require(ggplot2)
  require(ggrepel)
  require(ggView)
  if(neg) idx = grep("neg", colnames(rra))
  else idx = grep("pos", colnames(rra))
  gg = rra[, c(1, idx)]
  colnames(gg) = gsub("pos.|neg.", "", colnames(gg))
  gg = gg[!grepl("Zhang_|Ctrl", rownames(gg), ignore.case = TRUE), ]
  gg$Gene = rownames(gg)
  gg$lfc[gg$lfc< 0] = rnorm(sum(gg$lfc< 0), 0, 0.15)
  gg$rank = rank(-gg$lfc)
  gg = gg[order(gg$rank), ]
  gg$group = "Others"
  gg$group[gg$lfc>1] = "Selected"
  idx = rank(gg$lfc)>nrow(gg)-top | rownames(gg)%in%genelist
  gg$group[idx] = gg$Gene[idx]
  gg$group = factor(gg$group, levels = c(gg$Gene[idx], "Selected", "Others"))
  # Rank figure
  p = ggplot(gg, aes(rank, lfc, color = group, size = group))
  p = p + geom_jitter(alpha = 0.8)
  p = p + scale_color_manual(values = c(getCols(gg$Gene[idx], palette = palette), 
                                        "Selected"=color, "Others"="gray70"),
                             breaks = gg$Gene[idx])
  p = p + scale_size_manual(values = c(rep(1.5, sum(idx)), 0.4, 0.4),
                            breaks = gg$Gene[idx])
  p = p + xlim(0, max(gg$rank)) + ylim(min(gg$lfc), max(gg$lfc)+0.1)
  p = p + labs(x = "Rank", y = "Log2(Fold change)", color = NULL, size = NULL)
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 12, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=14),
                axis.text = element_text(colour="gray10"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + theme(legend.position=legend.pos, legend.key.height = unit(4.5, "mm"))
  p
}