#' @export
RankView2 <- function(gg, text = "Official", score = "LFC",
                      cutoff = 1, top = 5, genelist=c(),
                      color_pal = c("red", "gray70"),
                      legend.pos = c(0.7,0.7), title = NULL){
  require(ggplot2)
  require(ggrepel)
  colnames(gg)[colnames(gg)==score] = "LFC"
  # gg$LFC[gg$LFC< 0] = rnorm(sum(gg$LFC< 0), 0, 0.15)
  gg$Rank = rank(-gg$LFC)
  gg = gg[order(gg$Rank), ]
  gg$group = "Z2"
  gg$group[gg$LFC>cutoff] = "Z1"
  idx = gg$Rank<=top | gg[, text]%in%genelist
  gg$group[idx] = gg[idx, text]
  gg$group = factor(gg$group, levels = c(gg[idx, text], "Z1", "Z2"))

  pal = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
          "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
  if(length(color_pal)>length(levels(gg$group))){
    color_pal = color_pal[1:length(levels(gg$group))]
  }else if(length(color_pal)<length(levels(gg$group))){
    idx = length(levels(gg$group))-length(color_pal)
    print(idx)
    if(idx<6) color_pal = c(pal[sample(seq(2,10,2), idx)], color_pal)
    else color_pal = c(pal[sample(1:10, idx)], color_pal)
  }
  # Rank figure
  p = ggplot(gg, aes(Rank, LFC, color = group, size = group))
  p = p + geom_jitter(alpha = 0.6)
  p = p + scale_color_manual(values = color_pal,
                             breaks = setdiff(levels(gg$group), c("Z1", "Z2")))
  p = p + scale_size_manual(values = c(rep(2, sum(idx)), 0.4, 0.4),
                            breaks = setdiff(levels(gg$group), c("Z1", "Z2")))
  p = p + labs(y = "Log2(Fold change)", color = NULL, size = NULL, title = title)
  p = p + theme_bw(base_line_size = NA)
  p = p + theme(legend.key = element_rect(fill = "transparent"))
  p = p + theme(text = element_text(colour="black",size = 12, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=14),
                axis.text = element_text(colour="gray10"))
  p = p + theme(legend.position=legend.pos)
  p
}
