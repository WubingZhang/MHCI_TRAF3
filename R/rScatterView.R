#' @export
rScatterView <- function(genescore, ylimit = c(0, NA),
                         cutoff = 1, top = 5, genelist=c(),
                         color_pal = c("red", "gray70"),
                         labelsize = 4, title = NULL){
  gg = data.frame(Label = names(genescore), Score = genescore,
                  Rank = sample(1:length(genescore), length(genescore)),
                  stringsAsFactors = FALSE)
  if(!is.na(ylimit[1])) gg = gg[gg$Score>ylimit[1], ]
  if(!is.na(ylimit[2])) gg = gg[gg$Score<ylimit[2], ]
  gg$group = "Z2"
  gg$group[gg$Score>cutoff] = "Z1"
  idx = rank(gg$Score)>(nrow(gg)-top) | gg$Label%in%genelist
  gg$group[idx] = gg$Label[idx]
  tmp = seq(1, nrow(gg), length.out = sum(idx)+2)[-c(1,sum(idx)+2)]
  gg$Rank[idx] = sample(tmp, sum(idx))
  gg$group = factor(gg$group, levels = c(gg$Label[idx], "Z1", "Z2"))
  gg$Label = as.character(gg$group)
  gg$Label[gg$Label %in% c("Z1", "Z2")] = ""
  pal = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
  if(length(color_pal)!=length(levels(gg$group))){
    idx = length(levels(gg$group))-length(color_pal)
    if(idx<6) color_pal = c(pal[sample(seq(2,10,2), idx)], color_pal)
    else color_pal = c(pal[sample(1:10, idx)], color_pal)
  }
  p = ggplot(gg, aes(x=Rank, y=Score, color=group, size=group))
  p = p + geom_jitter(alpha = 0.6)
  p = p + scale_color_manual(values = color_pal, guide = "none")
  p = p + scale_size_manual(values = c(rep(2, sum(idx)), 0.4, 0.4), guide = "none")
  p = p + ggrepel::geom_text_repel(aes(label = Label), force = 0.1,
                                   fontface = 'bold', size = labelsize,
                                   segment.color = 'grey50', segment.size = 0.1,
                                   segment.alpha = 0)
  p = p + labs(color = NULL, size = NULL, title = title)
  p = p + theme_bw(base_line_size = NA)
  p = p + theme(plot.title = element_text(hjust = 0.5))
  p
}
