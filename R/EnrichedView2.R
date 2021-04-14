EnrcihedView2 <- function(enrichment, color="#08519c", size = 3,
                        legend.position = c(0.85,0.78), title = NULL){
  if(is(enrichment, "enrichResult")) enrichment = slot(enrichment, "result")
  if(is(enrichment, "gseaResult")) enrichment = slot(enrichment, "result")
  enrichment$LogFDR = -log10(enrichment$p.adjust)
  enrichment = enrichment[order(enrichment$p.adjust), ]
  enrichment$ID = factor(enrichment$ID, levels = enrichment$ID)
  idx = (max(enrichment$LogFDR)-enrichment$LogFDR) > enrichment$LogFDR
  enrichment$hjust = 1.1
  enrichment$hjust[idx] = -0.1
  p = ggplot(enrichment, aes(LogFDR, ID, label = Description))
  p = p + geom_point(size = size, color = color)
  p = p + xlim(0, NA)
  p = p + geom_text(aes(hjust = hjust))
  p = p + theme_bw() + theme(legend.position = legend.position, 
                             plot.title = element_text(hjust = 0.5))
  p = p + labs(x = expression(-log[10]*FDR), y = NULL, color = NULL, title = title)
  p
}

