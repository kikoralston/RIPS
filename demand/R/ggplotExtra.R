# **************************************************************** 
# ***** Custom functions to create special plots with ggplot *****
# **************************************************************** 


# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(..., grobs = list(...),
                                       ncol = length(grobs), nrow = 1, 
                                       position = c("bottom", "right")) {
  # Goal: share a legend between multiple plots that do not also share axes
  # Args:
  #   ...: ggplots
  #   grobs: list with ggplots
  #   ncol: number of columns in grid
  #   nrow: number of rows in grid
  #   position: position of legend
  #
  # Returns:
  #   ...
  plots <- grobs
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = 
                       arrangeGrob(do.call(arrangeGrob, gl), legend, ncol = 1,
                                   heights = unit.c(unit(1, "npc") - 
                                                      lheight, lheight)),
                     "right" = 
                       arrangeGrob(do.call(arrangeGrob, gl), legend, ncol = 2,
                                   widths = unit.c(unit(1, "npc") - 
                                                     lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
}
