ggGroup <- function(
  x = NULL, 
  y = NULL, 
  xlabel = NULL, 
  ylabel = NULL, 
  groupOrder = NULL,
  groupSort = FALSE,
  size = 1,  
  baseSize = 10,
  ridgeScale = 1, 
  ratioYX = NULL,
  alpha = 1,
  title = "", 
  pal = paletteDiscrete(values=x, set = "stallion"),
  addBoxPlot = TRUE,
  plotAs = "ridges",
  coord_flip = FALSE,
  ...
  ){

  .validInput(input = x, name = "x", valid = c("character"))
  .validInput(input = y, name = "y", valid = c("numeric"))
  .validInput(input = xlabel, name = "xlabel", valid = c("character", "null"))
  .validInput(input = ylabel, name = "ylabel", valid = c("character", "null"))
  .validInput(input = groupOrder, name = "groupOrder", valid = c("character", "null"))
  .validInput(input = groupSort, name = "groupSort", valid = c("boolean"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = ridgeScale, name = "ridgeScale", valid = c("numeric"))
  .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric", "null"))
  .validInput(input = alpha, name = "alpha", valid = c("numeric"))
  .validInput(input = title, name = "title", valid = c("character"))
  .validInput(input = pal, name = "pal", valid = c("character"))
  .validInput(input = addBoxPlot, name = "addBoxPlot", valid = c("boolean"))
  .validInput(input = plotAs, name = "plotAs", valid = c("character"))
  .validInput(input = coord_flip, name = "coord_flip", valid = c("boolean"))
  
  names(y) <- x
  dm <- stats::aggregate(y ~ names(y), FUN = mean)
  df <- data.frame(x, y)

  if(!is.null(groupOrder)){
    if(!all(x %in% groupOrder)){
      stop("Not all x values are present in groupOrder!")
    }
  }else{
    if(groupSort){
      groupOrder <- paste0(dm[,1])[order(dm[,2], decreasing= FALSE)]
    }else{
      if(tolower(plotAs) == "ridges"){
        groupOrder <- rev(gtools::mixedsort(unique(x)))
      }else{
        groupOrder <- gtools::mixedsort(unique(x))
      }
    }
  }

  df$x <- factor(df$x, groupOrder)
  
  p <- ggplot(df, aes(x = x, y = y, color = x)) +
      scale_color_manual(values = pal, guide = "none") + 
      scale_fill_manual(values = pal, guide = "none") +
      ggtitle(title)

  if(tolower(plotAs) == "ridges" | tolower(plotAs) == "ggridges"){
    if(!requireNamespace("ggridges", quietly = TRUE)){
      type <- "violin"
      message("ggridges is not available for plotting, continuing with geom_violin!")
      message("To install ggridges try : install.packages('ggridges')")
      p <- p + geom_violin(aes_string(fill="x"), alpha = alpha)
    }else{
      type <- "ridges"
      .requirePackage("ggridges", source = "cran")
      #p <- p + 
      #  stat_density_ridges(aes_string(x = "y", y = "x", fill = "x"), 
      #    quantile_lines = TRUE, quantiles = c(0.5), alpha = alpha, color = "black",
      #    scale = ridgeScale
      #  ) + scale_y_discrete(expand = c(0, 0))
      #   stat_density_ridges(
      #     aes_string(x = "y", y = "x", fill = "x"),
      #     quantile_lines = TRUE,
      #     alpha = alpha,
      #     geom = "density_ridges_gradient",
      #     calc_ecdf = TRUE,
      #     quantiles = c(0.5)
      # )
      val <- 1/length(unique(x))
      p <- p + geom_density_ridges(data = df,
        aes(x = y, y = x, color = x, fill = x), scale = ridgeScale,
        alpha = alpha, color = "black") + scale_y_discrete(expand = expansion(mult = c(0.01, val)))
    }
  }else{
    type <- "violin"
    p <- p + geom_violin(aes_string(x = "x", y = "y", color = "x", fill="x"), alpha = alpha)
  }
  
  if(addBoxPlot & type == "violin"){
    p <- p + geom_boxplot(size = size, outlier.size = 0, outlier.stroke = 0, fill = NA) 
  }

  if(type != "violin"){
    p <- p + theme_ArchR(baseSize = baseSize)
  }else{
    p <- p + theme_ArchR(xText90 = TRUE, baseSize = baseSize)
  }

  if(!is.null(ratioYX)){
    p <- p + coord_fixed(ratioYX, expand = TRUE)
  }

  if (!is.null(xlabel)) {
    if(type=="violin"){
      p <- p + xlab(xlabel)
    }else{
      p <- p + xlab(ylabel)
    }
  }
  
  if (!is.null(ylabel)) {
    if(type=="violin"){
      p <- p + ylab(ylabel)
    }else{
      p <- p + ylab(xlabel)
    }
  }
  
  p <- p + theme(legend.position = "bottom")

  if(!is.null(ratioYX)){
    attr(p, "ratioYX") <- ratioYX
  }
  
  if(coord_flip){
    p <- p + coord_flip()
  }


  return(p)

}
