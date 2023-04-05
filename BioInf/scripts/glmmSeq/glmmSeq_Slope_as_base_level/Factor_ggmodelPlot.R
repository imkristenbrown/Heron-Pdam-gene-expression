library(ggplot2)

Factor_ggmodelPlot <- function(object, geneName = NULL, x1var = NULL, x2var = NULL, 
          x2shift = NULL, xlab = NULL, ylab = geneName, plab = NULL, 
          title = geneName, logTransform = is(object, "GlmmSeq"), 
          shapes = 19, colours = "grey60", lineColours = "grey60", 
          markerSize = 1, fontSize = 12, alpha = 0.7, x2Offset = 5, 
          addPoints = TRUE, addModel = TRUE, modelSize = 4, modelColours = "blue", 
          modelLineSize = 1, modelLineColours = modelColours, addBox = FALSE, 
          ...) 
{
  if (!(is(object, "GlmmSeq") | is(object, "lmmSeq"))) {
    stop("object must be an output from glmmSeq or lmmSeq")
  }
  dfs <- formPlot(object, geneName, x1var, x2var, x2shift)
  df_long <- dfs[[1]]
  df_model <- dfs[[2]]
  xdiff <- dfs[[3]]
  x2shift <- dfs[[4]]
  modelData <- object@modelData
  maxX2 <- max(df_long$x2, na.rm = TRUE)
  df_long$x2 <- factor(df_long$x2)
  if (!is.null(x2var)) {
    x2labs <- levels(droplevels(factor(modelData[, x2var])))
  }
  xlim <- range(c(df_long$x, df_model$x), na.rm = TRUE)
  pval <- object@stats$pvals[geneName, , drop = FALSE]
  pval <- formatC(pval, digits = 2)
  lineColours <- rep_len(lineColours, maxX2)
  modelColours <- rep_len(modelColours, maxX2)
  colours <- rep_len(colours, maxX2)
  shapes <- rep_len(shapes, maxX2)
  if (addPoints) {
    p <- ggplot(df_long, aes_string(x = "x", y = "y", group = "id", 
                                    shape = "x2", color = "x2")) + geom_point(size = markerSize, 
                                                                              alpha = alpha, colour = colours[df_long$x2]) +
                                      geom_line(alpha = alpha) +
                                      geom_label(aes(label=id),position = position_dodge(width=.5),size=2)
  }
  else {
    p <- ggplot()
  }
  p <- p + theme_classic() + scale_shape_manual(values = shapes) + 
    labs(x = xlab, y = ylab, title = title) + theme(legend.position = "none", 
                                                    plot.margin = margin(7, 4, 14 + x2Offset, 4), text = element_text(size = fontSize), 
                                                    axis.text = element_text(colour = "black", size = fontSize - 
                                                                               1), ...) + scale_x_continuous(labels = modelData[,
                                                                                                                              x1var], breaks = if (x2shift < xdiff) {
                                                                                                                                modelData[, x1var]
                                                                                                                               }
                                                                                                            else {
                                                                                                              #modelData[, x1var] + (as.numeric(modelData[, x2var]) -1) * x2shift
                                                                                                              as.numeric(factor(modelData[, x1var])) + (as.numeric(modelData[, x2var])-1) * x2shift
                                                                                                              }) +
                                                                                              coord_cartesian(clip = "off") +
                                                                                                  scale_color_manual(values = lineColours)
  if (addBox) {
    p <- p + geom_boxplot(mapping = aes_string(x = "x", 
                                               y = "y", group = "x"), inherit.aes = FALSE, alpha = alpha * 
                            0.7, outlier.shape = NA, width = xdiff/6)
  }
  if (x2shift >= xdiff) {
    p <- p + geom_text(data = data.frame(label = x2labs, 
                                         #x = x2shift * (seq_along(x2labs) - 1) + xdiff/2, 
                                         x = x2shift * (seq_along(x2labs) - 1) + xdiff/2 + 1,
                                         y = min(c(modelData$LCI, df_long$y), na.rm = TRUE)), 
                       size = rel(4), mapping = aes_string(label = "label", 
                                                           x = "x", y = "y"), hjust = 0.5, nudge_x = 0, 
                       vjust = x2Offset, inherit.aes = FALSE) + theme(axis.title.x = element_text(vjust = -6))
  }
  if (addModel) {
    p <- p + annotate("line", x = df_model$x, y = df_model$y, 
                      group = df_model$group, size = modelLineSize, color = modelLineColours[as.numeric(df_model$group)]) + 
      annotate("errorbar", x = df_model$x, y = df_model$y, 
               color = modelLineColours[as.numeric(df_model$group)], 
               ymin = df_model$lower, ymax = df_model$upper, 
               width = xdiff/6, size = modelLineSize) + annotate("point", 
                                                                 x = df_model$x, y = df_model$y, shape = shapes[as.numeric(df_model$group)], 
                                                                 size = modelSize, color = modelColours[as.numeric(df_model$group)])
  }
  if (logTransform) 
    p <- p + scale_y_continuous(trans = "log10")
  if (is.null(plab)) 
    plab <- colnames(pval)
  ptext <- lapply(1:ncol(pval), function(i) {
    bquote("P"[.(plab[i])] * "=" * .(pval[, i]))
  })
  ptext <- bquote(.(paste(unlist(ptext), collapse = "*\", \"*")))
  p <- p + labs(subtitle = parse(text = ptext))
  return(p)
}


formPlot <- function(object, geneName, x1var, x2var, x2shift) {
  if (!x1var %in% colnames(object@modelData)) {
    stop("x1var must be a column name in object@modelData")}
  if (!is.null(x2var)) if (!x2var %in% colnames(object@modelData)) {
    stop("x2var must be a column name in object@modelData")}
  # if(ncol(object@modelData) > 2){
  #   stop("More than 2 variables in modelData")}
  
  maindata <- if (inherits(object, "GlmmSeq")) {
    object@countdata} else object@maindata
  if(! geneName %in% rownames(maindata)) {
    stop("geneName not found")}
  
  # Set up plotting data frame
  IDColumn <- object@vars$id
  id <- object@metadata[, IDColumn]
  y <- maindata[geneName, ]
  #x <- object@metadata[, x1var]
  x <- as.numeric(factor(object@metadata[, x1var]))
  xdiff <- diff(range(x, na.rm = TRUE))
  if (!is.null(x2var)) {
    x2 <- as.numeric(factor(object@metadata[, x2var]))
    nsegments <- length(unique(x)) -1
    if (is.null(x2shift)) {x2shift <- max(x, na.rm = TRUE) + xdiff / nsegments}
    x <- x + (x2-1) * x2shift
  } else {
    x2 <- 1
    x2shift <- -Inf
  }
  df_long <- data.frame(id, y, x, x2)
  
  # Set up model fit data
  modelData <- object@modelData
  preds <- object@predict[geneName, ]
  s <- nrow(modelData)
  modelx <- if (!is.null(x2var)) {
    #modelData[, x1var] + (as.numeric(modelData[, x2var])-1) * x2shift
    as.numeric(factor(modelData[, x1var])) + (as.numeric(modelData[, x2var])-1) * x2shift
  } else modelData[, x1var]
  df_model <- data.frame(x = modelx,
                         y = preds[1:s],
                         lower = preds[1:s +s],
                         upper = preds[1:s +s*2],
                         group = modelData[, x2var])
  if (is.null(x2var)) df_model$group <- 1
  
  return(list(df_long, df_model, xdiff, x2shift))
}
