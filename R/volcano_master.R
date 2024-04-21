#' Calculate Fold Change and p-values for differential expression analysis
#'
#' This function calculates the log2 fold change (FC) and p-values for differential
#' expression analysis between two groups of data.
#'
#' @param compare1 A matrix or data frame representing the first group of data.
#' @param compare2 A matrix or data frame representing the second group of data.
#' @param p.adj Logical. Should p-values be adjusted for multiple testing?
#' @param method The method to use for p-value adjustment. Options: 'fdr', 'bonferroni', 'holm', etc.
#'
#' @return A data frame containing the log2 fold change and -log10 transformed p-values
#' for each row (e.g., genes, features) in the input data.
#'
#' @examples
#' compare1 <- matrix(rnorm(100), ncol = 10)
#' compare2 <- matrix(rnorm(100), ncol = 10)
#' result <- FC_P(compare1, compare2, p.adj = TRUE)
#'
#' @export
FC_P <- function(compare1, compare2, p.adj = FALSE, method = 'fdr') {
  avg1 <- rowMeans(compare1)
  avg2 <- rowMeans(compare2)
  FC <- avg2 / avg1
  log2FC <- log2(FC)
  p_value <- c()
  for (row in 1:nrow(compare1)) {
    if (all(as.numeric(compare1[row,]) == as.numeric(compare2[row, ]))) {
      p <- 1
    }
    else if (length(unique(as.numeric(compare1[row,])))==1 &
             length(unique(as.numeric(compare2[row,]))) == 1) {
      p <- NA
    }
    else {
      p <- t.test(compare1[row, ], compare2[row, ])$p.value
    }
    p_value <- c(p_value, p)
  }
  if (p.adj) {
    p_value <- p.adjust(p_value, method = method)
  }
  fcp_list <- data.frame('Log2FC' = log2FC, '-log10Pvalue' = -log10(p_value),
                         row.names = rownames(compare1), check.names = FALSE)
  return(fcp_list)
}

#' Generate Volcano Plot
#'
#' This function generates a volcano plot for differential expression analysis results.
#'
#' @param data The data frame containing the results of differential expression analysis.
#' @param FC.threshold The fold change threshold for determining significant changes (default is 2).
#' @param P.threshold The significance threshold for p-values (default is 0.05).
#' @param change.label The labels for differentially expressed genes (default is c('Up', 'Down', 'Notsig')).
#' @param point.size The size of data points in the plot (default is 2).
#' @param point.color The colors for differentially expressed genes (default is c('lightsalmon2', 'cadetblue', 'grey')).
#' @param x_scale_mannual Logical value indicating whether to manually specify x-axis limits (default is FALSE).
#' @param x.scale The manual limits for the x-axis (default is NULL).
#' @param y_scale_mannual Logical value indicating whether to manually specify y-axis limits (default is FALSE).
#' @param y.scale The manual limits for the y-axis (default is NULL).
#' @param linetype The line type for significance thresholds (default is 4).
#' @param line.alpha The transparency level for significance thresholds (default is 0.4).
#' @param line.color The color for significance thresholds (default is 'grey34').
#' @param line.size The size of significance thresholds (default is 1).
#' @param annotation.label The name of species that need to be annotated (default is NULL).
#' @param annotation.color The color of the annotation points (default is '#C82423')
#' @param text.size The size of gene labels (default is 1.5).
#' @param max.overlap The maximum number of overlapping labels allowed (default is 40).
#' @param title The title for the plot (default is NULL).
#' @param interact Logical value indicating whteher to generate interactive volcano plot.
#' @return A list containing the plot object, data frame with plotted points, and omitted data points. If interact = TRUE, the html object of the interactive plot will be also returned.
#' @importFrom ggplot2 ggplot geom_point geom_hline geom_vline scale_x_continuous scale_y_continuous theme_bw
#' @importFrom dplyr filter
#' @importFrom ggrepel geom_text_repel
#' @export
volcano <- function(data, x.scale, y.scale, interact=FALSE, FC.threshold=2, P.threshold = 0.05, change.label=c('Up','Down','Notsig'),
                    point.size = 2, point.color = c('lightsalmon2','cadetblue','grey'),
                    x_scale_mannual = FALSE,  y_scale_mannual = FALSE,
                    linetype= 4, line.alpha = .4, line.color= 'grey34', line.size = 1,
                    annotation.label = NULL, annotation.color ='#C82423',
                    text.size = 2.5, max.overlap = 40, title= NULL){
  volc.data <- data
  volc.data_out <- volc.data %>% dplyr::filter(.data$Log2FC == Inf | .data$Log2FC == -Inf | .data$Log2FC == 'NaN' | .data$`-log10Pvalue` == NA) %>%
    as.data.frame()
  volc.data <- volc.data[setdiff(rownames(volc.data), rownames(volc.data_out)),]

  volc.data$Change<-ifelse(volc.data$Log2FC>=log(FC.threshold,2) & volc.data$'-log10Pvalue'>= -log(P.threshold,10),change.label[1],
                           ifelse(volc.data$Log2FC<=-log(FC.threshold,2) & volc.data$'-log10Pvalue'>= -log(P.threshold,10),change.label[2],change.label[3]))
  change <- volc.data$Change
  volc.data <- as.data.frame(lapply(volc.data,as.numeric),check.names = F,row.names = rownames(volc.data))
  volc.data$Change <- change

  volc.data2 <- volc.data
  volc.data2[change.label[1],]<-c(0,-1,change.label[1])
  volc.data2[change.label[2],]<-c(0,-1,change.label[2])
  volc.data2[change.label[3],]<-c(0,-1,change.label[3])
  change <- volc.data2$Change
  volc.data2 <- as.data.frame(lapply(volc.data2,as.numeric),check.names = F,row.names = rownames(volc.data2))
  volc.data2$Change <- factor(change, levels = change.label)

  if (x_scale_mannual == FALSE) {
    x.scale <- c(-1.5*max(abs(volc.data2$Log2FC)),1.5*max(abs(volc.data2$Log2FC)))
  }
  if (y_scale_mannual == FALSE) {
    y.scale <- c(0,1.5*max(volc.data2$`-log10Pvalue`))
  }

  anno_data <- volc.data2[annotation.label,]
  vs<-ggplot(setdiff(volc.data2, anno_data),aes(.data$Log2FC,.data$`-log10Pvalue`))+
    geom_point(size = point.size, aes(color=.data$Change))+
    geom_point(data = anno_data, aes(x = .data$Log2FC, y = .data$`-log10Pvalue`), size = point.size, color = annotation.color)+
    geom_text_repel(data = anno_data, aes(x = .data$Log2FC, y = .data$`-log10Pvalue`, label = rownames(anno_data)),
                                          size = text.size, max.overlaps = max.overlap)+
    scale_color_manual(values= point.color)+
    scale_x_continuous(limits= x.scale)+
    scale_y_continuous(limits = y.scale)+
    geom_hline(yintercept = -log10(P.threshold),linetype= linetype,alpha= line.alpha,color=line.color, size = line.size)+
    geom_vline(xintercept = c(-log(FC.threshold,2), log(FC.threshold,2)),linetype= linetype,alpha= line.alpha,color=line.color, size = line.size)+
    xlab(expression("log"[2]*" Fold Change"))+
    ylab(expression("-log"[10]*" P-value"))+
    labs(title = title)+

    theme_bw() +
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())

  if (interact){
    vs2<-vs+ geom_point_interactive(data = volc.data2, size=point.size,aes(color=.data$Change,tooltip = rownames(volc.data2), data_id = volc.data2$Log2FC))+
      geom_point(data = anno_data, aes(x = .data$Log2FC, y = .data$`-log10Pvalue`), size = point.size, color = annotation.color)
    vs2 = girafe(print(vs2),options = list(
          opts_hover_inv(css="opacity:0.3;"),
          opts_hover(css="fill:red")
    ))
  }else {
    vs2 <- NULL
  }
  return(list(plot = vs, html = vs2, row = volc.data2, omit = volc.data_out))
}
