#' Plot Abundance Data
#'
#' This function generates bar plots for abundance data.
#'
#' @param data The data frame containing abundance data.
#' @param group A data frame containing grouping information.
#' @param by A data frame specifying additional variables for grouping and plotting (default is an empty data frame).
#' @param summary The method to summarise data, containing "Mean" and "Median". The default value is "Mean".
#' @param error The method to calculate errors, containing "SD" and "SE". The default value is "SD".
#' @param .width The width of lines in the plot (default is 0.5).
#' @param .position_dodge The position adjustment parameter for dodging lines (default is 0.5).
#' @param errorbar.width The width of error bars (default is 0.5).
#' @param .xlab The label for the x-axis (default is 'group').
#' @param .ylab The label for the y-axis (default is 'abundance').
#' @param sigf The type of significance in plots, can be "No Signal", "Star" or "P-value". The default value is "P-value".
#' @param axis.title.size The size of axis title text (default is 10).
#' @param axis.title.x.vjust The vertical adjustment parameter for x-axis title (default is 0).
#' @param axis.title.y.vjust The vertical adjustment parameter for y-axis title (default is 0).
#' @param axis.text.size The size of axis text (default is 10).
#' @param axis.line.size The size of axis lines (default is 0.5).
#' @param axis.tick.length The length of axis ticks (default is 0.2).
#' @param legend.title The title for the legend (default is an empty string).
#' @param legend.color The color palette for the legend (default is a set of predefined colors).
#' @param .legend.direction The direction of the legend ('horizontal' or 'vertical', default is 'vertical').
#' @param .legend.position The position of the legend ('top', 'bottom', 'left', 'right', or NULL, default is 'right').
#' @param main.size The size of plot titles (default is 10).
#' @return A list of ggplot objects, each representing a line plot for abundance data of a species.
#' @importFrom  reshape2 melt
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @importFrom car leveneTest
#' @importFrom rcompanion scheirerRayHare
#' @export
abundance.lineplot <- function(data, group, by = data.frame(group = c(), variable.1 = c(), variable.2 = c()),
                               summary="Mean",error="SD",.width = 0.5, .position_dodge = 0.5, errorbar.width = .5,
                               .xlab = 'group', .ylab = 'abundance',sigf="P-value",
                               axis.title.size = 10, axis.title.x.vjust = 0, axis.title.y.vjust = 0 ,axis.text.size = 10, axis.line.size = .5, axis.tick.length = 0.2,
                               legend.title = '', legend.color = 'Set2',
                               .legend.direction = 'vertical', .legend.position = 'right', main.size = 10){
  two_way_test<-function(data){
    p_value=c()
    for (i in c(1:dim(data)[1])){
      dt=t(data[i,])
      dt=cbind(dt,unlist(lapply(strsplit(rownames(dt),split = '_'), function(x){x[1]})))
      dt=cbind(dt,unlist(lapply(strsplit(rownames(dt),split = '_'), function(x){x[2]})))
      colnames(dt)=c("value","gr.1","gr.2")
      colnames(dt)[2]="gr.1"
      colnames(dt)[3]="gr.2"
      dt=as.data.frame(dt)
      dt$group=sub(pattern = '_\\d$', replacement = '', x = rownames(dt))
      dt$value=as.numeric(dt$value)
      aov=TRUE
      for (gr in unique(dt$group)){
        if (unlist(shapiro.test(as.numeric(dt$value[dt$group==gr]))[2])<0.05 || leveneTest(dt$value~dt$group, center=median)$`Pr(>F)`[1]<0.05){
          aov=FALSE
        }
      }
      if (aov){
        p_value=c(p_value,as.numeric(unlist(summary(aov(value ~ gr.1 + gr.2, dt)))[13]))
      }else{
        p_value=c(p_value,as.numeric(unlist(scheirerRayHare(value ~ gr.1 + gr.2, dt))[13]))
      }
    }
    data=cbind(p_value,data)
    return(data)
  }

  p_to_stars <- function(value) {
    if (value <= 0.05 & value >0.01) {
      p="*"
    }else if (value <= 0.01 & value > 0.001) {
      p="**"
    }else if (value <= 0.001) {
      p='***'
    }else{
      p=""
    }
    return(p)
  }

  if (summary=="Mean"){
    .summary = function(x) mean(x)
  }else{
    .summary = function(x) median(x)
  }
  if (error=="SD"){
    .error = function(x) sd(x)
  }else{
    .error = function(x) sd(x) / sqrt(length(x))
  }

  plot_list <- list()
  if ('variable.2' %in% colnames(by)){
    pdata=two_way_test(data)
  }
  data.summary=abundance.summary(data,group,.summary,.error)
  data.summary <- right_join(data.summary, by)


  for (specie in unique(data.summary$species)) {
    cdt <- data.summary %>% filter(.data$species == specie)
    if (sigf=="P-value" & 'variable.2' %in% colnames(by)){
      plabel=paste0("P-value: ",round(pdata$p_value[rownames(pdata)==specie],6))
    }else if(sigf=="Star" & 'variable.2' %in% colnames(by)){
      plabel=p_to_stars(round(pdata$p_value[rownames(pdata)==specie],6))
    }else{
      plabel=""
    }
    if ('variable.2' %in% colnames(by)){
      myplot <- ggplot(cdt,mapping = aes(x = .data$variable.2, y = .data$summary, group = .data$variable.1, colour =.data$variable.1,shape=.data$variable.1)) +
        geom_point(position = position_dodge(.position_dodge))+
        geom_line(position = position_dodge(.position_dodge), size = 1.5) +
        geom_errorbar(aes(ymin = .data$summary - .data$error , ymax = .data$summary + .data$error ), position = position_dodge(.position_dodge),width=errorbar.width) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab(label = .xlab) +
        ylab(label = .ylab) +
        labs(title = paste0(specie),color=legend.title,shape=legend.title) +
        theme_bw() +
        ylim(min(0,max(cdt$summary+cdt$error)), max(cdt$summary+cdt$error))+
        theme(

          axis.title = element_text(size = axis.title.size),
          axis.text = element_text(size = axis.text.size, color = 'black'),
          axis.line = element_line(size = axis.line.size, color = 'black'),
          axis.title.x = element_text(vjust = axis.title.x.vjust),
          axis.title.y = element_text(vjust = axis.title.y.vjust),
          axis.ticks.length = unit(axis.tick.length, 'cm'),
          plot.title = element_text(size = main.size, hjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.title = element_text(size = 8)
        ) +
        annotate("text", x=max(cdt$variable.2),y=0.01*(max(cdt$summary)-min(cdt$summary)), label = plabel,
                 parse = TRUE)+

        theme(
          legend.direction = .legend.direction,
          legend.position = .legend.position
        )
    }else if ('variable.2' %in% colnames(by) == FALSE) {
      myplot <- ggplot(cdt,mapping = aes(x = .data$variable.1, y = .data$summary)) +
        geom_point(position = position_dodge(.position_dodge))+
        geom_line(position = position_dodge(.position_dodge), size = 1.5) +
        geom_errorbar(aes(ymin = .data$summary - .data$error , ymax = .data$summary + .data$error ), position = position_dodge(.position_dodge),width=errorbar.width) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab(label = .xlab) +
        ylab(label = .ylab) +
        labs(title = paste0(specie),color=legend.title,shape=legend.title) +
        theme_bw() +
        ylim(min(0,max(cdt$summary+cdt$error)), max(cdt$summary+cdt$error))+
        theme(

          axis.title = element_text(size = axis.title.size),
          axis.text = element_text(size = axis.text.size, color = 'black'),
          axis.line = element_line(size = axis.line.size, color = 'black'),
          axis.title.x = element_text(vjust = axis.title.x.vjust),
          axis.title.y = element_text(vjust = axis.title.y.vjust),
          axis.ticks.length = unit(axis.tick.length, 'cm'),
          plot.title = element_text(size = main.size, hjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.title = element_text(size = 8)
        ) +
        theme(
          legend.direction = .legend.direction,
          legend.position = .legend.position
        )

    }

    myplotl <- list(myplot)
    names(myplotl) <- specie
    plot_list <- append(plot_list, myplotl)
  }
  return(plot_list)
}
