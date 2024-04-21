#' Calculate Total Abundance Data
#'
#' This function calculates total abundance data based on specified groups.
#'
#' @param data The data frame containing abundance data.
#' @return A data frame containing total abundance of different lipid type, carbon number, and unsaturation.
#' @import tidyverse
#' @export
total.abundance<- function(data){
  spdt=sepclass(data,"all")
  ldt=pivot_longer(spdt,cols = -c(.data$lipid_type,.data$carbon_number,.data$unsaturation), names_to = "species", values_to = "abundance")
  lp_sum=ldt%>%group_by(.data$lipid_type,.data$species)%>%summarise(sum=sum(.data$abundance))
  cb_sum=ldt%>%group_by(.data$carbon_number,.data$species)%>%summarise(sum=sum(.data$abundance))
  sat_sum=ldt%>%group_by(.data$unsaturation,.data$species)%>%summarise(sum=sum(.data$abundance))
  res=data.frame()
  for (l in 1:3){
    if (l==1){
      dt=lp_sum
      name="lipid_type="
    }else if(l==2){
      dt=cb_sum
      name="carbon_number="
    }else{
      dt=sat_sum
      name="unsaturation_rate="
    }
    for (i in 1:(dim(dt)[1]/dim(data)[2])){
      res=rbind(res,dt$sum[(dim(data)[2]*(i-1)+1):(dim(data)[2]*i)])
      rownames(res)[dim(res)[1]]=paste0(name,dt[(dim(data)[2]*(i-1)+1),1],"")
    }
  }
  colnames(res)=colnames(data)
  return(res)
}

#' Summarize Abundance Data
#'
#' This function summarizes abundance data based on specified groups.
#'
#' @param data The data frame containing abundance data.
#' @param istotal Logical. If is true, the total summary table of lipid type, carbon number and unsaturation rate will be generated.
#' @param group A vector specifying the group membership for each sample.
#' @param .summary A function to summarize abundance data within each group (default is mean).
#' @param .error A function to compute error measures within each group (default is standard deviation).
#' @return A data frame summarizing abundance data by species, group, lipid type, carbon number, and unsaturation.
#' @importFrom reshape2 melt
#' @import dplyr
#' @export
abundance.summary <- function(data, group, istotal = FALSE, .summary = function(x) mean(x), .error = function(x) sd(x)) {
  if (istotal == TRUE) {
    data <- total.abundance(data = data)
  }
  data$species <- rownames(data)
  data <- sepclass(data, 'all')
  if (istotal == TRUE) {
    data$lipid_type <- unlist(lapply(strsplit(x = rownames(data), split = '='), function(x) x[1]))
  }
  datalong <- melt(data, id.vars = c('species', 'lipid_type', 'carbon_number', 'unsaturation'), variable.name = 'samplename', value.name = 'abundance')
  datalong$group <- group[datalong$samplename]
  data_summary <-  datalong %>% group_by(.data$species, .data$group, .data$lipid_type, .data$carbon_number, .data$unsaturation) %>%
    summarise(summary = .summary(.data$abundance), error = .error(.data$abundance)) %>% as.data.frame()
  return(data_summary)
}

#' Analyze Abundance Significance
#'
#' This function performs statistical analysis to determine the significance of abundance data.
#'
#' @param data The data frame containing abundance data.
#' @param group A vector specifying the group membership for each sample.
#' @param istotal Logical. If is true, statistics based on total summary table of lipid type, carbon number and unsaturation rate will be generated.
#' @param by A data frame specifying additional variables for the analysis.
#' @return A list containing statistical analysis results for each species.
#' @importFrom reshape2 melt
#' @import dplyr
#' @importFrom stats aov TukeyHSD
#' @importFrom broom tidy
#' @export
abundance.signif <- function(data, group, istotal = FALSE, by = data.frame(group = c(), variable.1 = c(), variable.2 = c())) {
  if (istotal == TRUE) {
    data <- total.abundance(data = data)
  }
  data$species <- rownames(data)
  data <- sepclass(data, 'all')
  if (istotal == TRUE) {
    data$lipid_type <- unlist(lapply(strsplit(x = rownames(data), split = '='), function(x) x[1]))
  }
  datalong <- melt(data, id.vars = c('species', 'lipid_type', 'carbon_number', 'unsaturation'), variable.name = 'group', value.name = 'abundance')
  datalong$group <- group[datalong$group]
  datalong <- right_join(datalong, by)

  signif_list <- list(initialize = NULL)
  list_id <- 1
  for (specie in unique(datalong$species)) {
    cdt <- datalong %>% filter(.data$species == specie)
    if ('variable.2' %in% colnames(cdt) == TRUE) {
      data_anova <- cdt %>% aov(formula = abundance ~ variable.1 + variable.2 + variable.1:variable.2) %>% tidy() %>% as.data.frame()
      data_tukey <- cdt %>% aov(formula = abundance ~ variable.1 + variable.2 + variable.1:variable.2) %>% TukeyHSD()
      tukey <- data_tukey$`variable.1:variable.2` %>% as.data.frame()
      same <- c()
      for (i in 1:nrow(tukey)) {
        wsame <- length(unique(unlist(lapply(strsplit(strsplit(rownames(tukey)[i],'-')[[1]],':'), function(x) x[1]))))==1
        same <- c(same, wsame)
      }

      tukey <- tukey[same,]
      signif_list_individual <- list(data_anova, tukey)
      names(signif_list_individual) <- c('anova','tukey')
    } else if ('variable.2' %in% colnames(cdt) == FALSE) {
      data_anova <- cdt %>% aov(formula = abundance ~ variable.1) %>% tidy() %>% as.data.frame()
      v1_tukey <- TukeyHSD(aov(formula = abundance ~ variable.1, data = cdt))[[1]]
      signif_list_individual <- list(data_anova, v1_tukey)
      names(signif_list_individual) <- c('anova', 'tukey')
    }
    signif_list[[list_id]] <- signif_list_individual
    names(signif_list)[[list_id]] <- specie
    list_id <- list_id + 1
  }
  return(signif_list)
}

#' Plot Abundance Data
#'
#' This function generates bar plots for abundance data.
#'
#' @param data.summary A data frame containing summarized abundance data.
#' @param by A data frame specifying additional variables for grouping and plotting (default is an empty data frame).
#' @param .width The width of bars in the plot (default is 0.5).
#' @param .position_dodge The position adjustment parameter for dodging bars (default is 0.5).
#' @param errorbar.type The type of error bars to be plotted ('up', 'down', or 'both', default is 'up').
#' @param errorbar.width The width of error bars (default is 0.5).
#' @param .xlab The label for the x-axis (default is 'group').
#' @param .ylab The label for the y-axis (default is 'abundance').
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
#' @return A list of ggplot objects, each representing a bar plot for abundance data of a species.
#' @importFrom reshape2 melt
#' @import dplyr
#' @import ggplot2
#' @export
abundance.plot <- function(data.summary, by = data.frame(group = c(), variable.1 = c(), variable.2 = c()),
                           .width = 0.5, .position_dodge = 0.5, errorbar.type = 'up', errorbar.width = .25,
                           .xlab = 'group', .ylab = 'abundance',
                           axis.title.size = 10, axis.title.x.vjust = 0, axis.title.y.vjust = 0 ,axis.text.size = 10, axis.line.size = .5, axis.tick.length = 0.2,
                           legend.title = '', legend.color = c("#A1A9D0", "#F0988C", "#B883D4", "#9E9E9E", "#CFEAF1", "#C4A5DE", "#F6CAE5", "#96CCCB"),
                           .legend.direction = 'vertical', .legend.position = 'right', main.size = 10){
  plot_list <- list()

  if (errorbar.type == 'up') {
    errordown <- 0; errorup <- 1
  }else if (errorbar.type == 'down'){
    errordown <- 1 ; errorup <- 0
  }else if (errorbar.type == 'both'){
    errordown <- 1 ; errorup <- 1
  }

  data.summary <- right_join(data.summary, by)

  for (specie in unique(data.summary$species)) {
  cdt <- data.summary %>% filter(.data$species == specie)
  if ('variable.2' %in% colnames(by)){
    myplot <- ggplot(cdt)+
              geom_bar(mapping = aes(x= .data$variable.1, y = .data$summary, group = .data$variable.2, fill = .data$variable.2), colour='black',
                            stat="identity", position=position_dodge(.position_dodge),width = .width)+
              scale_fill_manual(name =  legend.title, values = legend.color)+
              geom_errorbar(aes(x = .data$variable.1, group = .data$variable.2, ymin = .data$summary - .data$error * errordown, ymax = .data$summary + .data$error * errorup), color ='black' ,width = errorbar.width,
                           stat = 'identity', position = position_dodge(.position_dodge))+
              scale_y_continuous(expand = c(0,0))+
              xlab(label = .xlab)+
              ylab(label = .ylab)+
              labs(title = paste0(specie))+
              theme_bw() +
              theme(axis.title = element_text(size= axis.title.size),axis.text = element_text(size= axis.text.size,color='black'),axis.line = element_line(size = axis.line.size, color='black'),
                   axis.title.x = element_text(vjust = axis.title.x.vjust),  axis.title.y = element_text(vjust = axis.title.y.vjust), axis.ticks.length = unit(axis.tick.length, 'cm'),
                   plot.title = element_text(size = 10, hjust = 0.5),
                   panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                   legend.title = element_text(size=8))+
              theme(legend.direction = .legend.direction, legend.position = .legend.position)
  }else if ('variable.2' %in% colnames(by) == FALSE) {
    myplot <- ggplot(cdt)+
              geom_bar(mapping = aes(x= .data$variable.1, y = .data$summary), colour='black',
                       stat="identity", position=position_dodge(.position_dodge),width = .width)+
              scale_fill_manual(values = legend.color)+
              geom_errorbar(aes(x = .data$variable.1, ymin = .data$summary - .data$error * errordown, ymax = .data$summary + .data$error * errorup), color ='black' ,width = errorbar.width,
                            stat = 'identity', position = position_dodge(.position_dodge))+
              scale_y_continuous(expand = c(0,0))+
              xlab(label = .xlab)+
              ylab(label = .ylab)+
              labs(title = paste0(specie))+
              theme_bw() +
              theme(axis.title = element_text(size= axis.title.size),axis.text = element_text(size=axis.text.size,color='black'),axis.line = element_line(size = axis.line.size, color='black'),
                    axis.title.x = element_text(vjust = axis.title.x.vjust),  axis.title.y = element_text(vjust = axis.title.y.vjust), axis.ticks.length = unit(axis.tick.length, 'cm'),
                    plot.title = element_text(size = 10, hjust = 0.5),
                    panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                    legend.title = element_text(size=8))+
              theme(legend.direction = .legend.direction, legend.position = NULL)

  }

  myplotl <- list(myplot)
  names(myplotl) <- specie
  plot_list <- append(plot_list, myplotl)
  }
  return(plot_list)
}

#' Plot Abundance Data
#'
#' This function generates line plots for abundance data.
#'
#' @param data.summary A data frame containing summarized abundance data.
#' @param by A data frame specifying additional variables for grouping and plotting (default is an empty data frame).
#' @param .width The width of lines in the plot (default is 0.5).
#' @param .position_dodge The position adjustment parameter for dodging lines (default is 0.5).
#' @param errorbar.width The width of error bars (default is 0.5).
#' @param .xlab The label for the x-axis (default is 'group').
#' @param .ylab The label for the y-axis (default is 'abundance').
#' @param axis.title.size The size of axis title text (default is 10).
#' @param axis.title.x.vjust The vertical adjustment parameter for x-axis title (default is 0).
#' @param axis.title.y.vjust The vertical adjustment parameter for y-axis title (default is 0).
#' @param axis.text.size The size of axis text (default is 10).
#' @param axis.line.size The size of axis lines (default is 0.5).
#' @param axis.tick.length The length of axis ticks (default is 0.2).
#' @param legend.title The title for the legend (default is an empty string).
#' @param color.style The style of colors, which can be "Nature", "Science", "Lancet", "JCO", "D3", "IGV", "Star Trek", "Tron Legacy", "Rick and Morty", and "The Simpsons".
#' Default value is "Nature".
#' @param legend.color The color palette for the legend (default is a set of predefined colors).
#' @param .legend.direction The direction of the legend ('horizontal' or 'vertical', default is 'vertical').
#' @param .legend.position The position of the legend ('top', 'bottom', 'left', 'right', or NULL, default is 'right').
#' @param main.size The size of plot titles (default is 10).
#' @return A list of ggplot objects, each representing a line plot for abundance data of a species. In addition, the list includes a line plot for the sum of abundance of each class of lipid type, unsaturation, and carbon number.
#' @importFrom  reshape2 melt
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import ggsci
#' @importFrom car leveneTest
#' @importFrom rcompanion scheirerRayHare
#' @export
abundance.lineplot <- function(data.summary, by = data.frame(group = c(), variable.1 = c(), variable.2 = c()),
                               .width = 0.5, .position_dodge = 0, errorbar.width = .5,
                               .xlab = 'group', .ylab = 'abundance',
                               axis.title.size = 10, axis.title.x.vjust = 0, axis.title.y.vjust = 0 ,axis.text.size = 10, axis.line.size = .5, axis.tick.length = 0.2,
                               legend.title = '', legend.color = 'Set2', color.style="Nature",
                               .legend.direction = 'vertical', .legend.position = 'right', main.size = 10){
  plot_list <- list()
  data.summary <- right_join(data.summary, by)

  for (specie in unique(data.summary$species)) {
    cdt <- data.summary %>% filter(.data$species == specie)
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
        ylim(min(0,min(cdt$summary-cdt$error)), max(cdt$summary+cdt$error))+
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
        ylim(min(0,min(cdt$summary-cdt$error)), max(cdt$summary+cdt$error))+
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
    if (color.style=="Nature"){
      myplot <- myplot + scale_color_npg()
    }else if (color.style=="Science"){
      myplot <- myplot + scale_color_aaas()
    }else if (color.style=="Lancet"){
      myplot <- myplot + scale_color_lancet()
    }else if (color.style=="JCO"){
      myplot <- myplot + scale_color_jco()
    }else if (color.style=="D3"){
      myplot <- myplot + scale_color_d3()
    }else if (color.style=="IGV"){
      myplot <- myplot + scale_color_igv()
    }else if (color.style=="Star Trek"){
      myplot <- myplot + scale_color_startrek()
    }else if (color.style=="Tron Legacy"){
      myplot <- myplot + scale_color_tron()
    }else if (color.style=="Rick and Morty"){
      myplot <- myplot + scale_color_rickandmorty()
    }else{
      myplot <- myplot + scale_color_simpsons()
    }
    myplotl <- list(myplot)
    names(myplotl) <- specie
    plot_list <- append(plot_list, myplotl)
  }
  return(plot_list)
}
