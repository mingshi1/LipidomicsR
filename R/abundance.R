#' Summarize Abundance Data
#'
#' This function summarizes abundance data based on specified groups.
#'
#' @param data The data frame containing abundance data.
#' @param group A vector specifying the group membership for each sample.
#' @param .summary A function to summarize abundance data within each group (default is mean).
#' @param .error A function to compute error measures within each group (default is standard deviation).
#' @return A data frame summarizing abundance data by species, group, lipid type, carbon number, and unsaturation.
#' @importFrom reshape2 melt
#' @import dplyr
#' @export
abundance.summary <- function(data, group, .summary = function(x) mean(x), .error = function(x) sd(x)) {
  data$species <- rownames(data)
  data <- sepclass(data, 'all')
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
#' @param by A data frame specifying additional variables for the analysis.
#' @return A list containing statistical analysis results for each species.
#' @importFrom reshape2 melt
#' @import dplyr
#' @importFrom stats aov TukeyHSD
#' @importFrom broom tidy
#' @export
abundance.signif <- function(data, group, by = data.frame(group = c(), variable.1 = c(), variable.2 = c())) {
  data$species <- rownames(data)
  data <- sepclass(data, 'all')
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
      signif_list_individual <- list(data_anova, tukey)
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


