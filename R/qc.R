#' @title QCplot()
#' @description A function to exhibit the data quality between different samples, including
#' correlation heatmap, PCA, and quality boxplot.
#' @param data A dataframe storing concentration of lipids between different samples.
#' The column name should be the sample name and the row name should be the lipid type.
#' The class of column name and row name should be "character". The class of values should be
#' "numeric". The row names are recommended to be in a form like "PL(14:0/20:1)" or "LPL(16:1)".
#' It is highly recommended to input data after normalization.
#' @param ptype A vector to define picture types output. The input can include "heatmap",
#' "PCA" and "boxplot".
#'
#' @param group A vector defining which group the replicates belong to.
#' @param qcdt A dataframe containing internal labels and their abundance, which is used
#' to draw quality boxplot. If you don't want to draw the boxplot, the paramenter can
#' be ignored.
#' @param box.sample.name A vector containing the sample names, the length of "box.sample.name"
#' should equal the number of samples. This parameter can only change the sample name of
#' boxplot. The default values are the column names of the data input.
#' @param box.x The name of the x axis of the boxplot. The default value is "sample".
#' @param box.y The name of the x axis of the boxplot. The default value is "abundance".
#' @param box.title The picture title of the boxplot. The default value is "".
#' @param errorbar.show Whether show the errorbars of the boxplot or not. The default value is "".
#' @param group.col A vector containing the colors for each group. The length of "group.col" should equal the
#' number of the groups. If not input, the color will be default values.
#' @param outlie.col The color of outliers. The default value is NA (not show outliers).
#' @param outlie.shape The shape of outliers, which can be 1 - 25. The default value is NA (not show outliers).
#' @param heat.sample.name A vector containing the sample names, the length of "heat.sample.name"
#' should equal the number of samples. This parameter can only change the sample name of
#' heatmap. The default values are the column names of the data input.
#' @param heat.start.col The lightest color of the heatmap. The default value is "white".
#' @param heat.end.col The deepest color of the heatmap. The default value is "#3E8BCA".
#' @param heat.title The picture title of the heatmap. The default value is "Correlation Heatmap".
#' @param group.show Whether to show the classification of different groups of the heatmap. The default value is TRUE.
#' @param range.show Whether to show the range of PCA plot. The default value is TRUE.
#' @param range.alpha The transparency of the range in the PCA, only useful when range.show = TRUE. The default value is 0.25.
#' @param shape Whether to classify different groups by shape. The default value is TRUE.
#' @param pca.title The picture title of the PCA plot. The default value is "PCA Scores Plot".
#' @param point.size The size of points in the PCA plot. The default value is 3.5.
#' @param point.t.size The size of texts labeled on points. The default value is 1.5.
#' @param point.t.color The color of texts labeled on points. The default value is "grey25".
#' @param point.t.overlap To let the texts of points be shown without overlapping. The default value is 200.
#'
#' If you don't want to show texts labeled on points, please set "point.t.size=0" and "point.t.overlap=0".
#' @param combine  Whether to combine the three plots when ptype = "all". The default value is T. If combine = FALSE,
#' the three plots will be returned separately.
#'
#' @param marked Only useful when combine = T, it decides the labels on the top left of the picture. The default value
#' is  c("A", "B","C"). If you don't want to show, use "marked=c("", "","") ".
#' @param title.hjust Define the horizontal position of the picture title. The default value is 0.5.
#' @param title.vjust Define the vertical position of the picture title. The default value is 0.
#' @param title.size Define the size of the picture title. The default value is 15.
#' @param a.title.size Define the size of the axis title. The default value is 13.
#' @param a.text.size Define the size of the axis text. The default value is 8.
#' @param a.text.angle Define the angle of the X axis text of the boxplot, which can be 0 - 360. The default value is 45.
#' @param a.text.vjust Define the vertical position of the axis text. The default value is 1.
#' @param a.text.hjust Define the horizontal position of the axis text. The default value is 1.
#' @param l.text.size Define the size of the legend. The default value is 11.
#' @param l.title.size Define the size of the legend. The default value is 13.
#' @param margin Define the margin surrounding the plot area of each plot. It should be a vector whose length = 4. The default value is c(0.4,0.4,0.4,0.4).
#' @param unit The unit of the "margin" parameter, which can be "mm", "cm", "in", "pt", and "pc". The default value is "in".
#' @param cellsize The size of a cell in the heatmap. Default value is 8.
#' @param interactive Whether to get an interactive PCA plot or not. Default value is TRUE.
#' @return A correlation heatmap, PCA plot, quality boxplot, or a merged pictures containing
#' the above three.
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import ggiraph
#' @import ggforce
#' @importFrom stats quantile
#' @export
QCplot<-function(data,ptype,group,qcdt,box.sample.name=c("default"),box.x="sample",
                 box.y="abundance",box.title="",errorbar.show=TRUE,
                 group.col=c("default"),outlie.col=NA,outlie.shape=NA,heat.sample.name=c("default"),heat.start.col="white",heat.end.col="#3E8BCA",
                 heat.title="Correlation Heatmap",group.show=TRUE, range.show=TRUE,range.alpha=0.25,
                 shape=TRUE,pca.title="PCA Scores Plot",point.size=3.5,point.t.size=1.5,point.t.color="grey25",point.t.overlap=200,
                 marked= c("A", "B","C"),combine=TRUE,title.hjust=0.5,title.vjust=0,title.size=15,a.title.size=13,a.text.size=8,
                 a.text.angle=45,a.text.vjust=1,a.text.hjust=1,l.text.size=11,l.title.size=13,
                 margin=c(1,1,1,1),unit="in",cellsize=8,interactive=TRUE
){
  anno=data.frame(group=group)

  rownames(anno)=sample=colnames(data)
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col=col_vector[1:length(unique(anno$group))]

  plot_list <- list()
  ## boxplot
  if ("boxplot" %in% ptype) {
      if(length(box.sample.name)==1){
        box.sample.name=colnames(qcdt)
      }else{}
      if(length(group.col)==1){
        group.col=col
      }else{}
      colnames(qcdt)=box.sample.name
      qcdt_long=qcdt%>%
        pivot_longer(everything())
      qcdt_long$color=rep(anno$group,time=length(qcdt_long$name)/length(anno$group))
      if(errorbar.show){
        box=ggplot2::ggplot(qcdt_long,aes(x=as.factor(.data$name),y=.data$value,fill=.data$color))+
          scale_fill_manual(values = group.col)+
          stat_boxplot(geom = "errorbar",aes(ymin=quantile(.data$value)[1]),
                       width=0.8)+
          stat_boxplot(geom = "errorbar",aes(ymax=quantile(.data$value)[5]),
                       width=0.8)+
          geom_boxplot(outlier.shape = outlie.shape,outlier.colour = outlie.col)+
          stat_boxplot(aes(ymin=quantile(.data$value)[2],ymax=quantile(.data$value)[4]),outlier.shape = outlie.shape,outlier.colour = outlie.col)+

          labs(x=box.x,y=box.y,title = box.title,fill="")+
          theme_classic()+
          theme(plot.title = element_text(hjust = title.hjust,vjust=title.vjust,size = title.size),
                axis.text = element_text(size = a.text.size),axis.text.x=element_text(angle=a.text.angle,hjust=a.text.hjust,vjust=a.text.vjust),
                axis.title = element_text(size = a.title.size),
                legend.text = element_text(size = l.text.size),legend.title = element_text(size = l.title.size),
                plot.margin = unit(margin,unit))
      }else{
        box=ggplot2::ggplot(qcdt_long,aes(x=as.factor(.data$name),y=.data$value,fill=.data$color))+
          scale_fill_manual(values = group.col)+
          geom_boxplot(linetype="blank",outlier.shape = outlie.shape,outlier.colour = outlie.col)+
          stat_boxplot(aes(ymin=quantile(.data$value)[2],ymax=quantile(data$value)[4]),outlier.shape = outlie.shape,outlier.colour = outlie.col)+
          labs(x=box.x,y=box.y,title = box.title,fill="")+
          theme_classic()+
          theme(plot.title = element_text(hjust = title.hjust,vjust=title.vjust,size = title.size),
                axis.text = element_text(size = a.text.size),axis.text.x=element_text(angle=a.text.angle,hjust=a.text.hjust,vjust=a.text.vjust),
                axis.title = element_text(size = a.title.size),
                legend.text = element_text(size = l.text.size),legend.title = element_text(size = l.title.size),
                plot.margin = unit(margin,unit))
      }
    plot_list$boxplot <- box
}


  ## heatmap
  if ("heatmap" %in% ptype){
      if(length(heat.sample.name)==1){
        heat.sample.name=colnames(data)
      }else{}
      colnames(data)=heat.sample.name
      colnames(anno)=" "
      sampleDists <- dist(t(data))
      sampleDistMatrix <- as.matrix(sampleDists)
      # colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, heat.col)) )(255)
      colors <-colorRampPalette(c(heat.end.col, heat.start.col))(255)
      if(group.show){
        heat=pheatmap::pheatmap(sampleDistMatrix,
                                annotation_row = anno,
                                annotation_col = anno,
                                cluster_rows=FALSE,
                                cluster_cols = FALSE,
                                clustering_distance_rows = sampleDists,
                                clustering_distance_cols = sampleDists,
                                main = heat.title,
                                fontsize_row=a.text.size,fontsize_col=a.text.size,
                                fontsize=l.title.size,
                                cellheight = cellsize,cellwidth = cellsize,
                                col = colors)%>%ggplotify::as.ggplot()+
          theme(plot.title = element_text(hjust = title.hjust,size = title.size),
                legend.text = element_text(size = l.text.size),legend.title = element_text(size = l.title.size),
                plot.margin = unit(margin,unit))
      }else{
        heat=pheatmap::pheatmap(sampleDistMatrix,
                                cluster_rows=FALSE,
                                cluster_cols = FALSE,
                                clustering_distance_rows = sampleDists,
                                clustering_distance_cols = sampleDists,
                                main = heat.title,
                                fontsize_row=a.text.size,fontsize_col=a.text.size,
                                fontsize=l.title.size,
                                cellheight = cellsize,cellwidth = cellsize,
                                col = colors)%>%ggplotify::as.ggplot()+
          theme(plot.title = element_text(hjust = title.hjust,vjust=title.vjust,size = title.size),
                legend.text = element_text(size = l.text.size),legend.title = element_text(size = l.title.size),
                plot.margin = unit(margin,unit))
      }
    plot_list$heatmap <- heat
    dev.off()
}

  ## PCA
  if ("PCA" %in% ptype){
      pca1=prcomp(t(data),center = TRUE,scale. = TRUE)
      df1 <- pca1$x
      df1 <- as.data.frame(df1)
      df2 <- cbind(group,df1)
      colnames(anno)="group"
      summ1=summary(pca1)
      xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
      ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

      if(range.show){
        if(shape){
          p.pca1 <- ggplot(data = df1,aes(x = .data$PC1,y = .data$PC2,color = anno$group,shape = anno$group))+
            # stat_ellipse(aes(fill = anno$group),level = 0.95,
            #              type = "norm",geom = "polygon",alpha = range.alpha,color = NA)+
            ggforce::geom_mark_ellipse(aes(fill = anno$group,
                                           color = anno$group),alpha = range.alpha,color = NA) +
            geom_point(size = point.size)+
            labs(x = xlab1,y = ylab1,color = "",shape="",title = pca.title)+
            guides(fill = "none")+
            theme_bw()+
            scale_fill_manual(values = col)+
            scale_colour_manual(values = col)+
            scale_shape_manual(values = c(18:0))+
            theme(plot.title = element_text(hjust = title.hjust,vjust=title.vjust,size = title.size),
                  axis.text = element_text(size = a.text.size),axis.title = element_text(size = a.title.size),
                  legend.text = element_text(size = l.text.size),legend.title = element_text(size = l.title.size),
                  plot.margin = unit(margin,unit))+

            ggrepel::geom_text_repel(data=df1,aes(x = .data$PC1,y = .data$PC2,
                                                  label=rownames(df1)),colour = point.t.color,size=point.t.size,max.overlaps = point.t.overlap)
        }else{
          p.pca1 <- ggplot(data = df1,aes(x = .data$PC1,y = .data$PC2,color = anno$group))+
            ggforce::geom_mark_ellipse(aes(fill = anno$group,
                                           color = anno$group),alpha = range.alpha,color = NA) +
            # stat_ellipse(aes(fill = anno$group),level = 0.95,
            #              type = "norm",geom = "polygon",alpha = range.alpha,color = NA)+
            geom_point(size = point.size)+
            labs(x = xlab1,y = ylab1,color = "",title = pca.title)+
            guides(fill = "none")+
            theme_bw()+
            scale_fill_manual(values = col)+
            scale_colour_manual(values = col)+
            theme(plot.title = element_text(hjust = title.hjust,vjust=title.vjust,size = title.size),
                  axis.text = element_text(size = a.text.size),axis.title = element_text(size = a.title.size),
                  legend.text = element_text(size = l.text.size),legend.title = element_text(size = l.title.size),
                  plot.margin = unit(margin,unit))+
            ggrepel::geom_text_repel(data=df1,aes(x = .data$PC1,y = .data$PC2,
                                                  label=rownames(df1)),colour = point.t.color,size=point.t.size,max.overlaps = point.t.overlap)
        }
      }
      else{
        if(shape){
          p.pca1 <- ggplot(data = df1,aes(x = .data$PC1,y = .data$PC2,color = anno$group,shape = anno$group))+
            geom_point(size = point.size)+
            labs(x = xlab1,y = ylab1,color = "Color",shape="Shape",title = pca.title)+
            guides(fill = "none")+
            theme_bw()+
            scale_fill_manual(values = col)+
            scale_colour_manual(values = col)+
            scale_shape_manual(values = c(18:0))+
            theme(plot.title = element_text(hjust = title.hjust,vjust=title.vjust,size = title.size),
                  axis.text = element_text(size = a.text.size),axis.title = element_text(size = a.title.size),
                  legend.text = element_text(size = l.text.size),legend.title = element_text(size = l.title.size),
                  plot.margin = unit(margin,unit))+

            ggrepel::geom_text_repel(data=df1,aes(x = .data$PC1,y = .data$PC2,label=rownames(df1)),colour = point.t.color,size=point.t.size,max.overlaps = point.t.overlap)
        }else{
          p.pca1 <- ggplot(data = df1,aes(x = .data$PC1,y = .data$PC2,color = anno$group))+
            geom_point(size = point.size)+
            labs(x = xlab1,y = ylab1,color = "",title = pca.title)+
            guides(fill = "none")+
            theme_bw()+
            scale_fill_manual(values = col)+
            scale_colour_manual(values = col)+
            theme(plot.title = element_text(hjust = title.hjust,vjust=title.vjust,size = title.size),
                  axis.text = element_text(size = a.text.size),axis.title = element_text(size = a.title.size),
                  legend.text = element_text(size = l.text.size),legend.title = element_text(size = l.title.size),
                  plot.margin = unit(margin,unit))+

            ggrepel::geom_text_repel(data=df1,aes(x = .data$PC1,y = .data$PC2,
                                                  label=rownames(df1)),colour = point.t.color,size=point.t.size,max.overlaps = point.t.overlap)
        }
      }

      if (interactive == TRUE) {
        p.pca1 <- p.pca1+
                  geom_point_interactive(data=df2,aes(tooltip=rownames(df2),data_id =group),size=2)
        p.pca1 <- girafe(print(p.pca1),options = list(
                                                      opts_hover_inv(css="opacity:0.3;"),
                                                      opts_hover(css="fill:red")
                                                    ))
      }
      plot_list$PCA <- p.pca1
  }

  if ("boxplot" %in% ptype & "heatmap" %in% ptype & "PCA" %in% ptype) {
    if (combine == TRUE) {
      merged_plot=cowplot::plot_grid(plot_list$heatmap,plot_list$boxplot, plot_list$PCA, labels =marked,nrow = 2,ncol=2)
      return(merged_plot)
    }
  }
  return(plot_list)

}





