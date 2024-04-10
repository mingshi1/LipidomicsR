#' @title heatmap()
#' @description A function to plot a heatmap based on abundance of lipids. A series
#' of custom functions can be realized such as dividing groups.
#' @param data A dataframe storing absolute concentration or PL% of lipids between different samples. If not, use cleanXpert() or noridx.output() + normalization.calculator()
#' to get normalized data.
#' The column name should be the sample name and the row name should be the lipid type.
#' The class of column name and row name should be "character". The class of values should be
#' "numeric". The row names are recommended to be in a form like "PL(14:0/20:1)" or "LPL(16:1)".
#' @param group A vector defining which group the replicates belong to.
#' @param cluster_col A boolean variable controlling whether to perfrom clustering to
#' column variables (lipid abundance) or not. The default value is TRUE.
#' @param cluster_row A boolean variable controlling whether to perfrom clustering to
#' row variables (lipid abundance) or not. The default value is TRUE.
#' @param sel.group A vector containing the group you want to show in the heatmap.
#' The input can be like c("WT","KO"). Default value is "default".
#' @param type Can accept 3 values: "lipid", "CB", or "sat". Default value is "lipid".
#'
#'     If type="lipid“, the heatmap will divide rownames based on lipid types.
#'
#'     If pattern="CB“, the heatmap will divide rownames based on carbon number.
#'
#'     If pattern="sat“, the heatmap will divide rownames based on the number of double bonds of a lipid type.
#' @param sel.type A vector controlling which types to show. If you only want to check data of "PA", "PC“, and "PE".
#' You can set type="lipid", sel.type=c("PA","PC","PE")
#' @param sel.row A vector controlling which types to show. If you set it as c("LPC(16:0)","PC(14:0/16:1)",
#' "PC(18:1/18:1)","PE(18:0/20:1)"), only their abundance will be shown.
#' @param annotation_legend A boolean controlling whether to show the figure legend. The default value is TRUE.
#' @param cellwidth The width of a cell in the heatmap. Default value is 20.
#' @param cellheight The height of a cell in the heatmap. Default value is 15.
#' @param gaps_row To customize positions of row gaps. Default value is c(0).
#' @param gaps_col To customize positions of column gaps. Default value is c(0).
#'
#' Notice: gaps_row and gaps_col are only useful when cluster=FALSE.
#' @param constract The constract of heatmap, default is 8.5, value range from 0 to 10.
#' @param labels_row A vector contains the labels of each row of the heatmap.
#' Default value is row names of dataframe input. It can be input like ("PE(20:1/20:1)","PS(16:0/18:1)","","","","LPA(18:0)")
#' @param labels_col A vector contains the labels of each column of the heatmap.
#' Default value is column names of dataframe input.
#' @param title The title of heatmap. Default value is "".
#' @param show_rownames Whether to show row names or not. Default value is T.
#' @param show_colnames Whether to show column names or not. Default value is T.
#' @param cellcolor The color range of cells in the heatmap. It should be input in a vector with
#' three color values, such as c("blue","black","yellow").
#' @param legend Whether to show legends or not. Default value is FALSE.
#' @param border_color Useful when border=T. Default value is NA.
#' @param border Whether to show borders or not. Default value is TRUE.
#' @param cutree_rows Useful when cluster=T. If cutree_rows=T, the rows of heatmap will be divided according to
#' clustering results. Default value is TRUE.
#' @param cutree_cols Useful when cluster=T. If cutree_cols=T, the rows of heatmap will be divided according to
#' clustering results. Default value is TRUE.
#' @param rtitle Row title of the heatmap. Default value is "group".
#' @param ctitle Column title of the heatmap. Default value is " ".
#' @param fontsize_row Fontsize of row labels. Default value is 12.
#' @param fontsize_col Fontsize of column labels. Default value is 12.
#' @param fontsize Fontsize of all labels. Default value is 8.
#'
#' @return A heatmap that is color-coded by abundance of lipids.
#' @examples
#' WT_1=rnorm(n=10,mean=0.4,sd=0.1)
#'WT_2=rnorm(n=10,mean=0.4,sd=0.1)
#'WT_3=rnorm(n=10,mean=0.4,sd=0.1)
#'WT_4=rnorm(n=10,mean=0.4,sd=0.1)
#'KO_1=rnorm(n=10,mean=0.8,sd=0.1)
#'KO_2=rnorm(n=10,mean=0.8,sd=0.1)
#'KO_3=rnorm(n=10,mean=0.8,sd=0.1)
#'KO_4=rnorm(n=10,mean=0.8,sd=0.1)
#'WT_treat_1=rnorm(n=10,mean=0.1,sd=0.1)
#'WT_treat_2=rnorm(n=10,mean=0.1,sd=0.1)
#'WT_treat_3=rnorm(n=10,mean=0.1,sd=0.1)
#'WT_treat_4=rnorm(n=10,mean=0.1,sd=0.1)
#'KO_treat_1=rnorm(n=10,mean=0.6,sd=0.1)
#'KO_treat_2=rnorm(n=10,mean=0.6,sd=0.1)
#'KO_treat_3=rnorm(n=10,mean=0.6,sd=0.1)
#'KO_treat_4=rnorm(n=10,mean=0.6,sd=0.1)
#'data=data.frame(WT_1,WT_2,WT_3,WT_4,KO_1,KO_2,KO_3,KO_4,
#'                WT_treat_1,WT_treat_2,WT_treat_3,WT_treat_4,
#'                KO_treat_1,KO_treat_2,KO_treat_3,KO_treat_4)
#'rownames(data)=c("LPC(16:0)","PC(14:0/16:1)","PC(18:1/18:1)","PE(18:0/20:1)",
#'                 "PS(20:1/20:1)","PI(16:0/16:1)","PC(18:0/18:1)","PA(16:0/16:1)",
#'                 "LPE(18:0)","PE(O-18:1/18:0)")
#'group=rep(c("WT","KO","WT_treat","KO_treat"),each=4)
#' heatmap(data,group)
#' @export
heatmap<-function(data,group,cluster_row=TRUE,cluster_col=TRUE,sel.group="default", constract=8.5,
                  type="lipid",sel.type="default",sel.row=c("default"),annotation_legend=TRUE,
                  cellwidth = 20, cellheight = 15,gaps_row = c(0),gaps_col=c(0),
                  labels_row = c("default"),labels_col = c("default"),title="",show_rownames=TRUE,
                  show_colnames = TRUE,cellcolor=c("blue","black","yellow"),
                  legend = TRUE,border_color = NA,border = FALSE,cutree_rows = 1,
                  cutree_cols = 1,rtitle="group",ctitle=" ",fontsize_row=12,fontsize_col=12,fontsize=8){
  bk <- c(seq(-10.1+constract,-0.01,by=0.01),0,seq(0.01,10.1-constract,by=0.01))
  cellcolor=c(colorRampPalette(colors = cellcolor[1:2])(length(bk)/2),colorRampPalette(colors =cellcolor[2:3])(length(bk)/2))
  rn=rownames(data)
  data=apply(data,2,as.numeric)
  rownames(data)=rn
  dt=sepclass(data,"all")

  if(sel.group!="default"){
    df=dt[,colnames(dt)%in%sel.group]
    if(type=="lipid"){
      df=df[df$lipid_type%in%sel.type,]
    }else if(type=="CB"){
      df=df[df$carbon_number%in%sel.type,]
    }else if(type=="sat"){
      df=df[df$unsaturation%in%sel.type,]
    }else{
      return("Please input right 'sel.type' or 'type' parameter.")
    }
    if(sel.type!="default"){
      df=df[rownames(df)%in%sel.type,]
    }else{}
    df=df[,4:(length(group)+3)]
  }else{df=dt[,4:(length(group)+3)]}
  anno_row=data.frame(group=group)

  rownames(anno_row)=colnames(df)
  anno=dt[,1:3]
  if(type=="lipid"){
    anno=anno[order(anno[,1]),]%>%as.data.frame()
    anno_col=anno[,1]%>%as.data.frame()
    rownames(anno_col)=rownames(anno)
    colnames(anno_col)="lipid type"
  }else if(type=="CB"){
    anno=anno[order(anno[,2]),]%>%as.data.frame()
    anno_col=anno[,2]%>%as.data.frame()
    rownames(anno_col)=rownames(anno)
    colnames(anno_col)="carbon number"
  }else if(type=="sat"){
    anno=anno[order(anno[,3]),]%>%as.data.frame()
    anno_col=anno[,3]%>%as.data.frame()
    rownames(anno_col)=rownames(anno)
    colnames(anno_col)="unsaturation"
  }else{
    return("Please input right 'type' parameter.")
  }
  colnames(anno_col)="species"
  df=df[rownames(anno_col),]
  if(length(labels_row)==1){
    labels_row = rownames(df)
  }else{}
  if(length(labels_col)==1){
    labels_col = colnames(df)
  }else{}
  if (ctitle!=""){
    colnames(anno_col)=ctitle
  }else{}
  colnames(anno_col)=ctitle
  colnames(anno_row)=rtitle
  # if(cluster_row | cluster_col){
    heat=pheatmap::pheatmap(as.matrix(df),cluster_row = cluster_row,cluster_cols = cluster_col,annotation_legend=annotation_legend,
                            cellwidth = cellwidth, cellheight = cellheight,labels_row = labels_row,labels_col = labels_col,
                            main = title,show_rownames=show_rownames,show_colnames = show_colnames,
                            color = cellcolor,annotation_row = anno_col,
                            annotation_col =  anno_row,scale="row",gaps_row = gaps_row,gaps_col=gaps_col,
                            breaks=bk,
                            legend = legend,border_color =border_color,border = border,cutree_rows = cutree_rows,
                            cutree_cols = cutree_cols,fontsize_row=fontsize_row,fontsize_col=fontsize_col,
                            fontsize=fontsize)%>%ggplotify::as.ggplot()
  # }else{
  #   heat=pheatmap::pheatmap(as.matrix(df),cluster_row = cluster_row,cluster_cols = cluster_col,annotation_legend=annotation_legend,
  #                           cellwidth = cellwidth, cellheight = cellheight,gaps_row = gaps_row,
  #                           gaps_col=gaps_col,labels_row = labels_row,labels_col = labels_col,scale="row",
  #                           breaks=bk,cutree_rows = cutree_rows,
  #                           cutree_cols = cutree_cols,
  #                           main = title,show_rownames=show_rownames,show_colnames = show_colnames,
  #                           color = cellcolor,annotation_row = anno_col,annotation_col = anno_row,
  #                           legend = legend,border_color =border_color,border = border,fontsize_row=fontsize_row,fontsize_col=fontsize_col,
  #                           fontsize=fontsize)%>%ggplotify::as.ggplot()
  # }


  return(heat)
}

