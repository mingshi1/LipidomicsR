#' @title plotRadar()
#' @description A function to produce radar plot based on lipid types, carbon number, and number of double bonds.
#' @param data A dataframe storing absolute concentration or PL% of lipids between different samples. If not, use cleanXpert() or noridx.output() + normalization.calculator()
#' to get normalized data.
#' The column name should be the sample name and the row name should be the lipid type.
#' The class of column name and row name should be "character". The class of values should be
#' "numeric". The row names are recommended to be in a form like "PL(14:0/20:1)" or "LPL(16:1)".
#'
#' @param pattern Can accept 4 values: "lipid", "CB", "sat", or "all"
#'
#'     If pattern="lipid“, a new radar diagram based on lipid type will be saved, which was named as "lipid_RadarChart.pdf"
#'
#'     If pattern="CB“, a new radar diagram based on carbon number will be saved, which was named as "carbon_number_RadarChart.pdf"
#'
#'     If pattern="sat“, a new radar diagram based on the number of double bonds will be saved, which was named as "unsaturation_RadarChart.pdf"
#'
#'     If pattern="all“, all three diagrams will be saved.
#' @param group A vector defining which group the replicates belong to. Notice: the number of groups should be less than 17.
#' @param max The maximal absolute concentration or PL\% values of each class. The default value is 0.6.
#' @param min The minimal absolute concentration or PL\% values of each class. The default value is 0.
#' @param method The method to select the representative value from a group, which can be "median" or "mean".
#' If it equals "median", the median of the group samples will be chosen. Otherwise, the mean will be chosen to plot.
#' @param axislabcol The color of axis, default value is "grey".
#' @param plwd Defines the width of the data series line. Default value is 2.
#' @param plty Specifies the style of the data series line, which can be 1-6. Default value is 1.
#' @param cglcol Specifies the color of the gridlines. Default value is 1.
#' @param cglwd  Specifies the width of the gridlines. Default value is 1.
#' @param seg Defines the number of gridlines. Default value is 4,
#' which means 5 gridlines: "0\%", "25\%", "50\%", "75\%", and "100\%".
#' @param cglty  Specifies the grid line style, which can be 1-6. Default value is 3.
#' @param vlcex Specifies the size of the group label font. Default value is 1.
#' @param axistype Specifies the style of the axis, which can be 0-5. Default value is 1.
#' @param t.size The size of picture title. Default value is 15.
#' @param t.vjust The vertical position of picture title, which can be negative or positivew values. Default value is 0.
#' @param t.color The color of picture title. Default value is "black".
#' @param l.postion The position of legend, which can be "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" or "center". Default value is "topright".
#' @param l.bty Whether the legend box is drawn, "o" means drawn, and the default value is "n" not drawn.
#' @param lt.col The color of legend text. Default value is "grey25".
#' @param lt.cex The fontsize of legend text. Default value is 2.
#' @param l.cex The size of legend. Default value is 1.
#' @import dplyr
#' @import tidyselect
#' @return No return value, called for side effects, which is a radar diagram based on lipid type, carbon number, or unsaturation among different groups.
#' @examples
#' WT_1=rnorm(n=10,mean=0.4,sd=0.1)
#' WT_2=rnorm(n=10,mean=0.4,sd=0.1)
#' WT_3=rnorm(n=10,mean=0.4,sd=0.1)
#' WT_4=rnorm(n=10,mean=0.4,sd=0.1)
#' KO_1=rnorm(n=10,mean=0.8,sd=0.1)
#' KO_2=rnorm(n=10,mean=0.8,sd=0.1)
#' KO_3=rnorm(n=10,mean=0.8,sd=0.1)
#' KO_4=rnorm(n=10,mean=0.8,sd=0.1)
#' WT_treat_1=rnorm(n=10,mean=0.1,sd=0.1)
#' WT_treat_2=rnorm(n=10,mean=0.1,sd=0.1)
#' WT_treat_3=rnorm(n=10,mean=0.1,sd=0.1)
#' WT_treat_4=rnorm(n=10,mean=0.1,sd=0.1)
#' KO_treat_1=rnorm(n=10,mean=0.6,sd=0.1)
#' KO_treat_2=rnorm(n=10,mean=0.6,sd=0.1)
#' KO_treat_3=rnorm(n=10,mean=0.6,sd=0.1)
#' KO_treat_4=rnorm(n=10,mean=0.6,sd=0.1)
#' data=data.frame(WT_1,WT_2,WT_3,WT_4,KO_1,KO_2,KO_3,KO_4,
#'                 WT_treat_1,WT_treat_2,WT_treat_3,WT_treat_4,
#'                 KO_treat_1,KO_treat_2,KO_treat_3,KO_treat_4)
#' rownames(data)=c("LPC(16:0)","PC(14:0/16:1)","PC(18:1/18:1)","PE(18:0/20:1)",
#'                  "PS(20:1/20:1)","PI(16:0/16:1)","PC(18:0/18:1)","PA(16:0/16:1)",
#'                  "LPE(18:0)","PE(O-18:1/18:0)")
#' group=rep(c("WT","KO","WT_treat","KO_treat"),each=4)
#' plotRadar(data,"all",group) # This is the most simplified version
#' @export
plotRadar<-function(data,pattern,group,max=0.6,min=0,method="median",
                    axislabcol="grey",plwd=2,plty=1,
                    cglcol=1,seg=4,cglwd=1,cglty=3,vlcex=1,axistype = 1,
                    t.size=15,t.vjust=0,t.color="black",
                    l.postion="topright",l.bty = "n",
                    lt.col = "grey25", lt.cex = 2,l.cex=1){
  pat<-function(){
    if(pattern=="lipid"){
      df=sepclass(data,pattern)%>%pivot_longer(-.data$`lipid_type`,names_to = "species",values_to = "value")
      df$group=rep(group,time=length(data[,1]))
      if (method=="median"){
        dt1=df%>%group_by(.data$group,.data$lipid_type)%>%summarise(median=median(sum(.data$value)))
        dt=data.frame(lipid_type=unique(dt1$lipid_type))
        for (s in unique(group)) {
          dt=cbind(dt,dt1$median[dt1$group==s])
        }
      }else if(method=="mean"){
        dt1=df%>%group_by(.data$group,.data$lipid_type)%>%summarise(mean=mean(sum(.data$value)))
        dt=data.frame(lipid_type=unique(dt1$lipid_type))
        for (s in unique(group)) {
          dt=cbind(dt,dt1$mean[dt1$group==s])
        }
      }else{
        stop("Please input right 'method' value, which should be 'median' or 'mean'.")
      }
      rownames(dt)=dt$lipid_type
      dt=dt[,-1]
      colnames(dt)=unique(group)
      dt=apply(dt,2,function(x) apply(as.data.frame(x),1,function(y) y/sum(x)))
      range=data.frame(max=rep(max,dim(dt)[1]),min=rep(min,dim(dt)[1]))%>%cbind(dt)
      mydt=as.data.frame(t(range))
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector = c(RColorBrewer::brewer.pal(9, "Set1"),RColorBrewer::brewer.pal(8, "Dark2"))
      col=col_vector[1:(dim(mydt)[1]-2)]
      #      pdf(paste0(pattern,"_radarchart.pdf"))
      fmsb::radarchart(mydt,axislabcol=axislabcol,plwd=plwd,plty=plty,
                       cglcol=cglcol,seg=seg,cglwd=cglwd,cglty=cglty,vlcex=vlcex,
                       axistype = axistype,pty = 32, na.itp = FALSE,
                       pcol = scales::alpha(col, 2), pfcol =scales::alpha(col, 0.2))
      legend(l.postion,
             legend = rownames(mydt)[-c(1,2)],
             bty = l.bty, pch = 20, col = scales::alpha(col, 2),
             text.col = lt.col, pt.cex = lt.cex,cex=l.cex)
      title(main = "Radar plot based on Lipid Type",cex=t.size, line=t.vjust,col=t.color)
      #dev.off()
      #cat(paste0("\nNotice:\nPictures have been saved at path '",getwd(),"/lipid_RadarChart.pdf'.\n\n"))
    }else if(pattern=="CB"){
      df=sepclass(data,pattern)%>%pivot_longer(-.data$`carbon_number`,names_to = "species",values_to = "value")
      df$group=rep(group,time=length(data[,1]))
      if (method=="median"){
        dt1=df%>%group_by(.data$group,.data$carbon_number)%>%summarise(median=median(sum(.data$value)))

      }else if(method=="mean"){
        dt1=df%>%group_by(.data$group,.data$carbon_number)%>%summarise(median=mean(sum(.data$value)))
      }else{
        stop("Please input right 'method' value, which should be 'median' or 'mean'.")
      }
      dt=data.frame(carbon_number=unique(dt1$carbon_number))
      for (s in unique(group)) {
        dt=cbind(dt,dt1$median[dt1$group==s])
      }
      rownames(dt)=dt$carbon_number
      dt=dt[,-1]
      colnames(dt)=unique(group)
      dt=apply(dt,2,function(x) apply(as.data.frame(x),1,function(y) y/sum(x)))
      range=data.frame(max=rep(max,dim(dt)[1]),min=rep(min,dim(dt)[1]))%>%cbind(dt)
      mydt=as.data.frame(t(range))
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector = RColorBrewer::brewer.pal(9, "Set1")
      col=col_vector[1:(dim(mydt)[1]-2)]
      #      pdf(paste0(pattern,"_radarchart.pdf"))
      fmsb::radarchart(mydt,axislabcol=axislabcol,plwd=plwd,plty=plty,
                       cglcol=cglcol,seg=seg,cglwd=cglwd,cglty=cglty,vlcex=vlcex,
                       axistype = axistype,pty = 32, na.itp = FALSE,
                       pcol = scales::alpha(col, 2), pfcol =scales::alpha(col, 0.2))
      legend(l.postion,
             legend = rownames(mydt)[-c(1,2)],
             bty = l.bty, pch = 20, col = scales::alpha(col, 2),
             text.col = lt.col, pt.cex = lt.cex,cex=l.cex)
      title(main = "Radar plot based on carbon number",cex=t.size, line=t.vjust,col=t.color)
      # dev.off()
      # cat(paste0("\nNotice:\nPictures have been saved at path '",getwd(),"/carbon_number_RadarChart.pdf'.\n\n"))
    }else if(pattern=="sat"){
      df=sepclass(data,pattern)%>%pivot_longer(-.data$unsaturation,names_to = "species",values_to = "value")
      df$group=rep(group,time=length(data[,1]))
      if (method=="median"){
        dt1=df%>%group_by(.data$group,.data$unsaturation)%>%summarise(median=median(sum(.data$value)))
      }else if(method=="mean"){
        dt1=df%>%group_by(.data$grpu,.data$unsaturation)%>%summarise(median=mean(sum(.data$value)))
      }else{
        stop("Please input right method value, which should be 'median' or 'mean'.")
      }
      dt=data.frame(unsaturation=unique(dt1$unsaturation))
      for (s in unique(group)) {
        dt=cbind(dt,dt1$median[dt1$group==s])
      }
      rownames(dt)=dt$unsaturation
      dt=dt[,-1]
      colnames(dt)=unique(group)
      dt=apply(dt,2,function(x) apply(as.data.frame(x),1,function(y) y/sum(x)))
      range=data.frame(max=rep(max,dim(dt)[1]),min=rep(min,dim(dt)[1]))%>%cbind(dt)
      mydt=as.data.frame(t(range))
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector = RColorBrewer::brewer.pal(9, "Set1")
      col=col_vector[1:(dim(mydt)[1]-2)]
      #   pdf(paste0(pattern,"_radarchart.pdf"))
      fmsb::radarchart(mydt,axislabcol=axislabcol,plwd=plwd,plty=plty,
                       cglcol=cglcol,seg=seg,cglwd=cglwd,cglty=cglty,vlcex=vlcex,
                       axistype = axistype,pty = 32, na.itp = FALSE,
                       pcol = scales::alpha(col, 2), pfcol =scales::alpha(col, 0.2))
      legend(l.postion,
             legend = rownames(mydt)[-c(1,2)],
             bty = l.bty, pch = 20, col = scales::alpha(col, 2),
             text.col = lt.col, pt.cex = lt.cex,cex=l.cex)
      title(main = "Radar plot basewd on Number of Double Bonds",cex=t.size,line=t.vjust,col=t.color)
      # dev.off()
      # cat(paste0("\nNotice:\nPictures have been saved at path '",getwd(),"/unsaturation_RadarChart.pdf'.\n\n"))
    }else{
      stop("\nThe parameter 'pattern' is wrong, which should be 'lipid', 'CB', or 'sat'.\n\n")
    }
  }
  if(pattern=="all"){
    for (pattern in c("lipid","CB","sat")) {
      pat()
    }
  }else{
    pat()
  }
  return(NULL)
}


