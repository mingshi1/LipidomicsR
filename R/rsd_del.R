#' Calculate Relative Standard Deviation (RSD)
#'
#' This function calculates the relative standard deviation (RSD) based on the specified data range.
#'
#' @param data The data frame containing abundance data.
#' @param start The starting column index of the QC data range.
#' @param end The ending column index of the QC data range.
#' @param threshold The threshold value for RSD. Default is 0.2.
#' @param show.del Logical value indicating whether to show the deleted data. Default is FALSE.
#' @param del.zero Logical value indicating whether to delete rows with all QC being zero . Default is TRUE.
#' @return A data frame containing the calculated RSD values and the corresponding data.
#' @examples
#' qc_1=rnorm(n=5,mean=0.3,sd=0.2)
#' qc_2=rnorm(n=5,mean=0.3,sd=0.2)
#' qc_3=rnorm(n=5,mean=0.3,sd=0.2)
#' qc_4=rnorm(n=5,mean=0.3,sd=0.2)
#' qc_5=rnorm(n=5,mean=0.3,sd=0.2)
#' WT_1=rnorm(n=5,mean=0.3,sd=0.1)
#' WT_2=rnorm(n=5,mean=0.3,sd=0.1)
#' WT_3=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_1=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_2=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_3=rnorm(n=5,mean=0.3,sd=0.1)
#' data=data.frame(qc_1,qc_2,qc_3,qc_4,qc_5,WT_1,WT_2,WT_3,KO_1,KO_2,KO_3)
#' rownames(data)=c("LPC(16:0)","PC(14:0/16:1)","PC(18:1/18:1)","PE(18:0/20:1)","PS(20:1/20:1)")
#' rsd_calculator(data,1,5,show.del = TRUE)
#' @export

rsd_calculator<-function(data,start,end,threshold=0.2,show.del=FALSE,del.zero=TRUE){
  data=na.omit(data)
  data[data=='N/A' | data == 'na' | data == 'NA'] <- NA
  data[is.na(data)] <- 0
  cal_dt=data[,c(start:end)]
  rsd=c()
  for (i in 1:dim(data)[1]){
    dt=as.numeric(cal_dt[i,])
    rsd=c(rsd,sd(dt)/mean(dt))
  }
  cal_dt$rsd=rsd
  if (del.zero){

    data=data[!is.na(cal_dt$rsd),]
    cal_dt=cal_dt[!is.na(cal_dt$rsd),]
  }else{
    cal_dt[is.na(cal_dt)] <- 0
  }
  del=cbind(cal_dt$rsd[cal_dt$rsd>threshold],data[cal_dt$rsd>threshold,])
  colnames(del)[1]="RSD"
  res=cbind(cal_dt$rsd[cal_dt$rsd<=threshold],data[cal_dt$rsd<=threshold,])
  colnames(res)[1]="RSD"
  if (show.del){
    return(list(results=res,deleted_data=del))
  }else{
    return(del)
  }
}


