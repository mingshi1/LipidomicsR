#' @title delRep()
#' @description
#' A function to delete specific number of replicates, the replicates causing largest
#' standard deviation will be deleted.
#' @param data A dataframe storing concentration of lipids between different samples.
#' The column name should be the sample name and the row name should be the lipid type.
#' The class of column name and row name should be "character". The class of values should be
#' "numeric". The row names are recommended to be in a form like "PL(14:0/20:1)" or "LPL(16:1)".
#' @param n The whole number of replicates in one group.
#' @param m The number of replicates you want to delete.
#'
#' @return A new dataframe deleting replicates which cause the largest SD.
#' @examples
#' WT_1=rnorm(n=5,mean=0.3,sd=0.1)
#' WT_2=rnorm(n=5,mean=0.3,sd=0.1)
#' WT_3=rnorm(n=5,mean=0.3,sd=0.1)
#' WT_4=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_1=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_2=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_3=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_4=rnorm(n=5,mean=0.3,sd=0.1)
#' data=data.frame(WT_1,WT_2,WT_3,WT_4,KO_1,KO_2,KO_3,KO_4)
#' rownames(data)=c("LPC(16:0)","PC(14:0/16:1)","PC(18:1/18:1)","PE(18:0/20:1)","PS(20:1/20:1)")
#' n=4
#' m=1
#' delRep(data,n,m)
#' @export
delRep<-function(data,n,m){
  vec=c(1:n)
  cob=t(combn(vec, (n-m),simplify = T))
  spdt=data.frame(species=rownames(data))
  rownames(spdt)=spdt$species
  spdt=spdt[,-1]
  for (i in 1:(ncol(data)/n) ){
    dt=data[,(n*i+1-n):(n*i)]
    df=t(apply(dt,1,function(x) apply(cob,1,function(y) sd(as.data.frame(x)[as.vector(y),]))))
    cm=t(apply(df, 1, function(x) t(as.data.frame(cob[which.min(x),]))))
    sample=data.frame()

    for (i in 1:length(rownames(cm))) {
      sp=dt[i,as.vector(cm[i,])]
      colnames(sp)=colnames(dt)[1:(n-m)]
      sample=rbind(sample,sp)
    }
    rownames(sample)=rownames(cm)
    spdt=cbind(spdt,sample)
  }
  return(spdt)
  }
