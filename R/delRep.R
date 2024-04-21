#' @title delRep()
#' @description
#' A function to delete specific number of replicates, the replicates causing largest
#' standard deviation will be deleted.
#' @param data A dataframe storing concentration of lipids between different samples.
#' The column name should be the sample name and the row name should be the lipid type.
#' The class of column name and row name should be "character". The class of values should be
#' "numeric". The row names are recommended to be in a form like "PL(14:0/20:1)" or "LPL(16:1)".
#' @param group Vector. The group information, recommended to be generated with groupXpert()
#' @param m The number of replicates you want to delete.
#' @param method The method to find the worst replicates, can be "PCA" or "Euclidean". Default value is "PCA".
#' @param show.del Whether to show deleted columns. Default value is FALSE.
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
#' m=1
#' group=colnames(data)
#' names(group)=rep(c("WT","KO"),each=4)
#' delRep(data,group,m)
#' @export
delRep<-function(data,group,m,method="PCA",show.del=FALSE){
  spdt=data.frame(species=rownames(data))
  rownames(spdt)=spdt$species
  spdt=spdt[,-1]
  del=data.frame(species=rownames(data))
  rownames(del)=del$species
  del=del[,-1]
  for (i in unique(names(group))){
    name <- group[names(group) == i]
    n <- length(name)
    dt <- data[,name]
    if (method=="PCA"){
      pdt=as.data.frame(prcomp(dt)$rotation)
      c1=mean(pdt[,1])
      c2=mean(pdt[,2])
      pdt$res=0
      for (a in 1:n){
        pdt$res[a]=(pdt[a,1]-c1)^2+(pdt[a,2]-c2)^2

      }
      avt=order(pdt$res,decreasing = TRUE)
      res=dt[,-avt[1:m]]
      colnames(res)=name[1:(n-m)]
      dc=as.data.frame(dt[,avt[1:m]])
      colnames(dc)=name[avt[1:m]]
    }else{
      ddt=as.matrix(dist(t(dt)))
      avg_distances <- as.data.frame(colMeans(ddt))
      avt=order(avg_distances$`colMeans(ddt)`,decreasing = TRUE)
      res=dt[,-avt[1:m]]
      colnames(res)=name[1:(n-m)]
      dc=as.data.frame(dt[,avt[1:m]])
      colnames(dc)=name[avt[1:m]]
    }
    spdt=cbind(spdt,res)
    del=cbind(del,dc)
  }
  if (show.del){
    comb=list(results=spdt,deleted_data=del)
    return(comb)
  }else{
    return(spdt)
  }
}

