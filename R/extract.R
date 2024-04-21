#' @title sepclass()
#' @description
#' A function to identify lipid type and calculate the number of carbons and
#' unsaturation rate of lipids.
#' @param data  A dataframe storing concentration of lipids between different samples.
#' The column name should be the sample name and the row name should be the lipid type.
#' The class of column name and row name should be "character". The class of values should be
#' "numeric". The row names are recommended to be in a form like "PL(14:0/20:1)" or "LPL(16:1)".
#' @param pattern   Can accept 4 values: "lipid", "CB", "sat", or "all"
#'
#'     If pattern="lipid“, a new column named "lipid_type" will be added, which stores type of lipid like "PE", "LPC", and "TAG".
#'
#'     If pattern="CB“, a new column named "carbon_number" will be added, which stores the number of carbons. For example, the carbon_number of "PC(14:0/16:1)" is 14+16=30.
#'
#'     If pattern="sat“, a new column named "unsaturation" will be added, which stores the number of double bonds. For example, the unsaturation of "PC(14:1/16:1)" is 1+1=2.
#'
#'     If pattern="all“, all three columns will be added.
#' @param subtype A logic value to determine for a lipid like "PC(O-14:0/16:1)", "lipid_type" should be
#' "PC" (subtype=FALSE) or "PC(O)" (subtype=TRUE). Default value is FALSE.
#' @return A dataframe with new columns containing the class of lipid type, carbon number,
#' or unsaturation based on the original data input.
#' @examples
#' WT_1=rnorm(n=5,mean=0.3,sd=0.1)
#' WT_2=rnorm(n=5,mean=0.3,sd=0.1)
#' WT_3=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_1=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_2=rnorm(n=5,mean=0.3,sd=0.1)
#' KO_3=rnorm(n=5,mean=0.3,sd=0.1)
#' data=data.frame(WT_1,WT_2,WT_3,KO_1,KO_2,KO_3)
#' rownames(data)=c("LPC(16:0)","PC(O-14:0/16:1)","TAG56:2-FA20:1","PE(P-18:0/20:1)","PS(20:1/20:1)")
#' pattern="all" ## or "lipid", "CB", "sat"
#' sepclass(data,pattern)
#' @import stringr
#' @export
sepclass<-function(data,pattern,subtype=FALSE){
  l=gsub("\\-[A-Za-z].+","",rownames(data))
  dt=data.frame(x=rep(0,length(data[,1])))
  sub=str_extract(rownames(data), "\\((.*?)\\-")
  sub=gsub("\\-", "\\)", sub)
  sub[is.na(sub)]=""
  if(pattern=="lipid"){

    dt$lipid_type=gsub("\\(.+","",rownames(data))
    dt$lipid_type=gsub("[0-9].+","",dt$lipid_type)
    dt$lipid_type=gsub(" ","",dt$lipid_type)
    dt=as.data.frame(dt[,-1])
    colnames(dt)="lipid_type"
    if (subtype){
      dt$lipid_type=paste0(dt$lipid_type,sub)
    }
  }else if(pattern=="CB"){
    l1=stringr::str_extract_all(l,pattern = "[0-9]+\\:")
    cb=unlist(lapply(l1,function(x)
      sum(unlist(lapply(x,function(y) as.numeric(gsub(":","",y)))))))
    dt$carbon_number=cb
    dt=as.data.frame(dt[,-1])
    colnames(dt)="carbon_number"
  }else if(pattern=="sat"){
    l2=stringr::str_extract_all(l,pattern = ":[0-9]")
    sat=unlist(lapply(l2,function(x)
      sum(unlist(lapply(x,function(y) as.numeric(gsub(":","",y)))))))
    dt$unsaturation=sat
    dt=as.data.frame(dt[,-1])
    colnames(dt)="unsaturation"
  }else if(pattern=="all"){
    dt$lipid_type=gsub("\\(.+","",rownames(data))
    dt$lipid_type=gsub("[0-9].+","",dt$lipid_type)
    dt$lipid_type=gsub(" ","",dt$lipid_type)
    l1=stringr::str_extract_all(l,pattern = "[0-9]+\\:")
    cb=unlist(lapply(l1,function(x)
      sum(unlist(lapply(x,function(y) as.numeric(gsub(":","",y)))))))
    dt$carbon_number=cb
    l2=stringr::str_extract_all(l,pattern = ":[0-9]")
    sat=unlist(lapply(l2,function(x)
      sum(unlist(lapply(x,function(y) as.numeric(gsub(":","",y)))))))
    dt$unsaturation=sat
    dt=as.data.frame(dt[,-1])
    if (subtype){
      dt$lipid_type=paste0(dt$lipid_type,sub)
    }
  }else{
    stop("ERROR: The parameter 'pattern' is wrong, which should be 'lipid', 'CB', or 'sat'.")

  }
  data=cbind(dt,data)
  return(data)
}


