#' @title lEr
#' @description
#' Internal function to extract label data in lipidomicR
#'
#' @param data Data frame. Row lipidomics data. Should be imported with importer().
#'
#' @param Inlabel Vector. Name of Internal label.
#'
#' @param relative_as_default Logical, default as TRUE for automatically searching for internal standard data.
#'
#' @param relative.mannual Vector, the exact channel name.
#'
#' @return Data of internal label.
lEr  <- function(data,Inlabel=c('PE(17:0/17:0)','PC(17:0/17:0)','PS(14:0/14:0)'),
                 relative_as_default = TRUE, relative.mannual= NULL) {
  if (relative_as_default == TRUE) {
    data_cleaned <- channel.delete(data)
    lb<-data_cleaned[Inlabel,]
    lb<-as.data.frame(lapply(lb,as.numeric),row.names = rownames(lb),check.names = FALSE)
  }else {
    lb <- data[relative.mannual,]
  }
  return(lb)
}


#' @title lEa()
#'
#' @description
#' Internal function to extract label data in lipidomicR
#'
#' @param data Data frame. Row lipidomics data. Should be imported with importer().
#' @param Inchannel Vector. Exact channel name for isotope label data.
#' @param sample Vector. Column of sample added with isotope label. Default as 2:5.
#'
#' @return Data of internal label.
lEa<- function(data,Inchannel = c('PC(15:0/18:1(d7))+AcO_1',
                                  'PE(15:0/18:1(d7))_2',
                                  'PS(15:0/18:1(d7))_2',
                                  'PG(15:0/18:1(d7))_2',
                                  'PI(15:0/18:1(d7))_2',
                                  'PA(15:0/18:1(d7))_2',
                                  'LPC18:1(d7)+AcO',
                                  'LPE18:1(d7)'), sample = 2:5){
  lb<-data[Inchannel,sample]
  lb<-as.data.frame(lapply(lb,as.numeric),row.names = rownames(lb),check.names = FALSE)
  return(lb)
}


#' @title nor.absolute()
#' @description
#' A function to calculate parameters for absolute normalization.
#' @param data Data frame, row lipidomics data.
#' @param radio.data Data frame. Characteristic of the radio label data.
#' The row name must be the exact channel name of the label. Molecular mass should be
#' provided in a column named "Mass". Concentration should be provided in a column named "Concentration(mg/ml)".
#' @param sample Vector. Column of sample added with isotope label. Default as 2:5.
#'
#' @return Parameters for absolute normalization
nor.absolute  <- function(data, radio.data=NULL, sample = 2:5) {
  rad<-radio.data
  cM<-radio.data$Concentration.mg.ml./radio.data$Mass
  avg_rad<- apply(lEa(data, Inchannel = rownames(radio.data), sample = sample),1,'mean')
  unit<-cM/avg_rad
  rad$Concentration.mol.l <- cM
  rad$mol.peak_intensity <- unit
  return(rad)
}

#' @title nor.relative()
#' @description
#' A function to calculate parameters for relative normalization.
#'
#' @param data Data frame, row lipidomics data.
#' @param Inlabel Vector. Name of Internal label. Default as c('PE(17:0/17:0)','PC(17:0/17:0)','PS(14:0/14:0)')
#' @param normalize_to Vector. The column of samples that used as the standard for normalization.
#' @param relative_as_default Logical, default as TRUE for automatically searching for internal label data.
#' @param relative.mannual Vector, the exact channel name (if you want to define the channel of internal label mannually).
#'
#' @return Parameters for relative normalization
nor.relative <- function(data,
                         Inlabel=c('PE(17:0/17:0)','PC(17:0/17:0)','PS(14:0/14:0)'),
                         normalize_to = 1:5,
                         relative_as_default = TRUE,
                         relative.mannual = NULL) {
  rel <- lEr(data, Inlabel = Inlabel,
             relative_as_default = relative_as_default,relative.mannual = relative.mannual)
  rel_prop <- apply(rel[,normalize_to], 1, 'mean')/sum(apply(rel[,normalize_to], 1, 'mean'))
  idx<-0
  for (i in 1:length(rel_prop)) {
    ind<- rel[i,]*rel_prop[i]
    idx<-idx+ind
  }
  rownames(idx) <- c('normalization_index')
  return(idx)
}

#' @title noridx()
#'
#' @description
#' An integrated function that call nor.relative() and nor.absolute(), for simplifying.
#'
#' @usage
#' noridx(data, radio.data = NULL,
#' normalization.mode='both',
#' sample = 1:5, normalize_to = 2:5,
#' Inlabel=c('PE(17:0/17:0)','PC(17:0/17:0)','PS(14:0/14:0)'),
#' relative_as_default = TRUE,
#' relative.mannual = NULL)
#'
#' @param data Data frame, row lipidomics data.
#' @param normalization.mode Character. "absolute" tp output absolute normalization index
#' 'relative' to output relative normalization index. "both" to output both of them. Default as "both".
#' @param radio.data Data frame. Characteristic of the radio label data.
#' The row name must be the exact channel name of the label. Molecular mass should be
#' provided in a column named "Mass". Concentration should be provided in a column named "Concentration(mg/ml)".
#' @param sample Vector. Column of sample added with isotope label. Default as 2:5.
#' @param normalize_to Vector. The column of samples that used as the standard for normalization.
#' @param Inlabel Inlabel Vector. Name of Internal label. Default as c('PE(17:0/17:0)','PC(17:0/17:0)','PS(14:0/14:0)')
#' @param relative_as_default Logical, default as TRUE for automatically searching for internal label data.
#' @param relative.mannual Vector, the exact channel name (if you want to define the channel of internal label mannually).
#'
#' @return
#' 1. normalization.mode='both'. A list of data frames of normalization indexes of the two modes.
#' 2. normalization.mode='absolute' or 'relative'. A data frame of the respective normlaization index.
#' @export
noridx <- function(data, radio.data = NULL,
                   normalization.mode='both',
                   sample = 1:5, normalize_to = 2:5,
                   Inlabel=c('PE(17:0/17:0)','PC(17:0/17:0)','PS(14:0/14:0)'),
                   relative_as_default = TRUE,
                   relative.mannual = NULL) {
  if (normalization.mode == 'absolute') {
    return(list(absolute=nor.absolute(data = data, radio.data = radio.data, sample = sample)))}
  else if (normalization.mode =='relative'){
    return(list(relative=nor.relative(data = data,Inlabel = Inlabel, normalize_to = normalize_to, relative_as_default = relative_as_default,
                        relative.mannual = relative.mannual)))}
  else if (normalization.mode == "both") {
    x1 <- nor.absolute(data = data, radio.data = radio.data, sample = sample)
    x2 <- nor.relative(data = data,Inlabel = Inlabel, normalize_to = normalize_to, relative_as_default = relative_as_default,
                       relative.mannual = relative.mannual)
    return(list(absolute=x1,relative=x2))
  }else {
    warning('Error: unsupported normalization mode\n')
  }
}

#' @title absolute.calculator()
#'
#' @description
#' Internal function of lipidmicR::normalization_calculator()
#'
#' @param data Row data of lipidomics
#' @param absolute.dataset Data frame. Normalization index.
#'
#' @return
#' Return a data frame of absolute normalized lipidomic data.
absolute.calculator <- function(data,absolute.dataset){
  sep_data<-sepclass(data,'lipid')
  sep_index<-sepclass(absolute.dataset,'lipid')
  cnx<-data.frame()
  for (lipid in sep_index$lipid_type){
    cdt<-subset(data, subset = sep_data$`lipid_type`==lipid)
    cix<-sep_index[sep_index$`lipid_type`==lipid, 'mol.peak_intensity']
    cnx2 <- cdt*cix
    cnx<- rbind(cnx,cnx2)
  }
  return(cnx)
  }

#' @title relative.calculator()
#'
#' @description
#' Internal function of lipidmicR::normalization_calculator()
#'
#' @param data Row data of lipidomics
#' @param relative.dataset Data frame. Normalization index.
#'
#' @return
#' Return a data frame of relative normalized lipidomic data.
relative_calculator <- function(data, relative.dataset) {
  cno<-data.frame(matrix(nrow = nrow(data),ncol = ncol(data)))
  for (i in 1:ncol(relative.dataset)) {
    cdt<-data[,i]
    cix<- relative.dataset[,i]
    cno2<-cdt/cix
    cno[,i]<-cno2
  }
  rownames(cno)<-rownames(data)
  colnames(cno)<-colnames(data)
  return(cno)
  }

#' @title percent.calculator()
#' @description
#' A function to calculate the proportion of each lipid content.
#'
#' @param data Data frame. The row lipidomic data.
#' @param delete.pattern Pattern of characters that needs to be removed.
#'
#' @return
#' Return a data frame of normalized lipidomic data in the percentage of lipid content.
percent.calculator <- function(data, delete.pattern = c('_\\d','(\\+)AcO','_n', '\\(\\d+\\)')) {
  data <- channel.delete(data, delete.pattern = delete.pattern)
  sum <- colSums(data)
  for (i in 1:ncol(data)) {
    data[,i] <- data[,i]/sum[i]
  }
  return(data)
}

#' @title toGroup.calculator()
#'
#' @param data Data frame. The row lipidomic data.
#' @param group Vector. The group information, recommended to be generated with groupXpert()
#' @param to Character. The group to be normalized to.
#' @param delete.pattern Pattern of characters that needs to be removed.
#'
#' @return
#' A dataframe normalized to specified group.
toGroup.calculator <- function(data, group, to,
                               delete.pattern = c('_\\d','(\\+)AcO','_n', '\\(\\d+\\)')) {
  data <- channel.delete(data, delete.pattern = delete.pattern)
  each_avg <- rowMeans(data[group[names(group) == to]])
  for (lipid in names(each_avg)) {
    data[lipid,] <- data[lipid,]/each_avg[lipid]
  }
  return(data)
}



#' @title normalization_calculator()
#'
#' @description
#' A function uses normalization parameters to calculate normalized lipidomic data.
#'
#' @param data Data frame. Row lipidomic data.
#' @param normalization.mode Vector. Options can be 'absolute', 'relative', 'TSA', 'toGroup'.
#' @param normalization.index Data frame. normalization parameters, Suggested to be calculated by noridx.ouput()
#' @param delete.pattern Pattern of characters that needs to be removed.
#' @param group Vector. The group information, recommended to be generated with groupXpert()
#' @param to Character. The group to be normalized to.
#'
#' @return
#' A data frame of normalized data if either 'relative' or 'absolute' mode is used. A list if both of them are used.
#'
#' @export
normalization_calculator <- function(data,
                                     normalization.mode,
                                     normalization.index,
                                     delete.pattern = c('_\\d','(\\+)AcO','_n', '\\(\\d+\\)'),
                                     group,
                                     to
                                     ){
  ncl <- list()
  if ('absolute' %in% normalization.mode) {
    outputa <- absolute.calculator(data = data, absolute.dataset = normalization.index$absolute)
    ncl$absolute <- outputa
  }
  if ('relative' %in% normalization.mode) {
    outputr <- relative_calculator(data = data,relative.dataset = normalization.index$relative)
    ncl$relative <- outputr
  }
  if ('TSA' %in% normalization.mode) {
    outputp <- percent.calculator(data = data, delete.pattern = delete.pattern)
    ncl$percent <- outputp
  }
  if ('toGroup' %in% normalization.mode) {
    outputp <- toGroup.calculator(data = data, group = group, to = to)
    ncl$toGroup <- outputp
  }
  ncl <- lapply(ncl, function(x) channel.delete(x, delete.pattern = delete.pattern))
  return(ncl)
}
