#' @title importer()
#'
#' @description
#' Internal function to import data in lipidomicR, in order to unify data format.
#'
#' @param path Path of file loaded. The file should be in '.csv' format
#' @param header Logical. Whether to use the first row as header.
#' @param sep Character. The seperator of the file.
#'
#' @return A dataframe, with the first row set as header and the first column set as row name.
#' The data is unified to numeric.
#'
#' @export
importer <- function(path, header = TRUE, sep = ',') {
  data<-read.csv(path,header = header,row.names = 1,sep = sep)
  data[is.na(data)] <- 0
  data[data=='N/A' | data == 'na' | data == 'NA'] <- NA
  data2<- as.data.frame(lapply(data, as.numeric), row.names = rownames(data),check.names = FALSE)
  return(data2)
}

#' @title groupXpert()
#'
#' @param data Data.frame. Row data to be grouped.
#' @param sep Character. The separator used to separate variables from the number of repeats
#' @param as.name Logical. If is true, the group name will be used as the name of the return value. If is false, the sample name will be used as the name of the return value.
#' @param specify List. Used to set groups manually. Each sublist has a group name and a value based on the column range in which it is located. No duplicate values can appear in the sublist.
#'
#' @return A vector named as the colname and uses group name as the value.
#' @export
groupXpert <- function(data, sep = '_', as.name = TRUE, specify = NULL) {
  if (TRUE %in% duplicated(unlist(specify))) {
    stop('ERROR: duplicated column selected.')
  }
  if (NA %in% as.numeric(unlist(specify))) {
    stop('ERROR: column must represented in numeric.')
  }

  name <- unlist(lapply(colnames(data),
                  function(x){
                    x1 <- strsplit(x, split = sep)[[1]]
                    if (length(x1) == 1) {
                      return(x1)
                    }else if (length(x1) > 1) {
                      x1 <- x1[-length(x1)]
                      return(paste(x1, collapse = sep))
                    }
                    }))

  for (i in 1:length(specify)) {
    name[specify[[i]]] <- names(specify[i])
  }
  if (as.name == TRUE) {
    group <- colnames(data)
    names(group) <- name
  }else{
    group <- name
    names(group) <- colnames(data)
  }
  return(group)
}
