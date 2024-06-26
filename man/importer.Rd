% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importer_group.R
\name{importer}
\alias{importer}
\title{importer()}
\usage{
importer(path, header = TRUE, sep = ",")
}
\arguments{
\item{path}{Path of file loaded. The file should be in '.csv' format}

\item{header}{Logical. Whether to use the first row as header.}

\item{sep}{Character. The seperator of the file.}
}
\value{
A dataframe, with the first row set as header and the first column set as row name.
The data is unified to numeric.
}
\description{
Internal function to import data in lipidomicR, in order to unify data format.
}
