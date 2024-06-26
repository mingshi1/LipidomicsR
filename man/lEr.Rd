% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_cleaner.R
\name{lEr}
\alias{lEr}
\title{lEr}
\usage{
lEr(
  data,
  Inlabel = c("PE(17:0/17:0)", "PC(17:0/17:0)", "PS(14:0/14:0)"),
  relative_as_default = TRUE,
  relative.mannual = NULL
)
}
\arguments{
\item{data}{Data frame. Row lipidomics data. Should be imported with importer().}

\item{Inlabel}{Vector. Name of Internal label.}

\item{relative_as_default}{Logical, default as TRUE for automatically searching for internal standard data.}

\item{relative.mannual}{Vector, the exact channel name.}
}
\value{
Data of internal label.
}
\description{
Internal function to extract label data in lipidomicR
}
