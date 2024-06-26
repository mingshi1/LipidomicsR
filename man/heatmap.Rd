% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.R
\name{heatmap}
\alias{heatmap}
\title{heatmap()}
\usage{
heatmap(
  data,
  group,
  cluster_row = TRUE,
  cluster_col = TRUE,
  sel.group = "default",
  constract = 8.5,
  type = "lipid",
  sel.type = "default",
  sel.row = c("default"),
  annotation_legend = TRUE,
  cellwidth = 20,
  cellheight = 15,
  gaps_row = c(0),
  gaps_col = c(0),
  subtype = FALSE,
  labels_row = c("default"),
  labels_col = c("default"),
  title = "",
  show_rownames = TRUE,
  show_colnames = TRUE,
  cellcolor = c("blue", "black", "yellow"),
  legend = TRUE,
  border_color = NA,
  border = FALSE,
  cutree_rows = 1,
  cutree_cols = 1,
  rtitle = "group",
  ctitle = " ",
  fontsize_row = 12,
  fontsize_col = 12,
  fontsize = 8
)
}
\arguments{
\item{data}{A dataframe storing absolute concentration or PL% of lipids between different samples. If not, use cleanXpert() or noridx.output() + normalization.calculator()
to get normalized data.
The column name should be the sample name and the row name should be the lipid type.
The class of column name and row name should be "character". The class of values should be
"numeric". The row names are recommended to be in a form like "PL(14:0/20:1)" or "LPL(16:1)".}

\item{group}{A vector defining which group the replicates belong to.}

\item{cluster_row}{A boolean variable controlling whether to perfrom clustering to
row variables (lipid abundance) or not. The default value is TRUE.}

\item{cluster_col}{A boolean variable controlling whether to perfrom clustering to
column variables (lipid abundance) or not. The default value is TRUE.}

\item{sel.group}{A vector containing the group you want to show in the heatmap.
The input can be like c("WT","KO"). Default value is "default".}

\item{constract}{The constract of heatmap, default is 8.5, value range from 0 to 10.}

\item{type}{Can accept 3 values: "lipid", "CB", or "sat". Default value is "lipid".

    If type="lipid“, the heatmap will divide rownames based on lipid types.

    If pattern="CB“, the heatmap will divide rownames based on carbon number.

    If pattern="sat“, the heatmap will divide rownames based on the number of double bonds of a lipid type.}

\item{sel.type}{A vector controlling which types to show. If you only want to check data of "PA", "PC“, and "PE".
You can set type="lipid", sel.type=c("PA","PC","PE")}

\item{sel.row}{A vector controlling which types to show. If you set it as c("LPC(16:0)","PC(14:0/16:1)",
"PC(18:1/18:1)","PE(18:0/20:1)"), only their abundance will be shown.}

\item{annotation_legend}{A boolean controlling whether to show the figure legend. The default value is TRUE.}

\item{cellwidth}{The width of a cell in the heatmap. Default value is 20.}

\item{cellheight}{The height of a cell in the heatmap. Default value is 15.}

\item{gaps_row}{To customize positions of row gaps. Default value is c(0).}

\item{gaps_col}{To customize positions of column gaps. Default value is c(0).

Notice: gaps_row and gaps_col are only useful when cluster=FALSE.}

\item{subtype}{A logic value to determine for a lipid like "PC(O-14:0/16:1)", "lipid_type" should be
"PC" (subtype=FALSE) or "PC(O)" (subtype=TRUE). Default value is FALSE.}

\item{labels_row}{A vector contains the labels of each row of the heatmap.
Default value is row names of dataframe input. It can be input like ("PE(20:1/20:1)","PS(16:0/18:1)","","","","LPA(18:0)")}

\item{labels_col}{A vector contains the labels of each column of the heatmap.
Default value is column names of dataframe input.}

\item{title}{The title of heatmap. Default value is "".}

\item{show_rownames}{Whether to show row names or not. Default value is T.}

\item{show_colnames}{Whether to show column names or not. Default value is T.}

\item{cellcolor}{The color range of cells in the heatmap. It should be input in a vector with
three color values, such as c("blue","black","yellow").}

\item{legend}{Whether to show legends or not. Default value is FALSE.}

\item{border_color}{Useful when border=T. Default value is NA.}

\item{border}{Whether to show borders or not. Default value is TRUE.}

\item{cutree_rows}{Useful when cluster=T. If cutree_rows=T, the rows of heatmap will be divided according to
clustering results. Default value is TRUE.}

\item{cutree_cols}{Useful when cluster=T. If cutree_cols=T, the rows of heatmap will be divided according to
clustering results. Default value is TRUE.}

\item{rtitle}{Row title of the heatmap. Default value is "group".}

\item{ctitle}{Column title of the heatmap. Default value is " ".}

\item{fontsize_row}{Fontsize of row labels. Default value is 12.}

\item{fontsize_col}{Fontsize of column labels. Default value is 12.}

\item{fontsize}{Fontsize of all labels. Default value is 8.}
}
\value{
A heatmap that is color-coded by abundance of lipids.
}
\description{
A function to plot a heatmap based on abundance of lipids. A series
of custom functions can be realized such as dividing groups.
}
\examples{
WT_1=rnorm(n=10,mean=0.4,sd=0.1)
WT_2=rnorm(n=10,mean=0.4,sd=0.1)
WT_3=rnorm(n=10,mean=0.4,sd=0.1)
WT_4=rnorm(n=10,mean=0.4,sd=0.1)
KO_1=rnorm(n=10,mean=0.8,sd=0.1)
KO_2=rnorm(n=10,mean=0.8,sd=0.1)
KO_3=rnorm(n=10,mean=0.8,sd=0.1)
KO_4=rnorm(n=10,mean=0.8,sd=0.1)
WT_treat_1=rnorm(n=10,mean=0.1,sd=0.1)
WT_treat_2=rnorm(n=10,mean=0.1,sd=0.1)
WT_treat_3=rnorm(n=10,mean=0.1,sd=0.1)
WT_treat_4=rnorm(n=10,mean=0.1,sd=0.1)
KO_treat_1=rnorm(n=10,mean=0.6,sd=0.1)
KO_treat_2=rnorm(n=10,mean=0.6,sd=0.1)
KO_treat_3=rnorm(n=10,mean=0.6,sd=0.1)
KO_treat_4=rnorm(n=10,mean=0.6,sd=0.1)
data=data.frame(WT_1,WT_2,WT_3,WT_4,KO_1,KO_2,KO_3,KO_4,
               WT_treat_1,WT_treat_2,WT_treat_3,WT_treat_4,
               KO_treat_1,KO_treat_2,KO_treat_3,KO_treat_4)
rownames(data)=c("LPC(16:0)","PC(14:0/16:1)","PC(18:1/18:1)","PE(18:0/20:1)",
                "PS(20:1/20:1)","PI(16:0/16:1)","PC(18:0/18:1)","PA(16:0/16:1)",
                "LPE(18:0)","PE(O-18:1/18:0)")
group=rep(c("WT","KO","WT_treat","KO_treat"),each=4)
heatmap(data,group)
}
