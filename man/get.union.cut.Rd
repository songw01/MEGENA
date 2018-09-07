\name{get.union.cut}
\alias{get.union.cut}
\title{
Scale-thresholding of multiscale modules.
}
\description{
obtain a discrete, disjoint clustering results from multiscale MEGENA modules for a given alpha value. 
}
\usage{
get.union.cut(module.output,alpha.cut,output.plot = T,
plotfname = "validModules_alpha",module.pval = 0.05,remove.unsig = T)
}

\arguments{
  \item{module.output}{A direct output from "do.MEGENA". (i.e. MEGENA.output$module.output).}
  \item{alpha.cut}{Resolution parameter cut-off (i.e. alpha) value. alpha.cut = 1 corresponds to classical definition of "small-world" compactness.}
  \item{output.plot}{TRUE/FALSE to indicate outputting a .png file showing hierarchical structure with final outputted modules highlighted in red.}
  \item{plotfname}{.png file outputname.}
  \item{module.pval}{module significance p-value.}
  \item{remove.unsig}{TRUE/FALSE indicating to remove insignificant clusters.}
}
\details{
Returns a list object where each entry is a module.
}
\author{
Won-Min Song
}
\examples{
\dontrun{
rm(list = ls())
data(Sample_Expression)
ijw <- calculate.correlation(datExpr[1:100,],doPerm = 2)
el <- calculate.PFN(ijw[,1:3])
g <- graph.data.frame(el,directed = FALSE)
MEGENA.output <- do.MEGENA(g = g,remove.unsig = FALSE,doPar = FALSE,n.perm = 10)
get.union.cut(module.output = MEGENA.output$module.output,alpha.cut = 1,
output.plot = FALSE,plotfname = NULL,module.pval = 0.05,remove.unsig = TRUE)
}
}