\name{do.MEGENA}
\alias{do.MEGENA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MEGENA clustering + MHA
}
\description{
multiscale clustering analysis (MCA) and multiscale hub analysis (MHA) pipeline
}
\usage{
do.MEGENA(g,
do.hubAnalysis = TRUE,
mod.pval = 0.05,hub.pval = 0.05,remove.unsig = TRUE,
min.size = 10,max.size = 2500,
doPar = FALSE,num.cores = 4,n.perm = 100,singleton.size = 3,
save.output = FALSE)
}

\arguments{
  \item{g}{igraph object of PFN.}
  \item{do.hubAnalysis}{TRUE/FALSE indicating to perform multiscale hub analysis (MHA) in downstream. Default is TRUE.}
  \item{mod.pval}{cluster significance p-value threshold w.r.t random planar networks}
  \item{hub.pval}{hub significance p-value threshold w.r.t random planar networks}
\item{remove.unsig}{TRUE/FALSE indicating to remove insignificant clusters in MHA.}
  \item{min.size}{minimum cluster size}
  \item{max.size}{maximum cluster size}
  \item{doPar}{TRUE/FALSE indicating parallelization usage}
\item{num.cores}{number of cores to use in parallelization.}
\item{n.perm}{number of permutations to calculate hub significance p-values/cluster significance p-values.}
\item{singleton.size}{Minimum module size to regard as non-singleton module. Default is 3.}
\item{save.output}{TRUE/FALSE to save outputs from each step of analysis}
}
\details{
Performs MCA and MHA by taking PFN as input. Returns a list object containing clustering outputs, hub analysis outputs, and node summary table. 
}
\value{
A series of output files are written in wkdir. Major outputs are,
\item{module.output}{outputs from MCA}
\item{hub.output}{outputs from MHA}
\item{node.summary}{node table summarizing clustering results.}
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
}
}