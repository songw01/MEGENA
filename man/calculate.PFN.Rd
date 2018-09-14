\name{calculate.PFN}
\alias{calculate.PFN}
\alias{resume.calculate.PFN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PFN calculation
}
\description{
main function to calculate PFN a ranked list of edge pairs
}
\usage{
calculate.PFN(edgelist,max.skipEdges = NULL,maxENum = NULL,doPar = FALSE,
	num.cores = NULL, keep.track = TRUE, save.vertex.names=TRUE,
	vertex.names.file='vertex.names.RData')
resume.calculate.PFN(edgelist, net, max.skipEdges = NULL, doPar = FALSE,
	num.cores = NULL, keep.track = TRUE, vertex.names.file='vertex.names.RData')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edgelist}{three column edgelist: first two columns are topological edges, and the third column is the weight. Must be a data.frame object.}
  \item{max.skipEdges}{Maximum number of edges to be searched by planarity test without any inclusion to PFN. 
  If set NULL, it will be automatically set to number of cores x 1000. 
  It acts as a threhold to quicken PFN construction termination during PCP.}
  \item{maxENum}{maximum number of edges to include in final PFN. Default value is NULL, 
  which invokes maximal number of edges allowed in planar network.}
  \item{doPar}{TRUE/FALSE logical variable to choose parallelization.}
  \item{num.cores}{number of cores to use in parallelization.}
  \item{keep.track}{If TRUE, pfg_el.RData will be created in working folder. This file can be used later for restart in case PFN construction did not finish successfully. Default is TRUE.}
  \item{net}{a data.frame with intermediate network edges from a PFN job. The edges should have been converted to numeric index of vector \code{vertex.names} given by \code{vertex.names.file}.}
  \item{save.vertex.names}{whether to save vertex names that matches the gene indices in intermediate output.}
  \item{vertex.names.file}{file to save vertex names if \code{save.vertex.names} is \code{TRUE}.}
}
\details{
If \code{doPar = TRUE}, then num.cores are registered for PCP. 
}
\value{
 output is three column edgelist data.frame, third column being the weight.
}
\author{
Won-Min Song
}
\examples{
# test simplest case of planar network (a 3-clique).
a <- c(1,1,2);b <- c(2,3,3);w <- runif(3,0,1);
el <- cbind(a,b,w);el <- as.data.frame(el[order(el[,3],decreasing = TRUE),])
calculate.PFN(edgelist = el,max.skipEdges = Inf,doPar = FALSE,num.cores = NULL)
}