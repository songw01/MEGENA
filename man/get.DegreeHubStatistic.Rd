\name{get.DegreeHubStatistic}
\alias{get.DegreeHubStatistic}
\title{
calculate module degree statistics based on random triangulation model via T1 and T2 moves. 
}
\description{
calculation of module p-values. 
}
\usage{
get.DegreeHubStatistic(subnetwork,n.perm = 100,doPar = FALSE,n.core = 4)
}
\arguments{
  \item{subnetwork}{a planar network as an igraph object.}
  \item{n.perm}{number of random networks generated, constraint with number of links and nodes same to "subnetwork".}
  \item{doPar}{TRUE/FALSE to parallelize.}
  \item{n.core}{number of cores/threads to use.}
}
\details{
Hub significance calculation functionality. Make sure that, if doPar = TRUE, register cores using registerDoParallel() from doParallel package.
}
\value{
a data.frame table showing node-wise statistics. 
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

	out <- get.DegreeHubStatistic(subnetwork = g,n.perm = 100,doPar = FALSE,n.core = 4)
}
}