\name{compute.PFN.par}
\alias{compute.PFN.par}
\title{
Parallelized PFN computation
}
\description{
PFN construction by parallelized edge screening. 
}
\usage{
compute.PFN.par(sortedEdge,Ng,maxENum,Njob,Ncore,max.skipEdges = NULL,
keep.track = TRUE,initial.links = NULL)
}
\arguments{
  \item{sortedEdge}{3-column matrix for the input edgelist (e.g. - correlation pair list). Must be sorted by third column, which is usually weight vector.}
  \item{Ng}{integer. number of genes included in sortedEdge.}
  \item{maxENum}{Maximum number of edges to include in final PFN. The theoretical maximal number enforced by Euler's formula is 3(Ng-2).}
  \item{max.skipEdges}{Maximum number of edges to be counted before any valid edge to be included in PFN. This works as a termination condition to avoid exhaustive 
planarity testing over all edges provided in sortedEdge.}
  \item{Njob}{Number of edges to be passed to each core for parallelized edge screening.}
  \item{Ncore}{Number of cores to utilize.}
  \item{keep.track}{TRUE/FALSE logical. Indicate if the record of PFN construction is saved in temporary file "pfg_el.RData". Default is TRUE.}
  \item{initial.links}{If provided, PFN construction will restart by regarding these initial.links as already-built PFN.}
}
\details{
This is parallelized implementation of PFN construction, where it is possible to re-capture PFN construction by providing already computed edgelist into initial.links. 
Although provivded, this function itself may require careful caution and users are encouraged to use more user-friendly "calculate.PFN()" instead. 
}
\value{
A 3-column matrices, where first two columns are integer indices for vertices, and third is the weight vector.
}
\author{
Won-Min Song
}