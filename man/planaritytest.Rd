\name{planaritytest}
\alias{planaritytest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boyer-Myvold Planarity test of a network
}
\description{
wrapper function of _MEGENA_planaritytest. imports from Boost graph library, and test planarity of a network
}
\usage{
planaritytest(N, rows, cols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{must be an integer. number of nodes in the network.}
  \item{rows}{first column of edgelist. a vector of integers.}
  \item{cols}{second column of edgelist. a vector of integers.}
}
\details{
cbind(rows,cols) is equivalent to the two column edge list of the network. We assume that the network is undirected. 
}
\value{
 TRUE/FALSE is returned to indicate planarity. (TRUE -> network is planar). 
}
\author{
Won-Min Song
}
\examples{
# test simplest case of planar network (a 3-clique).
planaritytest(as.integer(3),c(1,1,2),c(2,3,3))
}