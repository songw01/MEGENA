\name{output.geneSet.file}
\alias{output.geneSet.file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
output gene signatures into .gmt file format
}
\description{
An interface function to output .gmt format gene signature file.
}
\usage{
output.geneSet.file(geneSet,outputfname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneSet}{a list object}
  \item{outputfname}{output file name}
}
\details{
Outputs each signature into a single line of lists in \code{outputfname}. 
}
\author{
Won-Min Song
}