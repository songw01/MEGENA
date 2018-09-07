\name{read.geneSet}
\alias{read.geneSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
.gmt file reader function
}
\description{
An interface function to read-in .gmt format gene signature file.
}
\usage{
read.geneSet(geneSet.file)
}
\arguments{
  \item{geneSet.file}{text file containing gene signatures in .gmt format}
  
}
\details{
Each line of lists in \code{geneset.file} is a single set of signature. 
}
\value{
loads signatures into a list object. 
}
\author{
Won-Min Song
}