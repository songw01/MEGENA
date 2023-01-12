\name{perform.ijPairs.FET}
\alias{perform.ijPairs.FET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fisher Exact Test (FET) analysis for user specific pairs of gene sets
}
\description{
  FET with given pairs
}
\usage{
  perform.ijPairs.FET(geneSets1,geneSets2,ij,background,adjust.FET.pvalue = T)
}
\arguments{
  \item{geneSets1}{a list object containing gene signatures as character vector in each list entry.}
  \item{geneSets2}{a list object containing gene signatures as character vector in each list entry.}
  \item{ij}{two column integer matrix. Each row contains a pair of indices to compare between geneSet1 and geneSet2. For instance, ij[1,] = c(1,2), then FET is performed between geneSets1[[1]] and geneSets2[[2]].}
  \item{background}{a character vector containing the background gene set.}
  \item{adjust.FET.pvalue}{If set TRUE, bonferroni correction is performed and output in corrected.FET.pvalue column.}
}
\details{
  Returns a data.frame entailing all comparisons. 
}
\value{
  In comparison to perform.AllPairs.FET, this function is designed to perform FET on a specified set of comparisons to save time. 
}
\author{
  Won-Min Song
}