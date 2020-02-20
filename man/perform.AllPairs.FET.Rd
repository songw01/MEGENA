\name{perform.AllPairs.FET}
\alias{perform.AllPairs.FET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fisher Exact Test (FET) analysis for ALL possible pairs of two gene set lists. 
}
\description{
  FET for all pairs.
}
\usage{
  perform.AllPairs.FET(geneSets1,geneSets2,background,adjust.FET.pvalue = T)
}
\arguments{
  \item{geneSets1}{a list object containing gene signatures as character vector in each list entry.}
  \item{geneSets2}{a list object containing gene signatures as character vector in each list entry.}
  \item{background}{a character vector containing the background gene set.}
  \item{adjust.FET.pvalue}{If set TRUE, bonferroni correction is performed and output in corrected.FET.pvalue column.}
}
\details{
  Returns a data.frame entailing all comparisons. 
}
\value{
  In comparison to perform.ijPairs.FET, this function is designed to perform FET on ALL comparisons. 
}
\author{
  Won-Min Song
}
\examples{
\dontrun{
rm(list = ls())
data(Sample_Expression)
data(FrequentSets.v6.2)
ijw <- calculate.correlation(datExpr[1:100,],doPerm = 2)
el <- calculate.PFN(ijw[,1:3])
g <- graph.data.frame(el,directed = FALSE)
MEGENA.output <- do.MEGENA(g = g,remove.unsig = FALSE,doPar = FALSE,n.perm = 10)
output.summary <- MEGENA.ModuleSummary(MEGENA.output,
mod.pvalue = 0.05,hub.pvalue = 0.05,
min.size = 10,max.size = 5000,
annot.table = NULL,id.col = NULL,symbol.col = NULL,
output.sig = TRUE)

FET.results = perform.AllPairs.FET(geneSets1 = output.summary$modules,geneSets2 = msigdb.sets[[2]],background = V(g)$name)
FET.results = FET.results[order(FET.results$FET_pvalue),]
head(FET.results)

}
}