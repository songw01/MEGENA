\name{MEGENA.ModuleSummary}
\alias{MEGENA.ModuleSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MEGENA module summary
}
\description{
Summarizes modules into a table. 
}
\usage{
MEGENA.ModuleSummary(MEGENA.output,
mod.pvalue = 0.05,hub.pvalue = 0.05,
min.size = 10,max.size = 2500,
annot.table = NULL,symbol.col = NULL,id.col = NULL, 
output.sig = TRUE)
}
\arguments{
  \item{MEGENA.output}{A list object. The output from "do.MEGENA()".}
  \item{mod.pvalue}{module compactness significance p-value, to identify modules with significant compactness.}
  \item{hub.pvalue}{node degree significance p-value to identify nodes with significantly high degree.}
  \item{min.size}{minimum module size allowed to finalize in the summary output.}
  \item{max.size}{maximum module size allowed to finalize in the summary output.}
  \item{annot.table}{Default value is NULL, indicating no mapping is provided between node names to gene symbols. 
  If provided, the mapping between node names (id.col) and gene symbol (symbol.col) are used.}
  \item{id.col}{column index of annot.table for node names.}
  \item{symbol.col}{column index of annot.table for gene symbols.}
  \item{output.sig}{Default value is TRUE, indicating significant modules are outputted.}
}
\details{
output$module.table contains many important information including module hierarchy, as indicated by 
}
\value{
A list object with the components:
 \item{modules}{Final set of modules obtained upon apply mod.pvalue for significance, min.size and max.size for module size thresholding.}
 \item{mapped.modules}{gene symbol mapped modules when "annot.table" is provided.}
 \item{module.table}{data.frame object for module summary table. Columns include: id, module.size, module.parent, module.hub, module.scale and module.pvalue.}
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
output.summary <- MEGENA.ModuleSummary(MEGENA.output,
mod.pvalue = 0.05,hub.pvalue = 0.05,
min.size = 10,max.size = 5000,
annot.table = NULL,id.col = NULL,symbol.col = NULL,
output.sig = TRUE)
}
}