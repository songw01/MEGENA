\name{module_convert_to_table}
\alias{module_convert_to_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
conversion of module list object to a data.frame table format
}
\description{
Summarizes module hub/hierarchy/membership into a data.frame table format. 
}
\usage{
module_convert_to_table(MEGENA.output,mod.pval = 0.05,
hub.pval = 0.05,min.size = 10,max.size)
}
\arguments{
  \item{MEGENA.output}{A list object. The output from "do.MEGENA()".}
  \item{mod.pval}{module compactness significance p-value, to identify modules with significant compactness.}
  \item{hub.pval}{node degree significance p-value to identify nodes with significantly high degree.}
  \item{min.size}{minimum module size allowed to finalize in the summary output.}
  \item{max.size}{maximum module size allowed to finalize in the summary output.}
  
}
\details{
the resulting data.frame contains the following essential columns: id, module.parent and module. 
If the co-expression network bears significant hubs, it will additionally have node.degree (connectivity), node.strength (sum of edge weights) and is.hub column to supplement hub information. 
}
\value{
A data.frame with the columns:
 \item{id}{gene name}
 \item{module.parent}{parent module id}
 \item{module}{module name.}
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
module.df = module_convert_to_table(MEGENA.output,mod.pval = 0.05,
hub.pval = 0.05,min.size = 10,max.size)
head(module.df)
}
}