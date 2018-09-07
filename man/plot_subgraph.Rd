\name{plot_subgraph}
\alias{plot_subgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
subnetwork plotting functionality. 
}
\description{
A modification of plot_module() function for more general subnetwork plotting purpose.
}
\usage{
plot_subgraph(module,hub = NULL,PFN,node.default.color = "black",
gene.set = NULL,color.code = "grey",show.legend = TRUE,
label.hubs.only = TRUE,hubLabel.col = "red",hubLabel.sizeProp = 0.5,show.topn.hubs = 10,
node.sizeProp = 13,label.sizeProp = 13,label.scaleFactor = 10,layout = "kamada.kawai")
}
\arguments{
  \item{module}{A character vector containing gene names to be subsetted.}
  \item{hub}{If provided, genes in hub will be highlighted as triangles in resulting figure.}
  \item{PFN}{igraph object retaining PFN topology.}
  \item{node.default.color}{Default node colors for those that do not intersect with signatures in gene.set.}
  \item{gene.set}{A list object containing signatures for customized coloring of nodes in resulting network plot.}
  \item{color.code}{A character vector with matched length to "gene.set", to specify colors for each signature.}
  \item{show.legend}{TRUE/FALSE for showing node legend on the bottom of the figure.}
  \item{label.hubs.only}{TRUE/FALSE to show labels for significant hub genes only, or all genes. Defauly is TRUE.}
  \item{hubLabel.col}{Label color for hubs. Default is "red"}
  \item{hubLabel.sizeProp}{A multiplicative factor to adjust hub label sizes with respect to node size values. Default is 0.5}
  \item{show.topn.hubs}{Maximal number of hubs to label on module subnetwork. Default is 10.}
  \item{node.sizeProp}{A multiplicative factor to adjust node sizes with respect to 90th percentile degree node size. Default is 13}
  \item{label.sizeProp}{A multiplicative factor to adjust node label sizes with respect to 90th percentile degree node size. Default is 13}
  \item{label.scaleFactor}{Overall scale factor to control the final size of node labels appearing in figure. Default is 10.}
  \item{layout}{Network layout algorithm to apply. Options are: "kamada.kawai", "fruchterman.reingold".}
  
}
\details{
Subnetwork plot functionality with application of "ggrepel" package for node labeling. The most effective way to control overall node label size is through label.scaleFactor. 
}
\value{
A list object holding ggplot object and node annotation table.
 }
\author{
Won-Min Song
}
\examples{
\dontrun{
rm(list = ls())
library(MEGENA)

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

pnet.obj <- plot_subgraph(module = output.summary$modules[[1]],
hub = c("CD3E","CD2"),PFN = g,node.default.color = "black",
gene.set = NULL,color.code = c("grey"),show.legend = TRUE,
label.hubs.only = TRUE,hubLabel.col = "red",hubLabel.sizeProp = 0.5,
show.topn.hubs = 10,node.sizeProp = 13,label.sizeProp = 13,
label.scaleFactor = 10,layout = "kamada.kawai")

# the plot	
pnet.obj[[1]]

# the annotation
pnet.obj[[2]]
}
}