\name{plot_module}
\alias{plot_module}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Module plotting function.
}
\description{
Extract subnetworks for modules and plot. 
}
\usage{
plot_module(output.summary,PFN,subset.module = NULL,col.names,
gene.set = NULL,color.code = "logFC",show.legend = TRUE,
label.hubs.only = TRUE,hubLabel.col = "red",hubLabel.sizeProp = 0.5,show.topn.hubs = 10,
node.sizeProp = 13,label.sizeProp = 13,label.scaleFactor = 10,label.alpha = 0.5,
layout = "kamada.kawai",output.plot = TRUE,out.dir = "modulePlot")
}
\arguments{
  \item{output.summary}{output from summary function, "MEGENA.ModuleSummary".}
  \item{PFN}{igraph object retaining PFN topology.}
  \item{subset.module}{A character vector for list of module names to plot. Default = NULL plots all modules in output.summary.}
  \item{col.names}{a character vector for list of colors to be used for coloring children modules.}
  \item{gene.set}{A list object containing signatures for customized coloring of nodes in resulting network plot.}
  \item{color.code}{A character vector with matched length to "gene.set", to specify colors for each signature.}
  \item{label.hubs.only}{TRUE/FALSE to show labels for significant hub genes only, or all genes. Defauly is TRUE.}
  \item{hubLabel.col}{Label color for hubs. Default is "red"}
  \item{show.legend}{TRUE/FALSE for showing node legend on the bottom of the figure.}
  \item{hubLabel.sizeProp}{A multiplicative factor to adjust hub label sizes with respect to node size values. Default is 0.5}
  \item{show.topn.hubs}{Maximal number of hubs to label on module subnetwork. Default is 10.}
  \item{node.sizeProp}{A multiplicative factor to adjust node sizes with respect to 90th percentile degree node size. Default is 13}
  \item{label.sizeProp}{A multiplicative factor to adjust node label sizes with respect to 90th percentile degree node size. Default is 13}
  \item{label.scaleFactor}{Overall scale factor to control the final size of node labels appearing in figure. Default is 10.}
  \item{label.alpha}{Transparency value ranging from 0 (transparent) to 1 (solid). Default is 0.5.}
  \item{layout}{Network layout algorithm to apply. Options are: "kamada.kawai", "fruchterman.reingold".}
  \item{output.plot}{logical value. output.plot = TRUE generates figure files under folder, "modulePlot".}
  \item{out.dir}{if output.plot = TRUE, then out.dir is created and resulting figures are exported to .png files to the folder.}
}
\details{
Subnetwork plot functionality with application of "ggrepel" package for node labeling. The most effective way to control overall node label size is through label.scaleFactor. 
}
\value{
A list object holding ggplot objects for plotted modules. 
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

pnet.obj <- plot_module(output = output.summary,PFN = g,subset.module = "comp1_2",
	layout = "kamada.kawai",label.hubs.only = FALSE,
	gene.set = list("hub.set" = c("CD3E","CD2")),color.code =  c("red"),
	output.plot = FALSE,out.dir = "modulePlot",col.names = c("grey","grey","grey"),
	hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)
	
pnet.obj
}
}