\name{plot_module_hierarchy}
\alias{plot_module_hierarchy}
\title{
Plot module hierarchy
}
\description{
 visualized module hierarchical structure.
}
\usage{
plot_module_hierarchy(module.table,plot.coord = NULL,
edge.color = "grey",node.color = "black",node.label.color = "black",
label.scaleFactor = 0.5,node.scaleFactor = 0.2,arrow.size = 0.015,
data.col = NULL,low.color = "blue",mid.color = "white",
high.color = "red",mid.value = 0.05)
}
\arguments{
  \item{module.table}{output from MEGENA.ModuleSummary. Specifically $module.table component of the output.}
  \item{plot.coord}{Two column coordinate matrix. rownames must be labelled according to module.table$id.}
  \item{edge.color}{Edge color to be shown.}
  \item{node.color}{If data.col = NULL, node.color is used to color nodes in figure.}
  \item{node.label.color}{Node label color.}
  \item{label.scaleFactor}{scale number to adjust node label sizes.}
  \item{node.scaleFactor}{scale number to adjust node sizes.}
  \item{arrow.size}{scale number to arrow size.}
  \item{data.col}{A character to specify data vector to color nodes in module.table.}
  \item{low.color}{If data.col != NULL, color to be used in lower value spectrum.}
  \item{mid.color}{If data.col != NULL, color to be used in middle value spectrum.}
  \item{high.color}{If data.col != NULL, color to be used in high value spectrum.}
  \item{mid.value}{If data.col != NULL, value to define middle value spectrum.}
}
\details{
Module hierarchy plotting functionality using ggplot2. 
}
\value{
A list containing output$hierarchy.obj = ggplot2 object, output$node.data = node attributes, output$edge.data = edge attributes.
}
\author{
Won-Min Song
}
\examples{
\dontrun{
	rm(list = ls())
	data(Sample_Expression)
	ijw <- calculate.correlation(datExpr,doPerm = 2)
	el <- calculate.PFN(ijw[,1:3])
	g <- graph.data.frame(el,directed = FALSE)
	MEGENA.output <- do.MEGENA(g = g,remove.unsig = FALSE,doPar = FALSE,n.perm = 10)
	output.summary <- MEGENA.ModuleSummary(MEGENA.output,
	mod.pvalue = 0.05,hub.pvalue = 0.05,
	min.size = 10,max.size = 5000,
	annot.table = NULL,id.col = NULL,symbol.col = NULL,
	output.sig = TRUE)

	module.table = output.summary$module.table
	colnames(module.table)[1] <- "id"
	output.obj <- plot_module_hierarchy(module.table = module.table,
	label.scaleFactor = 0.15,arrow.size = 0.005,node.label.color = "blue")

	print(output.obj[[1]])
}
}