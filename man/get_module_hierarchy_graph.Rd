\name{get_module_hierarchy_graph}
\alias{get_module_hierarchy_graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Module Hierarchy Visualization
}
\description{
  Utilizes various layout algorithms from ggraph to visualize module hierarchy with desired root. 
}
\usage{
  get_module_hierarchy_graph(htbl,max.depth = 5,anchor.mid = NULL,h.scale = NULL,is.circular = TRUE,layout = "dendrogram")
}
\arguments{
  \item{htbl}{Two column data.frame table capturing module hierarchy. First column is the parent module, and second column is an immediate child module.}
  \item{max.depth}{Maximum number of hierarchy depth to propagate. Default value 5 is often sufficient to capture the entire MEGENA hierarchy.}
  \item{anchor.mid}{Default is NULL. If specified, it is taken as the root module, and its module hierarchy is subsetted.}
  \item{h.scale}{When layout = "dendrogram" is used, it scales the heights of the dendrogram branching points to make the plot less busy.}
  \item{is.circular}{Default is TRUE. If TRUE, circular dendrogram is drawn.}
  \item{layout}{Options are "dendrogram" and "treemap". See ?ggraph for layout option descriptions.}
  \item{add_names}{If TRUE, it adds modules names at respective branching points.}
}
\details{
  Returns a list containing the results. 
}
\value{
  A list object with the components:
    \item{pobj}{ggplot object for the hierarchy graph.}
  \item{graph.obj}{igraph object for the hierarchy graph.}
  \item{anchor.mid}{The root module id.}
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
    tbl = output.summary$module.table
    htbl = tbl[,c("module.parent","module.id")]
    hplot = get_module_hierarchy_graph(htbl,max.depth = 5,anchor.mid = NULL,h.scale = NULL,is.circular = FALSE,layout = "dendrogram",add_names = TRUE)
    print(hplot$pobj)
  }
}