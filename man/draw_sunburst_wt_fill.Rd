\name{draw_sunburst_wt_fill}
\alias{draw_sunburst_wt_fill}
\title{
  Draw sunburst plot showing MEGENA module hierarchy.
}
\description{
  Sunburst plot and colored heatmaps 
}
\usage{
  draw_sunburst_wt_fill(module.df,
                        parent.col = "module.parent",id.col = "id",
                        min.angle = 5,
                        feat.col,
                        fill.type = "continuous",log.transform = TRUE,
                        fill.scale = NULL,
						border.col = "black",
						border.width = 0.25,
                        theme.adjust = NULL
                        )
}
\arguments{
  \item{module.df}{A data.frame table summarizing module information. Must contain module parent and child relation for hierarchy visualization.}
  \item{parent.col}{Character object, name for the parent module column in module.df.}
  \item{id.col}{Character object for the module id column in module df.}
  \item{min.angle}{Minimum angle that rectangles in the sunburst are labeled with respective module id.}
  \item{feat.col}{Chracter object, for the feature column in module.df to color the heatmaps.}
  \item{fill.type}{continuous/discrete, is the variable numeric (continuous) or factor (discrete)?}
  \item{log.transform}{TRUE/FALSE. do log10 transform for p-values?}
  \item{fill.scale}{A ggplot object to specify heatmap coloring scheme. Permissible functions are: scale_fill_gradient,scale_fill_gradient2,scale_fill_gradientn,scale_fill_manual.}
  \item{border.col}{sunburst border color. Default is black.}
  \item{border.width}{sunburst border width. Default is 0.25.}
  \item{theme.adjust}{A ggplot object to specify theme for plotting.}
}
\details{
  makes use of ggraph scheme to manipulate and draw sunburst plot in ggplot2 framework. fill.scale and theme.adjust provide flexibility to designate heatmap coloring schemes and figure aesthetics. 
}
\value{
  ggplot object for the figure
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
    
    # no coloring 
    sbobj = draw_sunburst_wt_fill(module.df = output.summary$module.table,
    feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
    sbobj
    
    # get some coloring (with log transform option)
    mdf= output.summary$module.table
    mdf$heat.pvalue = runif(nrow(mdf),0,0.1)
    
    sbobj = draw_sunburst_wt_fill(module.df = mdf,feat.col = "heat.pvalue",log.transform = TRUE,
    fill.type = "continuous",
    fill.scale = scale_fill_gradient2(low = "white",mid = "white",high = "red",
    midpoint = -log10(0.05),na.value = "white"), 
    id.col = "module.id",parent.col = "module.parent")
    sbobj
    
    # get discrete coloring done
    mdf$category = factor(sample(x = c("A","B"),size = nrow(mdf),replace = TRUE))
    sbobj = draw_sunburst_wt_fill(module.df = mdf,feat.col = "category",
    fill.type = "discrete",
    fill.scale = scale_fill_manual(values = c("A" = "red","B" = "blue")), 
    id.col = "module.id",parent.col = "module.parent")
    sbobj
  }
}