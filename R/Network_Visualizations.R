get_module_hierarchy_graph = function(htbl,max.depth = 5,anchor.mid = NULL,
                                      h.scale = NULL,is.circular = TRUE,add_names = TRUE,
                                      layout = "dendrogram")
{
  #### This function returns hierarchy in dendrogram format
  g = graph.data.frame(htbl,directed = TRUE)
  
  # subset for M7
  if (!is.null(anchor.mid))
  {
    v = names(ego(graph = g,order = max.depth,nodes = anchor.mid,mode = "out")[[1]])
    sg = induced.subgraph(graph = g,vids = v)
  }else{
    sg = g
  }
  
  kout = igraph::degree(sg,mode = "in")
  
  if (sum(kout == 0,na.rm = TRUE) > 1) stop("Hierarchy has more than 1 root.")
  
  
  if (layout == "dendrogram")
  {
    pobj = ggraph(graph = sg, layout = layout,circular = is.circular) 
    if (!is.null(h.scale)) pobj$data$y = pobj$data$y * h.scale
    pobj = pobj + geom_edge_elbow()
  }
  
  if (layout == "treemap")
  {
    pobj = ggraph(graph = sg, layout = layout) 
    if (!is.null(h.scale)) pobj$data$y = pobj$data$y * h.scale
    pobj = pobj + geom_edge_link(arrow = arrow(length = unit(2, 'mm')),alpha = 0.5)
  }
  
  if (add_names) pobj = pobj + geom_node_label(aes(label = name),size = 6)
  
  output = list(pobj = pobj,graph.obj = sg,anchor.mid = anchor.mid)
  
  return(output)
}