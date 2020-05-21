globalVariables(names = c("fill.value","sb.angle","x","y","name","rotation"))
draw_sunburst_wt_fill <- function(module.df,# module table
                                  ### main inputs to specify hierarchy
                                  parent.col = "module.parent",# module parent column
                                  id.col = "id",# module id columnm
                                  min.angle = 5,# minimum angle of section to label for module id
                                  # fill aspects
                                  feat.col,# feature column (character) to color sunburst
                                  fill.type = "continuous",# continuous/discrete, is the variable numeric or factor? 
                                  log.transform = TRUE,# TRUE/FALSE, do log10 transform for p-values?
                                  fill.scale = NULL,
                                  # colour aspects
                                  border.col = "black", # sunburst border color
                                  border.width = 0.25, # sunburst border line width
                                  theme.adjust = NULL)

{
  #require(ggraph)
  #############
  # make ggraph object
  hier.df = module.df[,c(parent.col,id.col)]
  htree = graph_from_data_frame(hier.df,vertices = union(hier.df[[1]],hier.df[[2]]))
  
  ## create initial sunburst
  sunb = ggraph(htree, 'partition', circular = TRUE)
  
  # filling aspect
  if (!is.null(feat.col))
  {
    if (fill.type == "continuous")
    {
      # process 
      vec = module.df[[which(colnames(module.df) == feat.col)]]
      fill.lab = feat.col
      if (!is.numeric(vec)) stop(paste("variable:",feat.col," is not numeric.",sep = "")) 
      if (log.transform) 
      {
        vec = -log10(vec); 
        fill.lab = paste("-log10(",feat.col,")",sep = "")
      }
      sunb$data$fill.value = vec[match(sunb$data$name,module.df[[which(colnames(module.df) == id.col)]])]
      sunb = sunb + geom_node_arc_bar(aes(fill = fill.value),size = border.width,colour = border.col) + 
        guides(fill = guide_colorbar(title = fill.lab)) + fill.scale
    }
    if (fill.type == "discrete")
    {
      vec = module.df[[which(colnames(module.df) == feat.col)]]
      fill.lab = feat.col
      if (!is.factor(vec)) stop(paste("variable:",feat.col," is not factor.",sep = "")) 
      sunb$data$fill.value = vec[match(sunb$data$name,module.df[[which(colnames(module.df) == id.col)]])]
      sunb = sunb + geom_node_arc_bar(aes(fill = fill.value),size = border.width,colour = border.col) + 
        guides(fill = guide_legend(title = fill.lab)) + fill.scale
    }
    
  }else{
    sunb = sunb + geom_node_arc_bar()
  }
  
  # add module id
  # add label rotation feature
  sunb$data$rotation = 180 - (90 + (sunb$data$start + sunb$data$end )/2 * (180/pi) ) - 90
  sunb$data$rotation[sunb$data$depth == 0] = 0
  sunb$data$sb.angle = abs(sunb$data$end - sunb$data$start) * (180/pi) 
  sunb$data$rotation[which(sunb$data$sb.angle < 10)] = sunb$data$rotation[which(sunb$data$sb.angle < 10)] + 90 
  sunb = sunb + geom_text(data = subset(sunb$data,sb.angle >= min.angle),mapping = aes(x = x,y = y,label = name,angle = rotation))
  
  # add theme
  if (!is.null(theme.adjust)) sunb = sunb + theme.adjust
  
  return(sunb)
}


