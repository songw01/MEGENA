#### module hierarchy plot functionality
globalVariables(names = c("node_size","node_value"))
plot_module_hierarchy <- function(module.table,
plot.coord = NULL,
edge.color = "grey",node.color = "black",node.label.color = "black",
label.scaleFactor = 0.5,node.scaleFactor = 0.2,arrow.size = 0.015,
data.col = NULL,low.color = "blue",mid.color = "white",high.color = "red",mid.value = 0.05)
{
	##### extract module hierarchy
	edgelist <- module.table[,c("module.parent","id")]
	module.hierarchy <- graph.data.frame(edgelist,directed = TRUE)
	#V(module.hierarchy)$label <- V(module.hierarchy)$name

	# identify roots
	k.out <- igraph::degree(module.hierarchy,mode = "out")
	k.in <- igraph::degree(module.hierarchy,mode = "in")

	root.vec <- intersect(names(k.out)[k.out > 0] ,names(k.in)[k.in == 0])

	if (length(root.vec) > 1)
	{
		final.edges <- rbind.data.frame(data.frame(module.parent = rep("root",length(root.vec)),id = root.vec),edgelist)
	}else{
		final.edges <- edgelist
		vec <- as.character(final.edges$module.parent)
		vec[vec == root.vec] <- "root"
		final.edges$module.parent <- vec;
	}

	module.hierarchy <- graph.data.frame(final.edges,directed = TRUE)

	##### get node features
	# get node layout
	if (is.null(plot.coord))
	{
		plot.coord = igraph::layout.fruchterman.reingold(module.hierarchy)
		colnames(plot.coord) = c("X1","X2");
		rownames(plot.coord) <- V(module.hierarchy)$name
	}else{
		plot.coord <- plot.coord[match(V(module.hierarchy)$name,rownames(plot.coord)),]
	}
	
	module.vsize = rep(NA,vcount(module.hierarchy));names(module.vsize) <- V(module.hierarchy)$name
	module.vsize[match(module.table$id,names(module.vsize))] <- module.table$module.size
	module.vsize[is.na(module.vsize)] <- max(module.table$module.size,na.rm = TRUE)
	module.vsize <- log10(module.vsize) * 0.5
	node.features <- data.frame(node.lab = gsub("^comp","c",V(module.hierarchy)$name),id = V(module.hierarchy)$name,node_size = module.vsize,as.data.frame(plot.coord),stringsAsFactors = FALSE)
	vert.features <- node.features

	if (!is.null(data.col))
	{
		if (is.character(data.col)) data.col <- match(data.col,colnames(module.table));
		vec <- rep(NA,nrow(node.features))
		vec[match(module.table$id,node.features$id)] <- module.table[[data.col]]
		node.features$node_value = vec;
		rm(vec)
	}
	##### get edge list
	edgelist <- get.edgelist(module.hierarchy)
	edges <- data.frame(node.features[match(edgelist[,1],node.features$id),c("X1","X2")],node.features[match(edgelist[,2],node.features$id),c("X1","X2")]);
	rownames(edges) <- NULL
	colnames(edges) <-  c("X1","Y1","X2","Y2")
	edges$midX  <- (edges$X1 + edges$X2) / 2
	edges$midY  <- (edges$Y1 + edges$Y2) / 2

	if (is.null(data.col))
	{
		hierarchy.obj <- ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),data=edges,size = 0.5,
		colour=edge.color,alpha = 0.5,arrow = arrow(length = unit(arrow.size, "npc"))) + 
		geom_point(data=node.features,mapping = aes(x = X1,y = X2,size = node_size * node.scaleFactor),colour = node.color) + 
		geom_text_repel(data = vert.features,aes(x = X1,y = X2,label = node.lab,size = node_size * label.scaleFactor),colour = node.label.color) +
		theme_bw() + 
		theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
		#panel.background = element_rect(fill = "white", colour = NA),
		axis.line=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(),
		panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		legend.position = "bottom",legend.direction = "horizontal")
	}else{
		hierarchy.obj <- ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),data=edges,size = 0.5,
		colour=edge.color,alpha = 0.5,arrow = arrow(length = unit(arrow.size, "npc"))) + 
		geom_point(data=node.features,mapping = aes(x = X1,y = X2,size = node_size * node.scaleFactor,colour = node_value)) + 
		scale_colour_gradient2(low = low.color,high = high.color,mid = mid.color,midpoint = mid.value) + 
		geom_text_repel(data = vert.features,aes(x = X1,y = X2,label = node.lab,size = node_size * label.scaleFactor),colour = node.label.color) +
		theme_bw() + 
		theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
		axis.line=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(),
		panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		legend.position = "bottom",legend.direction = "horizontal")
	}

	output <- list(hierarchy.obj = hierarchy.obj,node.data = node.features,edge.data = edges)
	return(output)
}
