
extract.module_with_hub <- function(output.summary)
{
	# extract hub-module matched results
	modules <- hubs <- vector("list",nrow(output.summary$module.table))
	names(modules) <- names(hubs) <- as.character(output.summary$module.table[[1]])
	for (i in 1:nrow(output.summary$module.table))
	{
		module.id <- as.character(output.summary$module.table[[1]][i]);
		modules[[i]] <- output.summary$modules[[which(names(output.summary$modules) == module.id)]];
		hubs[[i]] <- gsub("\\((.*)","",do.call('c',strsplit(as.character(output.summary$module.table$module.hub[i]),",")))
	}

	list(modules = modules,hubs = hubs)
}
globalVariables(names = c("X1","X2","Y1","Y2","id","n.cores","node.lab","node.stat"))
## gene set based module plotting
plot_module <- function(output.summary,PFN,subset.module = NULL,col.names,
gene.set = NULL,color.code = "logFC",show.legend = TRUE,
label.hubs.only = TRUE,hubLabel.col = "red",hubLabel.sizeProp = 0.5,show.topn.hubs = 10,
node.sizeProp = 13,label.sizeProp = 13,label.scaleFactor = 10,label.alpha = 0.5,
layout = "kamada.kawai",output.plot = TRUE,out.dir = "modulePlot")
{
	# get matched lists of modules and hubs
	matched.obj <- extract.module_with_hub(output.summary)
		
	# subset modules
	if (is.null(subset.module)) subset.module <- names(matched.obj[[1]])
	
	pnet <- vector("list",length(subset.module));names(pnet) <- subset.module
	
	if (output.plot) dir.create(out.dir);
	
	for (i in 1:length(subset.module))
	{
		cat(paste("Processing module:",subset.module[i],"\n",sep = ""))
		module <- matched.obj$modules[[match(subset.module[i],names(matched.obj$modules))]]
		hub <- matched.obj$hubs[[match(subset.module[i],names(matched.obj$hubs))]]
		modlen <- length(module)
		cat(paste(" - # of genes:",length(module),"\n"))
		cat(paste(" - # of hubs:",length(hub),"\n"))
		
		if (!is.null(show.topn.hubs)) hub <- hub[1:min(c(length(hub),show.topn.hubs))]
		
		### iterate for each i
		sub.PFN <- induced.subgraph(graph = PFN,vids = module)
		edgelist <- get.edgelist(sub.PFN)

		## create basic node features to show children split and key drivers
		# a. get children modules
		child.module <- as.character(output.summary$module.table[[1]][which(output.summary$module.table$module.parent == subset.module[i])])

		# b. assign node features
		node.col <- rep(paste("p:",subset.module[i],sep = ""),vcount(sub.PFN))
		node.shape <- rep("gene",vcount(sub.PFN));
		node.size <- log(igraph::degree(sub.PFN));names(node.size) <- V(sub.PFN)$name
		denom <- quantile(node.size,probs = 0.6);
		#if (denom == 0) denom = 1;
		node.size <- node.size/denom * node.sizeProp;node.size <- node.size[V(sub.PFN)$name]
				
		if (length(hub) > 0) node.size[which(V(sub.PFN)$name %in% hub)] <- node.size[which(V(sub.PFN)$name %in% hub)] * 1.5

		color.map <- c("grey")
		names(color.map)[1] <- paste("p:",subset.module[i],sep = "")
		if (length(child.module) > 0 & is.null(gene.set))
		{
		 child.mod <- matched.obj$modules[child.module]
		 child.hub <- matched.obj$hubs[child.module]
		 for (j in 1:length(child.mod))
		 {
		  node.col[match(child.mod[[j]],V(sub.PFN)$name)]  <- paste("c:",names(child.mod)[j],sep = "")
		  node.shape[match(child.hub[[j]],V(sub.PFN)$name)]  <- "child.hub"
		 }
		 
		 child.map <- col.names[1:length(child.mod)];names(child.map) <- paste("c:",names(child.mod),sep = "")
		 color.map <- c(color.map,child.map)
		}
		node.shape[which(V(sub.PFN)$name %in% hub)] <- "module.hub"

		## create co-ordinates 
		#if (is.null(plot.coord))
		#{
		if (layout == "kamada.kawai")
		{
		  plot.coord = igraph::layout_with_kk(graph = sub.PFN,kkconst = vcount(sub.PFN),maxiter = 50 * vcount(sub.PFN))
		}
		if (layout == "fruchterman.reingold")
		{
		  plot.coord = igraph::layout.fruchterman.reingold(sub.PFN)
		}
		
		colnames(plot.coord) = c("X1","X2");
		
		node.features <- data.frame(id = V(sub.PFN)$name,node.size = node.size,node.shape = factor(node.shape),child.module = factor(node.col),
		as.data.frame(plot.coord),stringsAsFactors = FALSE)
		rm(node.size,node.shape,node.col,plot.coord)
		
		if (length(grep("\\|",V(sub.PFN)$name)) > 0)
		{
		 nid <- gsub("\\|(.*)","",V(sub.PFN)$name)
		 node.features <- data.frame(node.lab = nid,id = V(sub.PFN)$name,node.features[,-1])
		 rm(nid)
		}else{
		 node.features <- data.frame(node.lab = V(sub.PFN)$name,id = V(sub.PFN)$name,node.features[,-1])
		}
		
		## combine with node attributes
		# create edges data.frame object for segment plotting
		edges <- data.frame(node.features[match(edgelist[,1],node.features$id),c("X1","X2")],node.features[match(edgelist[,2],node.features$id),c("X1","X2")]);
		rownames(edges) <- NULL
		colnames(edges) <-  c("X1","Y1","X2","Y2")
		edges$midX  <- (edges$X1 + edges$X2) / 2
		edges$midY  <- (edges$Y1 + edges$Y2) / 2

		# create ggplot object
		cat("- generating module subnetwork figure...\n")
		title.size <- 20;
		if (length(module) < 50) title.size = 14
		label.size.factor <- 0.75;
		if (length(module) < 50) label.size.factor <- 0.4
		node.frame.shape <- rep(1,nrow(node.features));node.frame.shape[which(node.features$node.shape == "child.hub")] <- 5;
		node.frame.shape[which(node.features$node.shape == "module.hub")] <- 2;
		
		do.set.color <- FALSE
		if (!is.null(gene.set)) 
		{
			intsct <- lapply(gene.set,function(x,y) intersect(x,y),y = as.character(node.features$id))
			if (any(sapply(intsct,length) > 0)) do.set.color = TRUE
		}
		if (!is.null(gene.set) & do.set.color)
		{
			stat.vec <- rep("NA",nrow(node.features))
			for (cc in 1:length(gene.set))
			{
				stat.vec[which(node.features$id %in% gene.set[[cc]])] <- names(gene.set)[cc]
			}
			node.features$node.stat = stat.vec;
			node.stat.map <- color.code;names(node.stat.map) <- names(gene.set)
			pnet[[i]] <- ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),data=edges, size = 0.5, colour="grey",
			alpha = 0.5) +
			geom_point(data=node.features,mapping = aes(x = X1,y = X2,colour = node.stat,size = node.size * 40,shape = node.shape)) +
			scale_shape_manual(values = c("gene" = 16,"child.hub" = 18,"module.hub" = 17)) + 
			scale_colour_manual(values = c("NA" = "grey",node.stat.map)) +
			labs(title = subset.module[i]) + 
			theme_bw() + 
			theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
						axis.line=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(),
						panel.grid.minor = element_blank(), panel.grid.major = element_blank(),plot.title = element_text(size = title.size),
						legend.position = "bottom",legend.direction = "horizontal") +
			guides(size = FALSE)
		}else{
		
			pnet[[i]] <- ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),data=edges, size = 0.5, colour="grey",alpha = 0.5) +
			geom_point(data=node.features,mapping = aes(x = X1,y = X2,colour = child.module,size = node.size * 40,shape = node.shape)) +
			scale_shape_manual(values = c("gene" = 16,"child.hub" = 18,"module.hub" = 17)) + 
			scale_colour_manual(values = color.map,guide = "legend") + 
			labs(title = subset.module[i]) + 
			theme_bw() + 
			theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
						axis.line=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(),
						panel.grid.minor = element_blank(), panel.grid.major = element_blank(),plot.title = element_text(size = title.size),
						legend.position = "bottom",legend.direction = "horizontal") +
			guides(size = FALSE)
		}
		
		if (!show.legend) pnet[[i]] <- pnet[[i]] + guides(colour = FALSE,shape = FALSE)
		
		if (label.hubs.only & length(hub) > 0)
		{
			pnet[[i]] <- pnet[[i]] + geom_text_repel(data = subset(node.features,id %in% hub),
			   aes(x = X1,y = X2,label = node.lab,size = node.size * hubLabel.sizeProp * label.scaleFactor),colour = hubLabel.col,alpha = label.alpha)
		}else{
			vert.features = node.features;
			vert.features$node.size = vert.features$node.size * label.sizeProp;
			pnet[[i]] <- pnet[[i]] + geom_text_repel(data = vert.features,
			   aes(x = X1,y = X2,label = node.lab,size = node.size * label.scaleFactor),colour = "black",alpha = label.alpha)
		}
		rm(edges,node.features,hub,module)
		
		# output to .tiff file
		if (output.plot)
		{
			cat(paste("- output figure to:",paste(out.dir,"/",subset.module[i],".png",sep = ""),"\n",sep = ""))
			if (modlen > 100)
			{
				if (modlen > 500)
				{
					png(filename = paste(out.dir,"/",subset.module[i],".png",sep = ""),
					width = 70 * round(sqrt(vcount(sub.PFN))) + 100,
					height = 70 * round(sqrt(vcount(sub.PFN))),res = 300)
					print(pnet[[i]]);
					dev.off()
				}else{
					png(filename = paste(out.dir,"/",subset.module[i],".png",sep = ""),
					width = 200 * round(sqrt(vcount(sub.PFN))),
					height = 200 * round(sqrt(vcount(sub.PFN))) + 150,res = 300)
					print(pnet[[i]]);
					dev.off()
				}
			}else{
				if (modlen > 10)
				{
					png(filename = paste(out.dir,"/",subset.module[i],".png",sep = ""),
					width = 700 * round(sqrt(vcount(sub.PFN))) + 600,
					height = 700 * round(sqrt(vcount(sub.PFN))),res = 300)
					print(pnet[[i]]);
					dev.off()
				}else{
					png(filename = paste(out.dir,"/",subset.module[i],".png",sep = ""),width = 1700 + 500,
					height = 1700,res = 300)
					print(pnet[[i]]);
					dev.off()

				}
			}
		}
	}	
	return(pnet)
}

plot_subgraph <- function(module,hub = NULL,PFN,node.default.color = "black",
gene.set = NULL,color.code = "grey",show.legend = TRUE,
label.hubs.only = TRUE,hubLabel.col = "red",hubLabel.sizeProp = 0.5,show.topn.hubs = 10,
node.sizeProp = 13,label.sizeProp = 13,label.scaleFactor = 10,layout = "kamada.kawai")
{
		
	if (is.null(hub)) hub <- c()
	modlen <- length(module)
	cat(paste(" - # of genes:",length(module),"\n"))
	cat(paste(" - # of hubs:",length(hub),"\n"))
	
	sub.PFN <- induced.subgraph(graph = PFN,vids = module)
	edgelist <- get.edgelist(sub.PFN)
	
	###### get node features
	# get node size according to degree
	node.size <- log(igraph::degree(sub.PFN));
	denom <- quantile(node.size,probs = 0.6);
	if (denom == 0) denom = 1;
	node.size <- node.size/denom * node.sizeProp;node.size <- node.size[V(sub.PFN)$name]
		
	if (length(hub) > 0) node.size[which(V(sub.PFN)$name %in% hub)] <- node.size[which(V(sub.PFN)$name %in% hub)] * 1.5
	
	# node shape
	node.shape <- rep("gene",vcount(sub.PFN));
	node.shape[V(sub.PFN)$name %in% hub] <- "hub"
		
	# get layout of network
	if (layout == "kamada.kawai")
	{
		plot.coord = layout.kamada.kawai(sub.PFN)
	}
	if (layout == "fruchterman.reingold")
	{
		plot.coord = layout.fruchterman.reingold(sub.PFN)
	}
	colnames(plot.coord) = c("X1","X2");
	
	# collect attributes
	node.features <- data.frame(id = V(sub.PFN)$name,node.size = node.size,node.shape = factor(node.shape),	as.data.frame(plot.coord),stringsAsFactors = FALSE)
	rm(node.size,node.shape,plot.coord)
	
	# add label. 
	if (length(grep("\\|",V(sub.PFN)$name)) > 0)
	{
		nid <- gsub("\\|(.*)","",V(sub.PFN)$name)
		node.features <- data.frame(node.lab = nid,id = V(sub.PFN)$name,node.features[,-1])
		rm(nid)
	}else{
		node.features <- data.frame(node.lab = V(sub.PFN)$name,id = V(sub.PFN)$name,node.features[,-1])
	}
	
	# create edges data.frame object for segment plotting
	edges <- data.frame(node.features[match(edgelist[,1],node.features$id),c("X1","X2")],node.features[match(edgelist[,2],node.features$id),c("X1","X2")]);
	rownames(edges) <- NULL
	colnames(edges) <-  c("X1","Y1","X2","Y2")
	edges$midX  <- (edges$X1 + edges$X2) / 2
	edges$midY  <- (edges$Y1 + edges$Y2) / 2

	# create ggplot object
	cat("- generating module subnetwork figure...\n")
	title.size <- 20;
	if (vcount(sub.PFN) < 50) title.size = 14
	label.size.factor <- 0.75;
	if (vcount(sub.PFN) < 50) label.size.factor <- 0.4
	node.frame.shape <- rep(1,nrow(node.features));
	node.frame.shape[which(node.features$node.shape == "hub")] <- 2;
		
	do.set.color <- FALSE
	if (!is.null(gene.set)) 
	{
		intsct <- lapply(gene.set,function(x,y) intersect(x,y),y = as.character(node.features$id))
		if (any(sapply(intsct,length) > 0)) do.set.color = TRUE
	}
	
	if (!is.null(gene.set) & do.set.color)
	{
		stat.vec <- rep("NA",nrow(node.features))
		for (cc in 1:length(gene.set))
		{
			stat.vec[which(node.features$id %in% gene.set[[cc]])] <- names(gene.set)[cc]
		}
		node.features$node.stat = stat.vec;
		node.stat.map <- color.code;names(node.stat.map) <- names(gene.set)
		pnet <- ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),data=edges, size = 0.5, colour="grey",
		alpha = 0.5) +
		geom_point(data=node.features,mapping = aes(x = X1,y = X2,colour = node.stat,size = node.size * 40,shape = node.shape)) +
		scale_shape_manual(values = c("gene" = 16,"hub" = 17)) + 
		scale_colour_manual(values = c("NA" = node.default.color,node.stat.map)) +
		#labs(title = subset.module[i]) + 
		theme_bw() + 
		theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
			axis.line=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(),
			panel.grid.minor = element_blank(), panel.grid.major = element_blank(),plot.title = element_text(size = title.size),
			legend.position = "bottom",legend.direction = "horizontal") +
		guides(size = FALSE)
	}else{
		pnet <- ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),data=edges, size = 0.5, colour="grey",alpha = 0.5) +
			geom_point(data=node.features,colour = node.default.color,mapping = aes(x = X1,y = X2,size = node.size * 40,shape = node.shape)) +
			scale_shape_manual(values = c("gene" = 16,"hub" = 17)) + 
			#scale_colour_manual(values = color.map,guide = "legend") + 
			#labs(title = subset.module[i]) + 
			theme_bw() + 
			theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
						axis.line=element_blank(), axis.ticks=element_blank(), axis.text = element_blank(),
						panel.grid.minor = element_blank(), panel.grid.major = element_blank(),plot.title = element_text(size = title.size),
						legend.position = "bottom",legend.direction = "horizontal") +
			guides(size = FALSE)
	}
		
	if (!show.legend) pnet <- pnet + guides(colour = FALSE,shape = FALSE)
		
	if (label.hubs.only & length(hub) > 0)
	{
		pnet <- pnet + geom_text_repel(data = subset(node.features,id %in% hub),
		   aes(x = X1,y = X2,label = node.lab,size = node.size * hubLabel.sizeProp * label.scaleFactor),colour = hubLabel.col)
	}else{
		vert.features = node.features;
		vert.features$node.size = vert.features$node.size * label.sizeProp;
		pnet <- pnet + geom_text_repel(data = vert.features,
		   aes(x = X1,y = X2,label = node.lab,size = node.size * label.scaleFactor),colour = "black")
	}
	rm(edges,hub,module)
	
	output <- list(pnet = pnet,node.features = node.features)
	return(output)
}