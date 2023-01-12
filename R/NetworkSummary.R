get.hub.summary <- function(MEGENA.output)
{
	all.genes <- Reduce("union",MEGENA.output$module.output$modules)
	all.genes <- Reduce("union",MEGENA.output$hub.output$hub.list)
	hub.member <- do.call(cbind,lapply(MEGENA.output$hub.output$hub.list,function(x,y) {out <- rep(0,length(y));out[which(y %in% x)] <- 1;return(out)},y = all.genes))
	hub.freq <- apply(hub.member,1,function(x) sum(x,na.rm = T))
	hub.summary <- apply(hub.member,1,function(x,y) paste(y[which(x == 1)],collapse = ","),y = colnames(hub.member))

	hub.Df <- data.frame(node = all.genes,as.data.frame(hub.member),frequency = hub.freq,scale.summary = hub.summary)
	return(hub.Df)
}

################################### module summary function
MEGENA.ModuleSummary <- function(MEGENA.output,
mod.pvalue = 0.05,hub.pvalue = 0.05,
min.size = 10,max.size = 2500,
annot.table = NULL,symbol.col = NULL,id.col = NULL, # attributes to use when annotation from id -> gene symbol is provided
output.sig = TRUE)
{
	module.output = MEGENA.output$module.output
	hub.output = MEGENA.output$hub.output

	# extract hubs
	run.hub = all(!is.null(MEGENA.output$hub.output)) & all(!is.na(MEGENA.output$hub.output))
	if (run.hub)
	{
		if (is.null(annot.table))
		{
			hub.summary <- sapply(hub.output$module.degreeStat,function(x,pval) {
														 tbl <- x[which(x$pvalue < hub.pvalue),];
														 tbl <- tbl[order(tbl$degree,decreasing = T),];
														 paste(paste(as.character(tbl$gene),"(",as.character(tbl$degree),")",sep = ""),collapse = ",")
														},pval = hub.pvalue)
		}else{
		 gene.vec <- as.character(annot.table[[symbol.col]]);
		 names(gene.vec) <- as.character(annot.table[[id.col]])
		 hub.summary <- sapply(hub.output$module.degreeStat,function(x,pval,gene.vec) {
																	 tbl <- x[which(x$pvalue < hub.pvalue),];
																	 tbl <- tbl[order(tbl$degree,decreasing = T),];
																	 matched.symbol <- gene.vec[match(as.character(tbl$gene),names(gene.vec))]
																	 tbl$gene <- paste(matched.symbol,as.character(tbl$gene),sep = "|")
																	 paste(paste(as.character(tbl$gene),"(",as.character(tbl$degree),")",sep = ""),collapse = ",")
																	},pval = hub.pvalue,gene.vec = gene.vec)
		}
	}else{
		hub.summary <- rep(NA,length(module.output$modules));names(hub.summary) <- names(module.output$modules)
	}
	
	# module pvalue
	module.pvalue <- module.output$module.pvalue;names(module.pvalue) <- names(module.output$modules)

	# module size 
	module.size <- sapply(module.output$modules,length)

	# module relation
	module.parent <- rep(NA,length(module.output$modules))
	names(module.parent) <- names(module.output$modules)
	module.parent[names(module.output$modules)[module.output$module.relation[,2]]] <- names(module.output$modules)[module.output$module.relation[,1]]

	if (all(!is.null(hub.output)) & all(!is.na(hub.output)))
	{
		# module by scale-groups
		module.scaleGroup <- lapply(hub.output$scale.summary.clusters,names)
		scale.group <- rep(NA,length(module.output$modules))
		names(scale.group) <- names(module.output$modules)
		for (i in 1:length(module.scaleGroup))
		{
		 scale.group[module.scaleGroup[[i]]] <- names(module.scaleGroup)[i]
		}

		module.id <- names(module.output$modules)
		module.table <- data.frame(module.id = module.id,module.size = module.size[module.id],module.parent = module.parent[module.id],
		module.hub = hub.summary[module.id],module.scale = scale.group[module.id],module.pvalue = module.pvalue[module.id])
	}else{
		module.id <- names(module.output$modules)
		module.table <- data.frame(module.id = module.id,module.size = module.size[module.id],module.parent = module.parent[module.id],
		module.hub = hub.summary[module.id],module.scale = rep(NA,length(module.id)),module.pvalue = module.pvalue[module.id])
	}
	

	if (output.sig)
	{
	 
	 sig.modules <- as.character(module.table$module.id[which(module.table$module.pvalue < mod.pvalue & module.table$module.size >= min.size &
	 module.table$module.size <= max.size)])
	 
	 modules <- module.output$modules[sig.modules]
	 
	 if (!is.null(annot.table))
	 {
	  gene.vec <- as.character(annot.table[[symbol.col]]);
	  names(gene.vec) <- as.character(annot.table[[id.col]])
	  n.modules <- lapply(modules,function(x,gene.vec) paste(gene.vec[match(x,names(gene.vec))],x,sep = "|"),gene.vec = gene.vec)
	  output <- list(modules = modules,mapped.modules = n.modules,module.table = module.table[which(module.table$module.id %in% sig.modules),])
	 }else{
	  output <- list(modules = modules,module.table = module.table[which(module.table$module.id %in% sig.modules),])
	 }
	 	 
	}else{
	 
	 modules <- module.output$modules
	 
	 if (!is.null(annot.table))
	 {
	  gene.vec <- as.character(annot.table[[symbol.col]]);
	  names(gene.vec) <- as.character(annot.table[[id.col]])
	  n.modules <- lapply(modules,function(x,gene.vec) paste(gene.vec[match(x,names(gene.vec))],x,sep = "|"),gene.vec = gene.vec)
	  output <- list(modules = modules,mapped.modules = n.modules,module.table = module.table[which(module.table$module.id %in% names(modules)),])
	 }else{
	  output <- list(modules = modules,module.table = module.table[which(module.table$module.id %in% names(modules)),])
	 }
	}

	return(output)
}


###### function to summarize MEGENA clusters and PFN informations.

node.summary <- function (PFN, module.output, hub.output, module.pval = 0.05,hub.pval = 0.05)
{
    group.summary <- hub.output$scale.cluster.info
    module.table <- data.frame(probe = V(PFN)$name)
    for (i in 1:nrow(group.summary)) {
        min.id <- paste(as.character(group.summary[[1]])[i],
            "(min.alpha=", group.summary[[2]][i], ")", sep = "")
        max.id <- paste(as.character(group.summary[[1]])[i],
            "(max.alpha=", group.summary[[3]][i], ")", sep = "")
		if (!is.infinite(group.summary[[2]][i]))
		{
			min.module <- get.union.cut(module.output = module.output,
				alpha.cut = group.summary[[2]][i], module.pval = module.pval,
				output.plot = F)
		}else{
			min.module <- module.output$modules[which(is.infinite(module.output$module.alpha))]
		}
        min.part <- rep(NA, vcount(PFN))
        for (j in 1:length(min.module)) min.part[which(V(PFN)$name %in%
            min.module[[j]])] <- names(min.module)[j]
        
		if (!is.infinite(group.summary[[3]][i]))
		{
			max.module <- get.union.cut(module.output = module.output,
				alpha.cut = group.summary[[3]][i], module.pval = module.pval,
				output.plot = F)
		}else{
			max.module <- module.output$modules[which(is.infinite(module.output$module.alpha))]
		}
		
        max.part <- rep(NA, vcount(PFN))
        for (j in 1:length(max.module)) max.part[which(V(PFN)$name %in%
            max.module[[j]])] <- names(max.module)[j]
        mat <- cbind(min.part, max.part)
        colnames(mat) <- c(min.id, max.id)
        module.table <- cbind.data.frame(module.table, as.data.frame(mat))
    }
    hubs <- lapply(hub.output$module.degreeStat, function(x,
        y) as.character(x[[1]][which(x[[3]] < y)]), y = hub.pval)
    sig.modules <- module.output$modules[which(module.output$module.pvalue <=
        module.pval)]
    node.topo <- get.nodeTopoInfo(PFN = PFN, modules = sig.modules,
        KD = hubs, label.pttrn = NULL, outputfname = NULL)
    all.module.table <- get.geneSet.nodeInfo(geneset = sig.modules,
        background = V(PFN)$name, prefix = NULL, abbrev.map = NULL,
        subset.pattern = NULL, remove.pattern = NULL)
    tbl <- combine.table(combine.table(node.topo, module.table),
        all.module.table)
    return(tbl)
}


################################

get.nodeTopoInfo <- function(PFN,modules,KD,label.pttrn = NULL,outputfname = NULL)
{
    
	#### designate module membership
	cat("Assigning module/KDA membership\n")
	moduleMatrix <- do.call(cbind,lapply(modules,function(x,y) {y %in% x},y = V(PFN)$name));
	
	mod.mat <- matrix("",ncol = ncol(moduleMatrix),nrow = nrow(moduleMatrix));mod.mat[which(moduleMatrix)] <- "YES";colnames(mod.mat) <- names(modules)
	#write.table(data.frame(probe = V(PFN)$name,as.data.frame(mod.mat)),file = "moduleMembership.txt",sep = "\t",row.names = F,col.names = T,quote = F)
	
	rownames(moduleMatrix) <- V(PFN)$name;
	module.vector <- apply(moduleMatrix,1,function(x,y) paste(y[x],collapse = ","),y = colnames(moduleMatrix))
	rm(moduleMatrix)

	#### designate KD membership
	
	KD.Matrix <- do.call(cbind,lapply(KD,function(x,y) {y %in% x},y = V(PFN)$name))
	rownames(KD.Matrix) <- V(PFN)$name;
	KD.vector <- apply(KD.Matrix,1,function(x,y) paste(y[x],collapse = ","),y = colnames(KD.Matrix))
	rm(KD.Matrix)

	#### get node degree and clustering coefficients
    cat("Calculating node topological properties\n")
	node.k <- igraph::degree(graph = PFN);node.k <- node.k[V(PFN)$name]
    node.str <- graph.strength(graph = PFN);node.str <- node.str[V(PFN)$name]
	node.wcf <- transitivity(graph = PFN,type = "weighted");
	#node.tcf <- transitivity(graph = PFN,type = "local");node.tcf <- node.tcf[V(PFN)$name]
	#node.btw <- betweenness(graph = PFN,directed = FALSE,weights = d.func(E(PFN)$weight));node.btw <- node.btw[V(PFN)$name];

	#### coerce into single node table
	if (is.null(label.pttrn))
	{
	 node.table <- data.frame(node = V(PFN)$name,module.membership = module.vector,KD.membership = KD.vector,node.degree = node.k,node.strength = node.str,node.weightedCF = node.wcf)#,node.betweenness = node.btw)
	}else{
	 node.label <- gsub(label.pttrn,"",V(PFN)$name)
	 hub.only.label <- rep("",length(node.label));hub.only.label[which(KD.vector != "")] <- node.label[which(KD.vector != "")]
	 node.table <- data.frame(node = V(PFN)$name,node.label = node.label,hub.label = hub.only.label,module.membership = module.vector,KD.membership = KD.vector,node.degree = node.k,node.strength = node.str,node.weightedCF = node.wcf)#,node.betweenness = node.btw)
	}

	if (!is.null(outputfname))
	{
	 #outputfname <- gsub("\\.txt","\\.nodeTopoInfo\\.txt",net.file)
	 cat(paste("output:",outputfname,"\n",sep = ""))
   	 write.table(node.table,file = outputfname,sep = "\t",row.names = F,col.names = T,quote = F)
	}
	
   	return(node.table)
}

######## general function to extract 
collect.table.features <- function(table.files,cols)
{
 ### table.files = files containing gene features. each file is named to label specific column. gene id is placed on the first column in each table. 
 ### cols = column names to collect 
 tbl.lst <- lapply(table.files,function(x,cc) {
                                        tbl <- read.delim(file = x,sep = "\t",header = T)
										tbl[,c(1,match(cc,colnames(tbl)))]
							   },cc = cols)
 
 for (i in 1:length(tbl.lst))
 {
  colnames(tbl.lst[[i]])[2:ncol(tbl.lst[[i]])] <- paste(names(table.files)[i],colnames(tbl.lst[[i]])[2:ncol(tbl.lst[[i]])],sep = "__")
 }
 
 coerce.manyTables(tbl.lst)
}

######### gene set process functions: outputs YES.
get.geneSet.nodeInfo <- function(geneset,background,prefix = NULL,abbrev.map = NULL,subset.pattern = NULL,remove.pattern = NULL)
{
	#### geneset = list object containing gene sets
	#### background = single character string of genes for background
	#### prefix = a character object to be used as prefix in the node table
	#### abbrev.map = a named character string to perform gsub(names(abbrev.map)[i],abbrev.map[i],names(geneset))
	#### subset.pattern = regular expression to identify subset of gene sets to coerce into gene table.
	#### remove.pattern = regular expression to remove in gene symbols 
	
	if (!is.null(remove.pattern))
	{
	 geneset <- lapply(geneset,function(x,y) unique(gsub(y,"",x)),y = remove.pattern)
	 background <- unique(gsub(remove.pattern,"",background))
	}
	
	if (!is.null(subset.pattern))
	{
	 geneset <- geneset[grep(subset.pattern,names(geneset))]
	}

	if (!is.null(abbrev.map))
	{
		for (i in 1:length(abbrev.map))
		{
		 names(geneset) <- gsub(names(abbrev.map)[i],abbrev.map[i],names(geneset))
		}
	}
	
	######### use abbreviations for long names

	comp.col <- names(geneset)
	gene.row <- background;

	node.matrix <- matrix("",nrow = length(gene.row),ncol = length(comp.col))
	for (i in 1:length(geneset))
	{
	 r <- which(gene.row %in% geneset[[i]])
	 node.matrix[r,i] <- "YES"
	}
	
	if (!is.null(prefix))
	{
	 colnames(node.matrix) <- paste(prefix,comp.col,sep = ".")
	}else{
	 colnames(node.matrix) <- comp.col
	}
	
	node.matrix <- data.frame(node = gene.row,as.data.frame(node.matrix))

	return(node.matrix)

}
