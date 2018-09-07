do.MEGENA <- function(g,
do.hubAnalysis = TRUE,
mod.pval = 0.05,hub.pval = 0.05,remove.unsig = TRUE,
min.size = 10,max.size = 2500,
doPar = FALSE,num.cores = 12,n.perm = 100,singleton.size = 3,
save.output = FALSE, oneplus=TRUE)
{
	# g = igraph object of PFN
	# mod.pval = module p-value
	# hub.pval = hub p-value
	
	if (doPar & getDoParWorkers() == 1 & num.cores > 1)
	{
		cl <- makeCluster(num.cores)
		registerDoParallel(cl)
		# check how many workers are there
		cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
	}
	
	###### do clustering
	cat("Commence multiscale clustering (multi-threaded)....\n") 
	module.output <- nested.kmeans.all(g = g,single.size = singleton.size, num.cores=num.cores);
	if (save.output) save(module.output,file = "multiscale_clusters.RData")
	
	# get the list of unsplit parent modules
	parent.idx <- which(is.infinite(module.output$module.alpha) &  is.na(module.output$module.pvalue))
	split.parent <- intersect(parent.idx,unique(module.output$module.relation[,1]))
	if (length(split.parent) > 0) parent.idx <- setdiff(parent.idx,split.parent)
	
	# get the list of significant split children
	child.modules <- module.output$modules[which(module.output$module.pvalue <= mod.pval & !(is.infinite(module.output$module.alpha) &  is.na(module.output$module.pvalue)))];
	
	# get the list of final modules
	sig.modules <- c(module.output$modules[parent.idx],child.modules)
	sig.modules <- sig.modules[which(sapply(sig.modules,length) >= min.size & sapply(sig.modules,length) <= max.size)]
	
	if (save.output) output.geneSet.file(sig.modules,"multiscale_significant.modules.txt")
	
	if (length(which(module.output$module.pvalue < mod.pval)) > 0 & do.hubAnalysis)
	{
		###### perform hub analysis and scale clustering 
		cat("Commence MHA...\n")
		hub.output <- get.multiScale.hubs(module.output,g,
		module.degreeStat = NULL,doPar = doPar,n.core = num.cores,remove.unsig = remove.unsig,
		alpha.range = NULL,n.perm = n.perm,pval = mod.pval,padjust.method = "bonferroni",output.figures = save.output)
		if (save.output) save(hub.output,file = "multiscale_hubAnalysis.RData")
		 
		node.table <- node.summary(PFN = g,module.output = module.output,hub.output = hub.output,module.pval = mod.pval,hub.pval = hub.pval)
		if (save.output) write.table(node.table,file = "multiscale_nodeSummary.txt",sep = "\t",row.names = F,col.names = T,quote = F)
		
		output <- list(module.output = module.output,hub.output = hub.output,node.summary = node.table);
	}else{
		output <- list(module.output = module.output,hub.output = NA,node.summary = NA)
	}
	
	return(output)
}
