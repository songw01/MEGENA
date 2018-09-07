module_convert_to_table <- function(MEGENA.output,mod.pval = 0.05,hub.pval = 0.05,min.size = 10,max.size)
{ 	
	summary.output <- MEGENA.ModuleSummary(MEGENA.output,mod.pvalue = mod.pval,hub.pvalue = hub.pval,
	min.size = min.size,max.size = max.size,annot.table = NULL,symbol.col = NULL,id.col = NULL,
	output.sig = TRUE)

	modules <- summary.output$modules
	#####
	is.annted <- any(sapply(modules,function(x) length(grep("\\|",x)) > 0))
	df <- mapply(FUN = function(x,y) data.frame(id = x,module = rep(y,length(x))),x = modules,y = as.list(names(modules)),SIMPLIFY = FALSE)
	df <- do.call('rbind.data.frame',df)
	if (is.annted)
	{
		df <- data.frame(id = df[[1]],gene.symbol = gsub("\\|(.*)","",as.character(df[[1]])),gene.symbol2 = gsub("^(.*)\\|","",as.character(df[[1]])),
		module = df[[2]])
	}
	
	# update module parent relationship
	if (!is.null(summary.output$module.table))
	{
		modtbl <- summary.output$module.table
		df <- cbind.data.frame(df[,-ncol(df)],module.parent = modtbl$module.parent[match(df$module,modtbl$module.id)],module = df[[ncol(df)]])
		rm(modtbl)
		colnames(df)[1] <- "id"
	}
	
	# update node statistics
	if (!is.null(MEGENA.output$node.summary))
	{
		colnames(df)[1] <- "id"
		# get hub summary
		hubs <- lapply(MEGENA.output$hub.output$module.degreeStat,function(x,hp) as.character(x[[1]])[which(x$pvalue < hp)],hp = hub.pval)
		hvec <- rep(NA,nrow(df))
		for (j in 1:length(hubs)) hvec[which(df[[1]] %in% hubs[[j]] & df$module == names(hubs)[j])] <- "hub"
		df <- cbind.data.frame(df[,-ncol(df)],node.degree = MEGENA.output$node.summary$node.degree[match(df$id,MEGENA.output$node.summary$id)],
		node.strength = MEGENA.output$node.summary$node.strength[match(df$id,MEGENA.output$node.summary$id)],is.hub = hvec,
		module = df[[ncol(df)]]);
		rm(MEGENA.output,summary.output)
	}
	rownames(df) <- NULL
	return(df)
}


module_rank=function(X){
#input:
# X is a matrix of discriminant values measuring the correlations between modules (in rows) and variables/features (in columns)
#output:
# a vector of module ranking scores
	X=as.matrix(X)
	G=apply(X,2,function(y) {z=rank(y);M=max(z);(M+1-z)/M} )
	apply(G,1,function(x) sum(log(x)))
}

flatten.table <- function(master.table,factor.col,out.cols)
{
 # assume that first column marks test names
 if (length(factor.col) == 1)
 {
  factor.vec <- factor(master.table[[factor.col]])
 }else{
  factor.vec <- factor(apply(do.call(cbind,lapply(master.table[factor.col],as.character)),1,function(x) paste(x,collapse = "~")))
 }
 split.table <- lapply(split(1:nrow(master.table),factor.vec),function(ii,m) m[ii,],m = master.table)

 common.id <- Reduce("union",lapply(split.table,function(x) as.character(x[[1]])))
 
 big.table <- list()
 for (out.col in out.cols)
 {
  aligned.table <- do.call(cbind,lapply(split.table,function(tbl,id,out.col) {
                                         vec <- rep(NA,length(id));names(vec) <- id;

										 if (is.numeric(tbl[[out.col]])) vec[as.character(tbl[[1]])] <- as.numeric(tbl[[out.col]])
										 if (is.factor(tbl[[out.col]])) vec[as.character(tbl[[1]])] <- as.character(tbl[[out.col]])
										 return(vec)
                                 },id = common.id,out.col = out.col))
  colnames(aligned.table) <- paste(names(master.table)[out.col],colnames(aligned.table),sep = "__")
  big.table <- c(big.table,list(aligned.table))
 }
 big.table <- data.frame(module = common.id,do.call(cbind.data.frame,big.table))
 return(big.table)
}

combine.table <- function(abl,bbl)
{
 common.id <- union(as.character(abl[[1]]),as.character(bbl[[1]]))
 abl.align <- do.call(cbind,lapply(abl[2:ncol(abl)],function(x,y,z) {
                                   names(x) <- z;
                                   vec <- rep(NA,length(y));names(vec) <- y;
								   if (is.numeric(x)) {vec[match(names(x),names(vec))] <- x;}else{
								   vec[match(names(x),names(vec))] <- as.character(x);}
								   
								   return(vec)
                          },y = common.id,z = as.character(abl[[1]])))
 bbl.align <- do.call(cbind,lapply(bbl[2:ncol(bbl)],function(x,y,z) {
                                   names(x) <- z;
                                   vec <- rep(NA,length(y));names(vec) <- y;
								   if (is.numeric(x)) {vec[match(names(x),names(vec))] <- x;}else{
								   vec[match(names(x),names(vec))] <- as.character(x);}
								   
								   return(vec)
                          },y = common.id,z = as.character(bbl[[1]])))
 out <- cbind.data.frame(data.frame(id = common.id),abl.align,bbl.align)
 return(out) 
}


coerce.manyTables <- function(table.lst)
{
 if (length(table.lst) > 1)
 {
  out.table <- table.lst[[1]]
  for (i in 2:length(table.lst))
  {
   out.table <- combine.table(out.table,table.lst[[i]])
  }
 }else{
  out.table <- table.lst[[1]];
 }
 
 return(out.table)
}

