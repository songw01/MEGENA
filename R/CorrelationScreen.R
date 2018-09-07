#require(arules)
#require(infotheo)
#require(WGCNA)

count.FD <- function(rho,thresh.vec)
{
	 # combine threshold and correlation values for efficiency
	 n.rho <- length(rho)
	 label.vec <- c(rep(1,length(rho)),rep(0,length(thresh.vec)));# label correlation value and threshold value
	 rho <- c(rho,thresh.vec)# combine correlation and threhold values for efficiency.
	 # sort correlation values
	 i <- order(rho)
	 rho <- rho[i]
	 label.vec <- label.vec[i]
	 
	 # count permuted correlation values above thresholds
	 j <- which(label.vec == 0)
	 output <- cbind(rho[j],(n.rho - cumsum(label.vec)[j])/n.rho)
	 colnames(output) <- c("rho.cutoff","FDR")
	 return(output)
} 

calculate.rho <- function(datExpr,n.perm,FDR.cutoff,estimator = "pearson",
                          rho.thresh = NULL,sort.el = TRUE, oneplus=TRUE,k=5,k_iter_max=10)
{
	if (is.null(rownames(datExpr))) rownames(datExpr) <- paste("g",1:nrow(datExpr),sep = "")
	gid <- rownames(datExpr)
	datExpr <- t(datExpr)
	if (estimator %in% c("pearson","spearman")){
		rho <- cor(datExpr,method = estimator)
	} else if (estimator %in% "mutualinformation"){
		datExpr_discrete = arules::discretizeDF(data.frame(datExpr), default=list("method"="cluster",
		                                                                          "centers"=k,
		                                                                          "iter.max"=k_iter_max))
		rho <- infotheo::mutinformation(datExpr_discrete)
	} else if (estimator %in% "bicor"){
		rho <- WGCNA::bicor(x=datExpr,quick=1)
	}
	if (oneplus){
		if (estimator %in% "mutualinformation"){
			message("warn:: Using 1+rho/2 with mutual information is unnecessary and not recommended.")
		}
		rho <- (1+rho)/2; ##IMPORTANT
	}
	
	if (is.null(rho.thresh)) rho.thresh <- seq(0,1,0.01)
	
	#### permute data matrix to calculate FDR
	nc <- nrow(datExpr)
	perm.ind <- lapply(1:n.perm,function(i,n) sample(1:n,n),n = nc)
	count.out <- vector("list",n.perm)
	for (i in 1:n.perm)
	{
		cat("i = ");cat(i);cat("\n");
		if (estimator %in% c("pearson","spearman")){
			random.rho <- cor(datExpr,datExpr[perm.ind[[i]],],method = estimator);
		} else if (estimator %in% "mutualinformation"){
			df1 = datExpr
			df2 = datExpr[perm.ind[[i]],]
			mi_df_2mat = matrix(data=NA,nrow=ncol(df1),ncol=ncol(df2))
			rownames(mi_df_2mat) = colnames(df1)
			colnames(mi_df_2mat) = colnames(df2)

			for(col1 in 1:length(colnames(df1))){
			  for(col2 in col1:length(colnames(df2))){
			    x_dis = arules::discretize(x=as.numeric(df1[,col1]), method="cluster", centers=k,iter.max=k_iter_max)
			    y_dis = arules::discretize(x=as.numeric(df2[,col2]), method="cluster", centers=k,iter.max=k_iter_max)
			    mi_df_2mat[col1,col2] = infotheo::mutinformation(x_dis,y_dis)
			    mi_df_2mat[col2,col1] = mi_df_2mat[col1,col2]
			  }
			}

			random.rho <- mi_df_2mat
		} else if (estimator %in% "bicor"){
			random.rho <- WGCNA::bicor(x=datExpr,y=datExpr[perm.ind[[i]],],quick=1)
		}
		if(oneplus){
			random.rho <- (1+random.rho)/2;
		}
		random.rho <- as.vector(random.rho[upper.tri(random.rho)]);
		
		count.out[[i]] <- count.FD(random.rho,rho.thresh)
		rm(random.rho)
	}
	PR = count.FD(as.vector(rho[upper.tri(rho)]),rho.thresh);PR = PR[,2]
	FPR = Reduce("+",lapply(count.out,function(x) x[,2]))/n.perm;FPR[1] <- 1;
	FDR = FPR/PR;FDR[which(FPR == 0)] <- 0;FDR[which(FDR > 1)] <- 1;
	FDR.table <- data.frame(cut.off = rho.thresh,FPR = FPR,PR = PR,FDR = FDR)

	# choose threshold respect to FDR.cutoff
	rho.cutoff <- min(FDR.table$cut.off[FDR.table$FDR < FDR.cutoff])

	# coerce into edgelist
	ij <- which(rho > rho.cutoff & upper.tri(rho),arr.ind = T)
	w <- apply(ij,1,function(xy,m) m[xy[1],xy[2]],m = rho)
	ijw <- cbind(ij,w);colnames(ijw) <- c("row","col","rho")

	ijw <- as.data.frame(ijw)
	ijw[[1]] <- gid[ijw[[1]]];ijw[[2]] <- gid[ijw[[2]]]

	if (sort.el) ijw <- ijw[order(ijw[,3],decreasing = T),]

	output <- list(signif.ijw = ijw,FDR = FDR.table)

	return(output)
}


calculate.rho.twoMat <- function(data.mat1,data.mat2,n.perm,FDR.cutoff,estimator = "pearson",
                                 rho.thresh = NULL,sort.el = TRUE, oneplus=TRUE,k=5,k_iter_max=10)
{
	if (ncol(data.mat1) != ncol(data.mat2)) stop("columns of data.mat1 and data.mat2 do not agree.")
	if (is.null(rownames(data.mat1))) rownames(data.mat1) <- paste("A",1:nrow(data.mat1),sep = "")
	if (is.null(rownames(data.mat2))) rownames(data.mat2) <- paste("B",1:nrow(data.mat2),sep = "")
	gid1 <- rownames(data.mat1);gid2 <- rownames(data.mat2)
	
	data.mat1 <- cbind(t(data.mat1))
	data.mat2 <- cbind(t(data.mat2))
	
	if (estimator %in% c("pearson","spearman")){
		rho <- cor(data.mat1,data.mat2,method = estimator);
	} else if (estimator %in% "mutualinformation"){
		df1 = data.mat1
		df2 = data.mat2
		mi_df_2mat = matrix(data=NA,nrow=ncol(df1),ncol=ncol(df2))
		rownames(mi_df_2mat) = colnames(df1)
		colnames(mi_df_2mat) = colnames(df2)

		for(col1 in 1:length(colnames(df1))){
		  for(col2 in col1:length(colnames(df2))){
		    x_dis = arules::discretize(x=as.numeric(df1[,col1]), method="cluster", centers=k,iter.max=k_iter_max)
		    y_dis = arules::discretize(x=as.numeric(df2[,col2]), method="cluster", centers=k,iter.max=k_iter_max)
		    mi_df_2mat[col1,col2] = infotheo::mutinformation(x_dis,y_dis)
		    mi_df_2mat[col2,col1] = mi_df_2mat[col1,col2]
		  }
		}

		rho <- mi_df_2mat
	} else if (estimator %in% "bicor"){
		rho <- WGCNA::bicor(x=data.mat1,y=data.mat2,quick=1)
	}
	if(oneplus){
		if (estimator %in% "mutualinformation"){
			message("warn:: Using 1+rho/2 with mutual information is unnecessary and not recommended.")
		}
		rho <- (1+rho)/2;
	}

	if (is.null(rho.thresh)) rho.thresh <- seq(0,1,0.01)
	
	#### permute data matrix to calculate FDR
	nc <- nrow(data.mat1)
	perm.ind <- lapply(1:n.perm,function(i,n) sample(1:n,n),n = nc)
	count.out <- vector("list",n.perm)
	for (i in 1:n.perm)
	{
		cat("i = ");cat(i);cat("\n");
		if (estimator %in% c("pearson","spearman")){
			random.rho <- cor(data.mat1,data.mat2[perm.ind[[i]],],method = estimator);
		} else if (estimator %in% "mutualinformation"){
			df1 = data.mat1
			df2 = data.mat2[perm.ind[[i]],]
			mi_df_2mat = matrix(data=NA,nrow=ncol(df1),ncol=ncol(df2))
			rownames(mi_df_2mat) = colnames(df1)
			colnames(mi_df_2mat) = colnames(df2)

			for(col1 in 1:length(colnames(df1))){
			  for(col2 in col1:length(colnames(df2))){
			    x_dis = arules::discretize(x=as.numeric(df1[,col1]), method="cluster", centers=k,iter.max=k_iter_max)
			    y_dis = arules::discretize(x=as.numeric(df2[,col2]), method="cluster", centers=k,iter.max=k_iter_max)
			    mi_df_2mat[col1,col2] = infotheo::mutinformation(x_dis,y_dis)
			    mi_df_2mat[col2,col1] = mi_df_2mat[col1,col2]
			  }
			}
			random.rho <- mi_df_2mat
		} else if (estimator %in% "bicor"){
			random.rho <- WGCNA::bicor(x=data.mat1,y=data.mat2[perm.ind[[i]],],quick=1)
		}
		if(oneplus){
			random.rho <- (1+random.rho)/2
		}
		random.rho <- as.vector(random.rho);
		
		count.out[[i]] <- count.FD(random.rho,rho.thresh)
		rm(random.rho)
	}
	PR = count.FD(as.vector(rho),rho.thresh);PR = PR[,2]
	FPR = Reduce("+",lapply(count.out,function(x) x[,2]))/n.perm;FPR[1] <- 1;
	FDR = FPR/PR;FDR[which(FPR == 0)] <- 0;FDR[which(FDR > 1)] <- 1;
	FDR.table <- data.frame(cut.off = rho.thresh,FPR = FPR,PR = PR,FDR = FDR)

	# choose threshold respect to FDR.cutoff
	rho.cutoff <- min(FDR.table$cut.off[FDR.table$FDR < FDR.cutoff])

	# coerce into edgelist
	ij <- which(rho > rho.cutoff & upper.tri(rho),arr.ind = T)
	w <- apply(ij,1,function(xy,m) m[xy[1],xy[2]],m = rho)
	ijw <- cbind(ij,w);colnames(ijw) <- c("row","col","rho")

	ijw <- as.data.frame(ijw)
	ijw[[1]] <- gid1[ijw[[1]]];ijw[[2]] <- gid2[ijw[[2]]]

	if (sort.el) ijw <- ijw[order(ijw[,3],decreasing = T),]

	output <- list(signif.ijw = ijw,FDR = FDR.table)

	return(output)
}

################################################

test.pairwiseCor <- function(data.mat1,data.mat2 = NULL,alternative = "two.sided",method = "pearson", 
                             oneplus=TRUE,k=5,k_iter_max=10)
{
 cat("##### Pairwise Correlation Analysis ######\n")
 
 if (is.null(data.mat2))
 {
  cor.pairs <- do.call(rbind,lapply(1:(nrow(data.mat1)-1),function(i,n) cbind(rep(i,(n-i)),(i+1):n),n = nrow(data.mat1)))
  data.mat2 <- data.mat1;
 }else{
  if (ncol(data.mat1) != ncol(data.mat2)) stop("data.mat1 and data.mat2 do not have same number of columns.");
  cor.pairs <- do.call(rbind,lapply(1:nrow(data.mat1),function(i,n) cbind(rep(i,n),1:n),n = nrow(data.mat2)))
 }
 
 cat(paste("method=",method,"\nN1=",nrow(data.mat1),"\nN2=",nrow(data.mat2),"\nTotal pairs=",nrow(cor.pairs),"\n",sep = ""))
 
 row.names <- rownames(data.mat1);
 col.names <- rownames(data.mat2);
 cat("Calculating correlation coefficient and respective p-value...\n")
 if (estimator %in% c("pearson","spearman")){
	 output <- apply(cor.pairs,1,function(ij,mat1,mat2,alternative,method) {
          out <- cor.test(x = mat1[ij[1],],y = mat2[ij[2],],alternative = alternative,method = method);
       	  out <- c(out$estimate,out$p.value);
       	  if(oneplus){
			  out$estimate = (1+out$estimate)/2;
			}
		  names(out) <- c("rho","p.value")
		  return(out)},
		  mat1 = data.mat1,mat2 = data.mat2,alternative = alternative,method = method)
 } else if (estimator %in% "bicor"){
 	 output <- apply(cor.pairs,1,function(ij,mat1,mat2,alternative,method) {
		  ## use WGCNA::bicor function as estimator here
	      out <- WGCNA::bicorAndPvalue(x = mat1[ij[1],],y = mat2[ij[2],]);
	   	  out <- list(estimate=out$bicor,p.value=out$p);
	   	  if(oneplus){
			  out$estimate = (1+out$estimate)/2;
			}
		  names(out) <- c("rho","p.value")
		  return(out)},
		  mat1 = data.mat1,mat2 = data.mat2,alternative = alternative,method = method)
 } else if (estimator %in% "mutualinformation"){
 	 output <- apply(cor.pairs,1,function(ij,mat1,mat2,alternative,method,k=k,k_iter_max=k_iter_max) {
			## assume bivariate normal for converting to p-values (first pass only)
			df1 = mat1[ij[1],]
			df2 = mat2[ij[2],]

			mi_df_2mat = matrix(data=NA,nrow=ncol(df1),ncol=ncol(df2))
			rownames(mi_df_2mat) = colnames(df1)
			colnames(mi_df_2mat) = colnames(df2)

			for(col1 in 1:length(colnames(df1))){
			  for(col2 in col1:length(colnames(df2))){
			    x_dis = discretize(x=as.numeric(df1[,col1]), method="cluster", centers=k,iter.max=k_iter_max)
			    y_dis = discretize(x=as.numeric(df2[,col2]), method="cluster", centers=k,iter.max=k_iter_max)
			    mi_df_2mat[col1,col2] = infotheo::mutinformation(x_dis,y_dis)
			    mi_df_2mat[col2,col1] = mi_df_2mat[col1,col2]
			  }
			}

			corrs = mi_df_2mat
			rho = sqrt(1-exp(1)**(-2*corrs)) #message("Calculating first-pass MI z-score difference using bivariate normal assumptions")
			df = min(ncol(df1),ncol(df2)) - 2
			pvals = matrix(2 * (1 - pt(abs(rho) * sqrt(df) /
			  sqrt(1 - rho^2), df)), nrow = nrow(rho))

			  out <- list(estimate=corrs,p.value=pvals);

			if(oneplus){
			  message("warn:: Using 1+rho/2 with mutual information is unnecessary and not recommended.")
			  out$estimate = (1+out$estimate)/2;
			}
			names(out) <- c("rho","p.value")
			return(out)},
			mat1 = data.mat1,mat2 = data.mat2,alternative = alternative,method = method)
 }
 if (nrow(output) != nrow(cor.pairs)) output <- t(output)
 output <- as.data.frame(output)
 output <- data.frame(row = cor.pairs[,1],col = cor.pairs[,2],output)
 
 output <- list(correlation = output,row.names = row.names,col.names = col.names);
 return(output)
}

test.pairwiseCor.par <- function(data.mat1,data.mat2 = NULL,n.cores,alternative = "two.sided",method = "pearson", 
                                 oneplus=TRUE,k=5,k_iter_max=10)
{
 cat("##### Pairwise Correlation Analysis ######\n")
 
 if (is.null(data.mat2))
 {
  cor.pairs <- do.call(rbind,lapply(1:(nrow(data.mat1)-1),function(i,n) cbind(rep(i,(n-i)),(i+1):n),n = nrow(data.mat1)))
  data.mat2 <- data.mat1;
 }else{
  if (ncol(data.mat1) != ncol(data.mat2)) stop("data.mat1 and data.mat2 do not have same number of columns.");
  cor.pairs <- do.call(rbind,lapply(1:nrow(data.mat1),function(i,n) cbind(rep(i,n),1:n),n = nrow(data.mat2)))
 }
 
 cat(paste("method=",method,"\nN1=",nrow(data.mat1),"\nN2=",nrow(data.mat2),"\nTotal pairs=",nrow(cor.pairs),"\n",sep = ""))
 
 row.names <- rownames(data.mat1);
 col.names <- rownames(data.mat2);
 cat("Calculating correlation coefficient and respective p-value...\n")
 # split pairs into n.cores
 dn <- ceiling(nrow(cor.pairs)/n.cores)
 fact <- factor(do.call(c,lapply(1:ceiling(nrow(cor.pairs)/dn),function(n,dn) rep(n,dn),dn = dn))[1:nrow(cor.pairs)])
 split.pairs <- lapply(split(1:nrow(cor.pairs),fact),function(ii,m) m[ii,],m = cor.pairs);rm(cor.pairs,fact,dn)
 
 output <- foreach(cpair = split.pairs) %dopar% {
 
                  out <- apply(cpair,1,function(ij,mat1,mat2,alternative,method,k=k,k_iter_max=k_iter_max) {
                  									if (method %in% c("pearson","spearman")){
	                                                    out <- cor.test(x = mat1[ij[1],],y = mat2[ij[2],],alternative = alternative,method = method);
								                    	out <- c(out$estimate,out$p.value);
								                    } else if (method %in% "bicor"){
								                    	out <- WGCNA::bicorAndPvalue(x = mat1[ij[1],],y = mat2[ij[2],]);
	   	  												out <- list(estimate=out$bicor,p.value=out$p);
	   	  											} else if (method %in% "mutualinformation"){
	   	  												## assume bivariate normal for converting to p-values (first pass only)
														df1 = mat1[ij[1],]
														df2 = mat2[ij[2],]

														mi_df_2mat = matrix(data=NA,nrow=ncol(df1),ncol=ncol(df2))
														rownames(mi_df_2mat) = colnames(df1)
														colnames(mi_df_2mat) = colnames(df2)

														for(col1 in 1:length(colnames(df1))){
														  for(col2 in col1:length(colnames(df2))){
														    x_dis = discretize(x=as.numeric(df1[,col1]), method="cluster", centers=k,iter.max=k_iter_max)
														    y_dis = discretize(x=as.numeric(df2[,col2]), method="cluster", centers=k,iter.max=k_iter_max)
														    mi_df_2mat[col1,col2] = infotheo::mutinformation(x_dis,y_dis)
														    mi_df_2mat[col2,col1] = mi_df_2mat[col1,col2]
														  }
														}

														corrs = mi_df_2mat
														rho = sqrt(1-exp(1)**(-2*corrs)) #message("Calculating first-pass MI z-score difference using bivariate normal assumptions")
														df = min(ncol(df1),ncol(df2)) - 2
														pvals = matrix(2 * (1 - pt(abs(rho) * sqrt(df) /
														  sqrt(1 - rho^2), df)), nrow = nrow(rho))

														  out <- list(estimate=corrs,p.value=pvals);
	   	  											}
							                    	if(oneplus){
							                    		if (method %in% "mutualinformation"){
							                    			message("warn:: Using 1+rho/2 with mutual information is unnecessary and not recommended.")
							                    		}
														out$estimate = (1+out$estimate)/2;
													}
													names(out) <- c("rho","p.value")
													return(out)},
													mat1 = data.mat1,mat2 = data.mat2,alternative = alternative,method = method)
				  if (nrow(out) != nrow(cpair)) out <- t(out)
				  out <- as.data.frame(out)
                  out <- data.frame(row = cpair[,1],col = cpair[,2],out)
				  return(out)
		   }

 output <- do.call("rbind.data.frame",output)
 output <- list(correlation = output,row.names = row.names,col.names = col.names);
 return(output)
}

calculate.correlation <- function(datExpr,doPerm = 100,doPar = FALSE,num.cores = 8,method = "pearson",
FDR.cutoff = 0.05,n.increment = 100,is.signed = FALSE,
output.permFDR = TRUE,output.corTable = TRUE,saveto = NULL, oneplus=TRUE)
{
 # Input
 # datExpr = expression matrix (row = probe,column = sample)
 # doPerm = number of permutations
 # num.cores = number of cores to call when parallel computing 
 # FDR.cutoff = FDR cut-off
 # n.increment = number of increments to segment (rho_min,rho_max).
 # is.signed = FALSE (unsigned correlation), TRUE (signed correlation)
 # output.permFDR = TRUE (output permutation indices into .txt file)
 # output.corTable = TRUE (output final correlation into .txt file)
 # saveto = designated save folder
 if(oneplus){
 	message("Using (1+r)/2 rho calculation")
 }
 if (doPerm == 0)
 {
  if (!doPar)
  {
   cor.output <- test.pairwiseCor(datExpr,method = method, oneplus=oneplus);
  }else{
   cor.output <- test.pairwiseCor.par(datExpr,n.cores = num.cores,method = method, oneplus=oneplus);
  }
  # extract significant correlation
  vertex.names <- cor.output$row.names
  edgelist <- cor.output[[1]];
  edgelist <- cbind(edgelist,fdr.q.value = p.adjust(edgelist$p.value,"fdr"));
  edgelist <- edgelist[which(edgelist$fdr.q.value < FDR.cutoff),]
  
  if (is.signed)
  {
   edgelist <- edgelist[order(edgelist[,3],decreasing = T),]
   edgelist <- data.frame(row = vertex.names[edgelist[[1]]],col = vertex.names[edgelist[[2]]],as.data.frame(edgelist[,3:ncol(edgelist)]));
  }else{
   sign <- rep("",nrow(edgelist));sign[which(edgelist[,3] < 0)] <- "negative";sign[which(edgelist[,3] > 0)] <- "positive";
   edgelist[,3] <- abs(edgelist[,3])
   edgelist <- edgelist[order(edgelist[,3],decreasing = T),]
   edgelist <- data.frame(row = vertex.names[edgelist[[1]]],col = vertex.names[edgelist[[2]]],as.data.frame(edgelist[,3:ncol(edgelist)]),sign = sign);
   
   # output results to files
   
  }
 }else{
  
  rho.output <- calculate.rho(datExpr,n.perm = doPerm,FDR.cutoff = FDR.cutoff,estimator = method,
                              rho.thresh =  seq(0,1,1/n.increment),sort.el = TRUE, oneplus=oneplus)
  
  if (output.permFDR)
  {
   if (!is.null(saveto))
   {
    write.table(rho.output$FDR,file = paste(saveto,"correlation_FDR_table.txt",sep = "/"),sep = "\t",row.names = F,col.names = T,quote = F)
   }else{
    write.table(rho.output$FDR,file = "correlation_FDR_table.txt",sep = "\t",row.names = F,col.names = T,quote = F)
   }
  }
 }
 
 if (output.corTable)
 {
  cat("- outputting correlation results...\n");
  if (!is.null(saveto))
  {
    write.table(rho.output$signif.ijw,file = paste(saveto,"Data_Correlation.txt",sep = "/"),sep = "\t",row.names = F,col.names = T,quote = F)
  }else{
    write.table(rho.output$signif.ijw,file = "Data_Correlation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
  }
 }
 edgelist <- rho.output$signif.ijw;
 
 return(edgelist)
}

##########################
