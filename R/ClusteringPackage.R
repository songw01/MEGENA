########################### code chunks for random planar network generation and hub statistics
#require(igraph)
#require(Matrix)
#require(doParallel)
#require(fpc)

get.network.metric <- function(g,d.func = function(x) {1-x})
{
 dg <- g;igraph::E(dg)$weight <- d.func(igraph::E(dg)$weight);
 igraph::shortest.paths(dg)
}

####### gene set interface
read.geneSet <- function(geneSet.file)
{
 gene.list <- readLines(geneSet.file)
 gene.list <- strsplit(gene.list,"\t")
 
 names(gene.list) <- sapply(gene.list,function(x) x[1])
 gene.list <- lapply(gene.list,function(x) x[2:length(x)])
 return(gene.list)
}

output.geneSet.file <- function(geneSet,outputfname)
{
 if (!is.list(geneSet)) stop("geneSet is not a list.")
 if (is.null(names(geneSet))) stop("names(geneSet) is not defined properly.")
  
 sink(outputfname)
 cat(paste(paste(paste(names(geneSet),"\t",sapply(geneSet,function(x) paste(x,collapse = "\t")),sep = ""),collapse = "\n"),"\n",sep = ""))
 sink()

 return(0)
}

geneSet.to.Membership <- function(geneSet,genes = NULL)
{
 if (is.null(genes)) genes <- do.call(c,geneSet)
 if (any(table(genes) > 1)) stop("There are genes repeatedly appearing")
 if (is.null(names(geneSet))) names(geneSet) <- 1:length(geneSet)
 
 out <- rep("",length(genes));names(out) <- genes
 for (i in 1:length(geneSet))
 {
  out[which(names(out) %in% geneSet[[i]])] <- names(geneSet)[i]
 }
 return(out)
}

membership.to.geneSet <- function(membership)
{
 split(names(membership),factor(membership))
}

get.subnetworks <- function(PFN,modules)
{
 lapply(modules,function(x,PFN) igraph::induced.subgraph(g = PFN,vids = x),PFN = PFN)
}

######## network similarity packages

get.connectedModule <- function(genes,g)
{
 gg <- igraph::induced.subgraph(g = g,vids = genes);
 if (!igraph::is.connected(gg))
 {
  cls.out <- igraph::clusters(gg)
  out <- split(igraph::V(gg)$name,cls.out$membership);
 }else{
  out <- list(genes)
 }
 return(out)
}

get.LPI <- function(subNet,eta = 0.01,n.order = 3)
{
 if (n.order < 3) n.order <- 3;
 
 A <- igraph::get.adjacency(subNet)
 LP.index <- A %*% A
 for (n in 3:n.order)
 {
  LP.index <- LP.index + eta^(3-n) * (LP.index %*% A);
 }
 #LP.index <- as.matrix(LP.index)
 return(LP.index)
}

module.flowMatrix <- function(modules,SRW)
{
 if (is.matrix(SRW))
 {
  intraFlow <- sapply(modules,function(x,m) { 
                                       mm <- m[x,x];
									   if (length(x) > 1)
									   {
									    out <- sum(mm[upper.tri(mm)],na.rm = T)
									   }else{
									    out <- mm;
									   }
									   return(out)},m = SRW) 
 }else{
  intraFlow <- sapply(modules,function(x,m) { 
                                       mm <- m[x,x];
									   if (length(x) > 1)
									   {
									    out <- sum(Matrix::tril(mm),na.rm = T)
									   }else{
									    out <- mm;
									   }
									   return(out)},m = SRW) 
 }
 
 module.size <- sapply(modules,length)
 #moduleMatrix <- Matrix::Matrix(do.call(cbind,lapply(modules,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec);},y = rownames(SRW))))
 moduleMatrix <- do.call(cbind,lapply(modules,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec);},y = rownames(SRW)))
 
 flowMatrix <- as.matrix((t(moduleMatrix) %*% SRW) %*% moduleMatrix)
 diag(flowMatrix) <- intraFlow
  
 return(flowMatrix)
}

evaluate.flow <- function(modules,S)
{
 if (length(modules) > 1)
 {
  #S <- get.LPI(subNet)
  flow.mat <- module.flowMatrix(modules,S)
  norm.mat <- diag(1/diag(flow.mat)) %*% flow.mat
  sym.mat <- (norm.mat + t(norm.mat))/2
  out <- mean(sym.mat[upper.tri(sym.mat)],na.rm = T)
 }else{
  out <- NA;
 }
 return(out)
}

####### DBI
get.moduleMedian <- function(module,S) 
{
 if (length(module) > 1)
 {
  vec <- colMeans(S[module,module],na.rm = T)
  out <- module[which(vec == min(vec))[1]]
 }else{
  out <- module;
 }
 return(out)
}

get.intraModule.distance <- function(module,S)
{
 if (length(module) > 1)
 {
  subset.S <- S[module,module]
  median.gene <- get.moduleMedian(module,S)
  out <- sum(subset.S[median.gene,])/(length(module)-1)
 }else{
  out <- 0;
 }
 return(out)
}

get.DB.index <- function(modules,D,g)
{
 if (length(modules) > 1)
 {
  centroids <- sapply(modules,get.moduleMedian,S = D)
  intraDist <- sapply(modules,get.intraModule.distance,S = D)
  interDist <- D[centroids,centroids];
  for (i in 1:length(modules))
  {
   interDist[i,i] <- NA;
  }
 
  intraDistMat <- matrix(rep(intraDist,length(modules)),nrow = length(modules))
  intraDistMat <- intraDistMat + t(intraDistMat);
 
  out <- mean(apply(intraDistMat/interDist,1,function(x) max(x,na.rm = T)))
 }else{
  out <- Inf;
 }
 return(out)
}


####### nested kmeans pipeline

evaluate.boundary <- function(adj,memMatrix,S)
{
 dj <- adj;
 #dj[which(dj != 0)] <- 1;
 bndMatrix <- Matrix::t(Matrix::t(memMatrix) %*% dj) * (memMatrix - 1)
 bnd.nodes <- rownames(bndMatrix)[which(apply(bndMatrix,1,function(x) length(which(x < 0))) >= 1)]
 
 new.memMatrix <- memMatrix;
 module.membership <- S[bnd.nodes,] %*% memMatrix
 new.memMatrix[bnd.nodes,] <- Matrix::Matrix(t(apply(module.membership,1,function(x) {vec <- rep(0,length(x));vec[which(x == max(x))[1]] <- 1;return(vec)})))
 #new.memMatrix[bnd.nodes,] <- t(apply(module.membership,1,function(x) {vec <- rep(0,length(x));vec[which(x == max(x))[1]] <- 1;return(vec)}))

 
 return(new.memMatrix)
}

iterative.boundaryUpdate <- function(g,modules,S,max.iter = 10)
{
 adj <- igraph::get.adjacency(g,sparse = TRUE,attr = NULL)
 memMatrix <- Matrix::Matrix(0,nrow = vcount(g),ncol = length(modules))
 for (i in 1:length(modules)) memMatrix[match(modules[[i]],V(g)$name),i] <- 1
 rownames(memMatrix) <- igraph::V(g)$name
 
 do.update <- TRUE
 n.iter <- 0;
 while (do.update)
 {
  new.memMatrix <- evaluate.boundary(adj,memMatrix,S)
  comp.mat <- Matrix::t(memMatrix) %*% new.memMatrix;
  
  if (all((Matrix::diag(comp.mat) - Matrix::colSums(memMatrix)) == 0)) do.update = FALSE
  memMatrix <- new.memMatrix;
  n.iter <- n.iter + 1;
  if (n.iter < max.iter) do.update = FALSE
 }
 
 out <- apply(memMatrix,2,function(x,y) y[which(x != 0)],y = rownames(memMatrix))
 out <- out[which(sapply(out,length) > 0)]
 return(out)
}

get.median <- function(module,S)
{
 if (length(module) > 1)
 {
  x <- colSums(S[module,module],na.rm = T);
  out <- module[which(x == max(x,na.rm = T))][1];
 }else{
  out <- module
 }
 return(out)
}


iplanarityTesting <- function(epair,rows,cols,N)
{
 nrows <- c(rows,epair[1]);ncols <- c(cols,epair[2]);
 out <- planaritytest(as.integer(N),nrows,ncols)
 return(out);
}

#### random triangulation network package
random.T2 <- function(nv)
{
 el <- t(matrix(c(1,2,1,3,2,3),ncol = 3))
 tri <- rbind(c(1,2,3))
 N <- 3;
 
 if (nv < 4) nv = 4;
 
 while (N < nv)
 {
  N <- N+1
  tri.i <- sample(1:nrow(tri),1)
  new.el <- cbind(tri[tri.i,],rep(N,3))
  new.tri <- cbind(t(combn(tri[tri.i,],2)),rep(N,3))
  
  el <- rbind(el,new.el)
  if (N == 4)
  {
   tri <- rbind(tri,new.tri)
  }else{
   tri <- rbind(tri[-tri.i,],new.tri)
  }
 }
 
 output <- list(edgelist = el,surf.triangle = tri)
 
 return(output)
}

random.T1 <- function(el,surf.tri,n1)
{
 # make sure el is sequential
 v <- sort(Reduce("union",apply(el,2,unique)))
 clean.el <- matrix(match(el,v),ncol = 2)
 clean.surf.tri <- matrix(match(surf.tri,v),ncol = 3)
 
 N <- 0;
 surf.mat <- Matrix::Matrix(apply(clean.surf.tri,1,function(x) {null.vec <- rep(0,length(v));null.vec[x] <- 1;return(null.vec)}))
 adj.mat <- Matrix::Matrix(0,nrow = length(v),ncol = length(v))
 
 for (i in 1:nrow(clean.el))
 {
  adj.mat[clean.el[i,1],clean.el[i,2]] <- 1;
  adj.mat[clean.el[i,2],clean.el[i,1]] <- 1;
 }
 
 while (N < n1)
 {
  ei <- sample(1:nrow(clean.el),1)
  ij <- clean.el[ei,]

  adj.tri <- which(colSums(surf.mat[ij,]) == 2)
  
  if (length(adj.tri) == 2)
  {
   # get new edge and new triangle after T1 move
   cell.v <- which(rowSums(surf.mat[,adj.tri]) > 0)
   
   new.ij <- setdiff(cell.v,ij)
   
   if (adj.mat[new.ij[1],new.ij[2]] == 0)
   {
    # update edge and triangle
	clean.el <- rbind(clean.el[-ei,],new.ij) # edge update
	
	# triangle update
    surf.mat[,adj.tri[1]] <- 0
    surf.mat[c(new.ij,ij[1]),adj.tri[1]] <- 1;
	surf.mat[,adj.tri[2]] <- 0
    surf.mat[c(new.ij,ij[2]),adj.tri[2]] <- 1;
		
    # adjacency update
	adj.mat[ij[1],ij[2]] <- 0;adj.mat[ij[2],ij[1]] <- 0
	adj.mat[new.ij[1],new.ij[2]] <- 1;adj.mat[new.ij[2],new.ij[1]] <- 1;
	rm(new.ij)
	N <- N+1
   }
   
  }else{
   print(adj.tri)
   stop("error with surf.tri: two adjacent triangles must be recognized.")
  }
  rm(ei,ij,adj.tri)
  
  #if ((floor(N/n1 * 100) %% 10) == 0) {cat(paste("T1 complete:",ceiling(N/n1 * 100),"\\% done\n",sep = ""))}
 }
 
 final.edge <- cbind(v[clean.el[,1]],v[clean.el[,2]]);
 final.tri <- apply(surf.mat,2,function(x,v) v[which(x == 1)],v = v)
 output <- list(edgelist = final.edge,surf.triangle = final.tri)
 
 return(output)
 
}

match.connected <- function(el,ne)
{
 do.test = ne > 2 & (nrow(el) > ne)
 if (do.test)
 {
	 dn <- nrow(el) - ne;
	 nt <- 0
	 tel <- rbind(el);
	 
	 while (nt < dn)
	 {
	  nel <- rbind(tel[-sample(1:nrow(tel),1),])
	  if (igraph::is.connected(igraph::graph.edgelist(nel,directed = F)))
	  {
	   tel <- nel;
	   nt <- nt + 1;
	   rm(nel)
	  }
	 }
 }else{
  tel <- el
 }
 return(tel)
}

get.randomPlanar <- function(n,n1 = NULL,ne = NULL,verbose = TRUE)
{
 #if (is.null(n1)) n1 <- 3 * n;
 if (verbose) cat("Performing T2 moves...\n")
 T2.output <- random.T2(n)
 
 if (!is.null(n1))
 {
  if (verbose) cat("Performing T1 moves...\n")
  T1.output <- random.T1(el = T2.output[[1]],surf.tri = T2.output[[2]],n1 = n1)
  output.el <- T1.output$edgelist
 }else{
  T1.output <- list(NULL)
  output.el <- T2.output[[1]];
 }
 
 
 if (!is.null(ne))
 {
  #output.el <- output.el[sample(1:nrow(output.el),ne),]
  output.el <- match.connected(output.el,ne)
 }
 
 output <- list(T2.network = T2.output[[1]],T1.network = T1.output[[1]],edge.match = output.el)
 
 return(output)
}
globalVariables(names = c("ind","gene","is.hub","cpair"))
get.DegreeHubStatistic <- function(subnetwork,n.perm = 100,doPar = FALSE,n.core = 4)
{
 v <- V(subnetwork)$name;
  
 if (!doPar)
 {
  cat("permutation no.:")
  k.random <- c()
  for (n in 1:n.perm)
  {
     cat(n);cat(",")
 	 random.out <- get.randomPlanar(n = vcount(subnetwork),n1 = NULL,ne = ecount(subnetwork),verbose = FALSE)
	 random.net <- graph.edgelist(random.out[[3]],directed = F)
	 k.random <- c(k.random,degree(random.net))
  }
  cat("\n")
 }else{
  cat("permutation no.:")
  dn <- ceiling(n.perm/n.core);
  split.fact <- do.call(c,lapply(1:ceiling(n.perm/dn),function(n,dn) rep(n,dn),dn = dn))
  split.fact <- factor(split.fact[1:n.perm])
  split.ind <- split(1:n.perm,split.fact)

  k.random <- foreach (ind = split.ind, .combine = 'c') %dopar%
   {
    k.rand <- c()
    for (n in ind)
	{
	 cat(n);cat(",")
 	 random.out <- get.randomPlanar(n = igraph::vcount(subnetwork),n1 = NULL,ne = igraph::ecount(subnetwork),verbose = FALSE)
	 random.net <- igraph::graph.edgelist(random.out[[3]],directed = F)
	 k.rand <- c(k.rand,igraph::degree(random.net));
	}
	return(k.rand)
   }
   cat("\n")
 }
 
 h.out <- hist(k.random,breaks = 1:(max(k.random)+1) - 0.5,plot = F)
 FDR <- 1 - cumsum(h.out$counts/length(k.random));
 FDR <- data.frame(h.out$mids,FDR)
 k.cutoff <- min(FDR[[1]][FDR[[2]] < 0.01])
 
 k <- degree(subnetwork)
 pval <- vector("numeric",length(k))
 for (i in 1:length(k))
 {
  j <- which(FDR[[1]] >= k[i])
  if (length(j) > 0)
  {
   pval[i] <- FDR[[2]][j[1]]
  }else{
   pval[i] <- 0
  }
 }

 hub.df <- data.frame(gene = V(subnetwork)$name,degree = k,pvalue = pval,pvalue.BH = p.adjust(pval,"BH"))
 
 return(hub.df)
}

######################################################################################

####################################### code chunks for updated efficient clustering

########## multiscale-clustering in a single run
compact.indiv <- function(v,D,alpha)
{
 if (length(v) > 1)
 {
  d <- D[match(v,rownames(D)),match(v,colnames(D))]
  SPD.mu <- mean(d[upper.tri(d)],na.rm = T)
  out <- SPD.mu/(log(length(v))^alpha)
 }else{
  out <- Inf;
 }
 return(out)
}

get.random.D <- function(nv,w,d.func,name.vec)
{
 T2.output <- random.T2(nv);
 #T2.output <- get.randomPlanar(nv,n1 = NULL,ne = length(w))
 el <- cbind(T2.output[[1]],sample(w,nrow(T2.output[[1]]),replace = T));colnames(el) <- c("row","col","weight")
 g.rand <- graph.data.frame(as.data.frame(el),directed = F);
 #V(g.rand)$name <- name.vec
 D.rand <- get.network.metric(g.rand,d.func);
 return(D.rand)
}

evaluate.significance <- function(score.o,nv,D.rand,alpha,n.perm = 100)
{
 comp.perm <- rep(NA,n.perm)
 sample.ind <- lapply(1:n.perm,function(i,nm,nv) sample(nm,nv),nm = rownames(D.rand),nv = nv)
 
 perm.p <- compact.indiv(1:nrow(D.rand),D.rand,alpha)
 comp.perm <- sapply(sample.ind,function(v,d,a) compact.indiv(v,d,a),d = D.rand,a = alpha)
 comp.perm.rel <- comp.perm/perm.p
 
 FDR <- sum(comp.perm.rel <= score.o,na.rm = T)/length(comp.perm.rel)
 
 pvalue.norm <- pnorm(-abs((score.o - mean(comp.perm.rel,na.rm = T))/sd(comp.perm.rel,na.rm = T)))
 output <- list(parent.comp.rel = score.o,pval = FDR,pval.norm = pvalue.norm,permuted.data = comp.perm.rel)
 return(output)
}

######### clustering main functions

pam.split <- function(g,D,k.max = Inf,n.skip = 20)
{
 if (k.max > (igraph::vcount(g)-1)) k.max <- igraph::vcount(g)-1

 if (k.max >= 2)
 {
  S <- get.LPI(g,eta = 0.01,n.order = 3)
  D.dist <- as.dist(D[match(igraph::V(g)$name,rownames(D)),match(igraph::V(g)$name,colnames(D))]);# make sure D is specific for this network, g.
    
  k <- 2;
  n.sk <- 0;
  qual.o <- Inf;
  optimal.module <- NULL
  val <- c();
  module.history <- list()
  cat("k=");
  while (k <= k.max & n.sk <= n.skip)
  {
   cat(paste(k,",",sep = ""))
   k.out <- cluster::pam(x = D.dist,k = k,diss = TRUE,do.swap = FALSE,cluster.only = TRUE)
   new.modules <- split(names(k.out),factor(k.out))
   new.modules <- iterative.boundaryUpdate(g,new.modules,S)
   #qual.i <- evaluate.compactQuality(new.modules,D);
   if (length(new.modules) > 1)
   {
    qual.i <- evaluate.flow(new.modules,S);
   }else{qual.i <- Inf}
   #qual.i <- get.DB.index(modules = new.modules,D = D,g = g)
   val <- c(val,qual.i)
   if (qual.i < qual.o)
   {
    qual.o <- qual.i;
	optimal.module <- new.modules;
	n.sk <- 0;
   }else{
    n.sk <- n.sk + 1;
   }
   k <- k + 1;
   module.history <- c(module.history,list(new.modules))
   rm(new.modules)
  }
  cat("\n")
  names(module.history) <- paste("k = ",2:(length(val)+1),sep = "")
  output <- list(optimal.module = optimal.module,eval.stat = data.frame(k = 2:(length(val)+1),flow.val = val),module.history = module.history)
 }else{
  output <- list(optimal.module = list(igraph::V(g)$name),eval.stat = NULL)
 }
 return(output)
}

pam.cleanSplit <- function(g,D,k.max = 10,n.skip)
{
 if (igraph::vcount(g) > 1)
 {
  #SimMat <- get.LPI(g)
  out <- pam.split(g,D,k.max = k.max,n.skip = n.skip)
  modules <- out[[1]]
  #if (length(out[[1]]) > 1)
  #{
  # modules <- iterative.boundaryUpdate(g,out[[1]],SimMat)
  #}else{
  # modules <- out[[1]];
  #}
 }else{
  modules <- list(igraph::V(g)$name);
 }
 return(modules)
}




nested.kmeans <- function(sub.g,
k.max = Inf,n.skip = 20,n.perm = 100,module.pvalue = 0.05,sig.method = c("pval.perm","pval.norm"),
d.func = function(x) {1-x},n.singleton = 3,alpha.range = seq(0.01,10,0.01))
{
 if (!is.connected(sub.g)) {stop("input network is not connected.")}
 
 cat("Calculating distance metric and similarity...\n")
 D <- get.network.metric(sub.g,d.func = d.func);
 S <- get.LPI(sub.g);
  
 ##### while-loop to construct 
 
 # define modules.keep: modules that no longer evaluated as conditions are met. 
 modules.keep <- list(V(sub.g)$name)
 modules.singleton <- list()
 do.test <- 1;
 do.update <- TRUE
 
 alpha.defined <- c(Inf)
 pvalue.calc <- c(NA)
 comp.score <- c(NA)
 #modules.parent <- list()
 #parent.index <- c()
 subset.relation <- matrix(0,ncol = 2,nrow = 0)
  
 n.iter <- 1;
 while (do.update)
 {
  #test.update <- list()
  test.i <- which(do.test == 1)
  cat(paste("iteration:",n.iter,"\n",sep = ""))
  cat(paste("- #. tested:",length(test.i),"\n",sep = ""))
  
  for (i in 1:length(test.i))
  {
   
   do.test[test.i[i]] <- 0;
   alph <- alpha.range[which(alpha.range <= alpha.defined[test.i[i]])]
   
   D.rand <- get.random.D(nv = length(modules.keep[[test.i[i]]]),w = E(sub.g)[modules.keep[[test.i[i]]] %--% modules.keep[[test.i[i]]]]$weight,d.func = d.func,name.vec = V(sub.g)$name);
   #print(sub.g)
   comp.p <- compact.indiv(modules.keep[[test.i[i]]],D,alpha = alpha.range) # parent 
   
   sub.net <- igraph::induced.subgraph(sub.g,modules.keep[[test.i[i]]])
   cat("- ")
   pam.out <- pam.cleanSplit(sub.net,D,k.max,n.skip);
   pam.out <- do.call(c,lapply(pam.out,get.connectedModule,sub.net));
   # remove singletons for downsream analysis
   i.single <- which(sapply(pam.out,length) <= n.singleton)
   modules.singleton <- c(modules.singleton,pam.out[i.single]);
   pam.out <- pam.out[setdiff(1:length(pam.out),i.single)]
   
   cat(paste("- #. of split:",length(pam.out),"\n",sep = ""))
   
   if (length(pam.out) > 1)
   {
    cat("- assess improvements over compactness\n")
    # evaluate exceeding point
    comp.c <- lapply(pam.out,function(v,d,alph) compact.indiv(v,d,alph),d = D,alph = alpha.range)
    comp.rel <- lapply(comp.c,function(x,y) x/y,y = comp.p)
	
	# check if the split is valid: is any sub-module more compact at some alpha? 
	def.alpha <- rep(NA,length(comp.rel))
	comp.sig <- rep(1,length(comp.rel))
	comp.sc <- rep(NA,length(comp.rel))
	
	
    for (j in 1:length(comp.rel))
    {
	 if (length(comp.rel[[j]]) > 0)
	 {
		 if (any(comp.rel[[j]] < 1))
		 {
		  alpha.i <- max(which(comp.rel[[j]] < 1))
		  def.alpha[j] <- min(c(alpha.range[alpha.i],alpha.defined[test.i[i]]))
		  stat.out <- evaluate.significance(score.o = comp.rel[[j]][alpha.i],nv = length(pam.out[[j]]),D.rand = D.rand,alpha = alph[alpha.i],n.perm = n.perm)
		  if (sig.method == "pval.perm") {pval = stat.out$pval}else{pval = stat.out$pval.norm}
		  comp.sig[j] <- pval
		  comp.sc[j] <- comp.rel[[j]][alpha.i]
		 }
	 }
    }
	
	subset.relation <- rbind(subset.relation,cbind(rep(test.i[i],length(pam.out)),(length(modules.keep)+1):(length(modules.keep) + length(pam.out))))
	modules.keep <- c(modules.keep,pam.out)
	
	
	if (all(is.na(def.alpha) & comp.sig >= module.pvalue))
	{
	 # if there is no meaningful split in terms of alpha and pvalue,
	 # then there is no need to accept the split.
	 do.test <- c(do.test,rep(0,length(pam.out)))
	}else{
	 update.test <- rep(0,length(pam.out))
	 update.test[which(!is.na(def.alpha) & comp.sig < module.pvalue)] <- 1;
	 do.test <- c(do.test,update.test)
	 def.alpha[is.na(def.alpha)] <- alpha.defined[test.i[i]];
	}
	alpha.defined <- c(alpha.defined,def.alpha)
	pvalue.calc <- c(pvalue.calc,comp.sig);
	comp.score <- c(comp.score,comp.sc)
   }
   
  }
  if (all(do.test == 0)) do.update <- FALSE
  n.iter <- n.iter + 1
 } 
 
 alpha.defined[is.na(alpha.defined)] <- 0;
 names(modules.keep) <- 1:length(modules.keep);
 output <- list(modules = modules.keep,singletons = modules.singleton,module.relation = subset.relation,module.alpha = alpha.defined,module.pvalue = pvalue.calc,module.compRel = comp.score)
 return(output)
}


####### adaptor function to process many connected components at once
nested.kmeans.all <- function(g,single.size = 3,sig.method = "pval.perm",
k.max = Inf,n.skip = 20,n.perm = 100,module.pvalue = 0.05,
d.func = function(x) {1-x},alpha.range = seq(0.01,10,0.01))
{
	con.comp <- get.connectedModule(V(g)$name,g)

	output.mod <- vector("list",length(con.comp))
	for (i in 1:length(con.comp))
	{
	 gg <- igraph::induced.subgraph(g,con.comp[[i]])
	 output.mod[[i]] <- nested.kmeans(sub.g = gg,sig.method = sig.method,n.singleton = single.size,k.max = k.max,n.skip = n.skip,
	 module.pvalue = module.pvalue,d.func = d.func,alpha.range = alpha.range)
	}

	###### combine outputs from all components

	singletons <- list()
	modules <- list()
	module.relation <- matrix(0,nrow = 0,ncol = 2)
	module.alpha <- c()
	module.pvalue <- c()
	module.compRel <- c()
	for (i in 1:length(output.mod))
	{
	 new.mod <- output.mod[[i]]$modules
	 names(new.mod) <- paste("c",i,"_",1:length(new.mod),sep = "")
	 modules <- c(modules,new.mod)
	 
	 mod.rel <- output.mod[[i]]$module.relation;
	 mod.rel <- cbind(names(new.mod)[mod.rel[,1]],names(new.mod)[mod.rel[,2]])
	 module.relation <- rbind(module.relation,mod.rel)
	 
	 singletons <- c(singletons,output.mod[[i]]$singletons)
	 module.alpha <- c(module.alpha,output.mod[[i]]$module.alpha)
	 module.pvalue <- c(module.pvalue,output.mod[[i]]$module.pvalue)
	 module.compRel <- c(module.compRel,output.mod[[i]]$module.compRel)
	}

	j <- which(sapply(modules,length) > single.size)
	modules <- modules[j]
	module.alpha <- module.alpha[j]
	module.pvalue <- module.pvalue[j]
	module.compRel <- module.compRel[j]
	names(singletons) <- NULL

	integer.relation <- cbind(match(module.relation[,1],names(modules)),match(module.relation[,2],names(modules)))

	output <- list(modules = modules,singletons = singletons,module.relation = integer.relation,module.alpha = module.alpha,
	module.pvalue = module.pvalue,module.compRel = module.compRel)

	return(output)
}

##########


get.union.cut <- function(module.output,alpha.cut,output.plot = T,plotfname = "validModules_alpha",module.pval = 0.05,remove.unsig = T)
{
 ####### obtain the partition based on the rule: if any sub-cluster is realized significant, accept the split. 
 # this condition enforces the following sub-conditions: 
 # i) a parent cluster must be defined at alpha.cut, 
 # ii) if a split exists, at least one sub-cluster is significant,
 # iii) if no valid split exists, then parent itself is
 
 gh <- graph.data.frame(as.data.frame(module.output$module.relation[,c(2,1)]),directed = T);
 gh <- gh + vertex(setdiff(1:length(module.output$modules),unique(sort(module.output$module.relation[]))))
 is.valid <- module.output$module.alpha[as.integer(V(gh)$name)] >= alpha.cut;
 names(is.valid) <- V(gh)$name
 
 module.o <- V(gh)$name[igraph::degree(gh,mode = "out") == 0]
 module.keep <- c()
 
 while (length(module.o) > 0)
 {
  update.o <- c();
  for (i in 1:length(module.o))
  {
   child <- setdiff(V(gh)$name[neighborhood(graph = gh,nodes = module.o[i],order = 1,mode = "in")[[1]]],module.o[i])
  
   # first check: if any child exists
   if (length(child) > 1)
   {
    # second check: if any child is significant
    if (any(is.valid[child]))
    {
     update.o <- c(update.o,child[is.valid[child]])
	 module.keep <- c(module.keep,child[!is.valid[child]])
    }else{
	 module.keep <- c(module.keep,module.o[i]);
    }
   }else{
    module.keep <- c(module.keep,module.o[i]);
   }
  }
  module.o <- update.o;
  rm(update.o)
 }
 
 if (!remove.unsig)
 {
  cut.module <- module.output$modules[as.integer(module.keep)]
 }else{
  cut.module <- module.output$modules[intersect(as.integer(module.keep),which(module.output$module.pvalue < module.pval))]
 }
 
 if (output.plot)
 {
  vert.size <- log10(sapply(module.output$modules,length)[as.integer(V(gh)$name)]);vert.size <- vert.size/max(vert.size) * 7
  col.vec <- rep("white",vcount(gh));col.vec[V(gh)$name %in% module.keep] <- "red";col.vec[as.integer(V(gh)$name) %in% which(module.output$module.pvalue >= module.pval)] <- "grey"
  xy <- layout.kamada.kawai(gh)
  
  if (!is.null(plotfname))
  {
   outputfname <- paste(plotfname,"__",alpha.cut,".png",sep = "")
   png(outputfname,width = 1000,height = 1000);
  }
  plot(gh,vertex.color = col.vec,vertex.label = V(gh)$name,vertex.label.color = "black",vertex.size = vert.size,edge.arrow.size = 0.6,vertex.label.cex = 1.2,layout = xy,edge.color = "black");
  if (!is.null(plotfname))
  {
   dev.off()
  }
 }
 
 return(cut.module)
 
}


######## multi scale hub evaluation

### multiscale degree profile
get.DegreeProfiles <- function(module.output,network,module.pval = 0.05,remove.unsig = TRUE,alpha.range = NULL)
{
	if (is.null(alpha.range)) alpha.range <- sort(unique(setdiff(module.output$module.alpha,Inf)))
	cat("Calculating within-module degree profiles.....\n")
	all.genes <- Reduce("union",module.output$modules)
	connect.matrix <- matrix(0,nrow = length(all.genes),ncol = length(alpha.range));
	rownames(connect.matrix) <- all.genes;colnames(connect.matrix) <- alpha.range;
	for (i in 1:length(alpha.range))
	{
	 cut.output <- get.union.cut(module.output,alpha.cut = alpha.range[i],output.plot = F,plotfname = "validModules_alpha",module.pval = module.pval,remove.unsig = remove.unsig)
	 module.str <- lapply(cut.output,function(x,g) {gg <- induced.subgraph(g,x);igraph::graph.strength(gg)},g = network) 
	 names(module.str) <- NULL
	 module.str <- do.call(c,module.str);
	 connect.matrix[,i] <- module.str[all.genes]
	}
	connect.matrix[is.na(connect.matrix)] <- 0;

	return(connect.matrix)
}

### cluster similar scales

cluster.scales <- function(connect.matrix,K.max = NULL,save.output = FALSE)
{
    #if (is.null(K.max)) K.max = ncol(connect.matrix)-1;
	
	restrict.kmax = is.matrix(connect.matrix) & is.null(K.max)
	if (restrict.kmax) K.max <- min(c((ncol(connect.matrix)-1),K.max))
	
	if (K.max > 20) K.max = 20; # restrict cluster scales to be 20 at max.
	
	vec <- rep(1,ncol(connect.matrix));names(vec) <- colnames(connect.matrix)
	output <- list(clusters = vec);
	cat(paste("K.max:",K.max,"\n",sep = ""))
	if (K.max > 1)
	{
		#require(fpc);require(cluster)
		cat("Cluster scales based on degree profiles...\n")
		D <- dist(t(connect.matrix),"euclidean")
		k.cls <- list()
		cat("k = ")
		for (k in 2:K.max)
		{
		 cat(k);cat(",")
		 k.cls <- c(k.cls,list(pam(x = D,k = k,cluster.only = TRUE)))
		}
		cat("\n")

		#k.cls <- lapply(k.output,function(x) x$cluster)
		
		stat.out <- lapply(k.cls,function(cls,d) cluster.stats(d = d,clustering = cls),d = D)
		stat.extract <- c("avg.silwidth","pearsongamma","dunn","sindex")
		#stat.extract <- c("avg.silwidth","pearsongamma","dunn")
		stat.out <- do.call(rbind,lapply(stat.out,function(x,y) unlist(x[y]),y = stat.extract))
		stat.out <- rbind(stat.out)
		
		if (nrow(stat.out) > 1)
		{
			if (K.max > 2)
			{
				majority.vote <- rowMeans(apply(stat.out,2,function(x) {
															   #out <- log(x/max(x,na.rm = T))
															   out <- log(rank(x)/length(x))
															   #out <- rank(x)/length(x)
															   return(out)
															  }))
				names(majority.vote) <- 2:K.max
			
				#X11();
				if (save.output)
				{
					png("majority_vote.png",width = 600,height = 600)
					barplot(majority.vote,xlab = "k",ylab = "mean(log(normalized rank))")
					dev.off()
				}
				#X11();
				if (save.output)
				{
					png("cluster_stats.png",width = 900,height = 900)
					par(mfrow = c(2,2),mar = c(3,5,1,1))
					for (c in 1:ncol(stat.out))
					{
					 if (!all(is.nan(stat.out[,c]) | is.infinite(stat.out[,c]) | is.na(stat.out[,c]))) plot(2:K.max,stat.out[,c],xlab = "k",ylab = colnames(stat.out)[c])
					}
					dev.off()
				}
				optimal.i <- max(which(majority.vote == max(majority.vote,na.rm = T)))
				optimal.cls <- k.cls[[optimal.i]]
				
				if (save.output)
				{
					png("scalecluster_heatmap.png",width = 700,height = 700)
					heatmap(as.matrix(D),Colv = "Rowv",scale = "none",col = heat.colors(50),ColSideColors = as.character(optimal.cls))
					dev.off()
				}
			}else{
				majority.vote <- 0
				optimal.cls <- k.cls[[1]]
			}
			
		}else{
			majority.vote <- 0
			optimal.cls <- k.cls[[1]]
		}
				
		## correct for non-continuous alpha values
		optimal.cls <- optimal.cls[order(as.numeric(names(optimal.cls)))];
		unique.cls <- unique(optimal.cls);
		cls.id <- max(unique.cls) + 1;
		for (i in 1:length(unique.cls))
		{
		 j <- which(optimal.cls == unique.cls[i])
		 len <- length(j)
		 dj <- max(j) - min(j) + 1
		 
		 if (dj > len)
		 {
		  bin.vec <- optimal.cls == unique.cls[i]
		  pred <- FALSE;
		  for (j in 1:length(bin.vec))
		  {
		   if (!pred & bin.vec[j]) optimal.cls[j] <- cls.id;
		   if (pred & bin.vec[j]) optimal.cls[j] <- cls.id;
		   
		   if (pred & !bin.vec[j]) {cls.id = cls.id + 1}
		   pred <- bin.vec[j]
		  }
		 }
		}
		
		##### get centroid scales as representative scales
		
		#if (summary.method == "centroid")
		#{
		# summary.alpha <- sapply(split(names(optimal.cls),factor(optimal.cls)),function(alph,d) {val <- rowSums(cbind(d[alph,alph]));max(as.numeric(alph[which(val == min(val,na.rm = T))]))},d = as.matrix(D))
		#}
		#if (summary.method == "max")
		#{
		# summary.alpha <- sapply(split(names(optimal.cls),factor(optimal.cls)),function(x) max(as.numeric(x)))
		#}
		#if (summary.method == "min")
		#{
		# summary.alpha <- sapply(split(names(optimal.cls),factor(optimal.cls)),function(x) min(as.numeric(x)))
		#}
		#if (summary.method == "median")
		#{
		# summary.alpha <- sapply(split(names(optimal.cls),factor(optimal.cls)),function(x) median(as.numeric(x)))
		#}
		
		#X11();
		
		output <- list(clusters = optimal.cls,cluster.across.K = k.cls,quality.score = stat.out,summary.score = majority.vote)
	}
	return(output)
}

###### overall functions

get.multiScale.hubs <- function(module.output,network,
module.degreeStat = NULL,doPar = FALSE,n.core = 4,
alpha.range = NULL,n.perm = 100,pval = 0.05,remove.unsig = TRUE,
padjust.method = "bonferroni",use.pvalue = "perm.pvalue.correct",alpha.summary.method = "centroid",
output.figures = FALSE)
{
	##### calculate modulewise degree statistics
	if (is.null(module.degreeStat))
	{
		cat("Calculating hub significance.....\n")
		subnet.lst <- lapply(module.output$modules,function(x,g) induced.subgraph(g,x),g = network)
		
		#if (doPar & (getDoParWorkers() == 1)) set.parallel.backend(n.core)
		module.degreeStat <- vector("list",length(subnet.lst));names(module.degreeStat) <- names(subnet.lst)
		for (i in 1:length(subnet.lst))
		{
		 module.degreeStat[[i]] <- get.DegreeHubStatistic(subnet.lst[[i]],n.perm = n.perm,doPar = doPar,n.core = n.core)
		}
	}
	
	##### calculate loglikelihood of hubs
	if (is.null(alpha.range)) alpha.range <- sort(unique(module.output$module.alpha))
	

	all.genes <- Reduce("union",module.output$modules)
	loglik.matrix <- matrix(0,nrow = length(all.genes),ncol = length(alpha.range));
	rownames(loglik.matrix) <- all.genes;colnames(loglik.matrix) <- alpha.range;
	for (i in 1:length(alpha.range))
	{
	 if (!is.infinite(alpha.range[i]))
	 {
	  cut.output <- get.union.cut(module.output,alpha.cut = alpha.range[i],output.plot = F,plotfname = "validModules_alpha",module.pval = pval,remove.unsig = remove.unsig)
	 }else{
	  cut.output <- module.output$modules[which(is.infinite(module.output$module.alpha))]
	 }
	 sigModule.stat <- module.degreeStat[names(cut.output)]
	 gene.pvalue <- lapply(sigModule.stat,function(x) {vec <- -log(x$pvalue);names(vec) <- as.character(x[[1]]);return(vec)});names(gene.pvalue) <- NULL
	 gene.pvalue <- do.call(c,gene.pvalue);
	 gene.pvalue <- gene.pvalue[intersect(names(gene.pvalue),all.genes)]
	 loglik.matrix[match(names(gene.pvalue),all.genes),i] <- gene.pvalue
	}
	loglik.matrix[is.na(loglik.matrix)] <- 0;
	loglik.matrix[is.infinite(loglik.matrix)] <- max(setdiff(loglik.matrix,Inf),na.rm = T);

	##### identify scales that cluster together in first PCs
	cat("Identifying similar scales....\n- ")
	#connect.matrix <- get.DegreeProfiles(module.output,network,module.pval = pval,remove.unsig = remove.unsig,alpha.range = alpha.range)
	connect.matrix <- get.DegreeProfiles(module.output,network,module.pval = pval,remove.unsig = FALSE,alpha.range = alpha.range)
	scale.clusters <- cluster.scales(connect.matrix,K.max = NULL,save.output = output.figures)
	output <- NULL
	if (!is.null(scale.clusters))
	{
		cls <- scale.clusters$clusters;
		cls <- cls[match(names(cls),colnames(loglik.matrix))]
			
		cat(paste("- identified: ",length(unique(cls)),"\n",sep = ""))
		split.col <- split(1:ncol(loglik.matrix),factor(cls))
		split.alpha <- split(alpha.range,factor(cls))
		scale.alpha.summary <- data.frame(id = paste("S",order(sapply(split.alpha,min)),sep = ""),min.alpha = sapply(split.alpha,min),max.alpha = sapply(split.alpha,max))
		names(split.alpha) <- as.character(scale.alpha.summary$id)
		names(split.col) <- as.character(scale.alpha.summary$id)
		#names(split.alpha) <- sapply(split.alpha,function(x) paste(min(x),max(x),sep = "-"))
		#names(split.col) <- sapply(split.alpha,function(x) paste(min(x),max(x),sep = "-"))
		split.col <- split.col[order(scale.alpha.summary$min.alpha)]
		
		scaleCluster.modules <- vector("list",length(split.alpha));names(scaleCluster.modules) <- names(split.alpha)
		for (i in 1:length(split.alpha))
		{
		 name.lst <- lapply(split.alpha[[i]],function(alph,mod.out,pval,remove.unsig) names(get.union.cut(mod.out,alpha.cut = alph,output.plot = F,plotfname = "validModules_alpha",module.pval = pval,remove.unsig = remove.unsig)),
		 mod.out = module.output,pval = pval,remove.unsig = remove.unsig)
		 name.lst <- Reduce("union",name.lst)
		 scaleCluster.modules[[i]] <- module.output$modules[name.lst]
		}
		##### evaluate significant hubs in each clustered scale
		cat("Identifying hub genes significant in each scale level...\n")

		split.loglik <- lapply(split.col,function(i,m) cbind(m[,i]),m = loglik.matrix)
		hub.stat <- vector("list",length(split.loglik));names(hub.stat) <- names(split.col)
		#n.perm = 100
		for (i in 1:length(split.loglik))
		{
		 
		 loglik.m <- split.loglik[[i]]
		 real.stat <- 2 * rowSums(loglik.m,na.rm = T)
		 rand.stat <- c()
		 for (n in 1:n.perm)
		 {
		  rand.matrix <- matrix(sample(loglik.m,length(loglik.m)),ncol = ncol(loglik.m))
		  #rand.matrix <- cbind(apply(loglik.m,2,function(x) sample(x,length(x))))
		  rand.stat <- c(rand.stat,2 * rowSums(rand.matrix,na.rm = T))
		 }
		 
		 h.out <- hist(rand.stat,breaks = seq(0,max(c(max(rand.stat,na.rm = T),max(real.stat,na.rm = T))) * 1.1,0.01),plot = F)
		 prop <- h.out$counts/sum(h.out$counts)
		 inv.prop <- 1-cumsum(prop)
		 pvalue.table <- data.frame(x = h.out$mids,pvalue = inv.prop)
		 
		 stat.pvalue <- sapply(real.stat,function(x,tbl) tbl[[2]][min(which(tbl[[1]] > x))],tbl = pvalue.table)
		 fisher.pval <- pchisq(real.stat, df = 2*ncol(loglik.m), lower.tail = FALSE)
		 hub.stat[[i]] <- data.frame(gene = rownames(loglik.matrix),as.data.frame(loglik.m),combined.stat = real.stat,fisher.pvalue = fisher.pval,fisher.pvalue.correct = p.adjust(fisher.pval,padjust.method),perm.pvalue = stat.pvalue,perm.pvalue.correct = p.adjust(stat.pvalue,padjust.method))
		 rm(h.out,loglik.m,real.stat)
		}

		scalewise.hubs <- lapply(hub.stat,function(x,y) as.character(x[[1]][which(x[[which(names(x) == y)]] <= pval)]),y = use.pvalue)

		output <- list(hub.stat = hub.stat,hub.list = scalewise.hubs,loglik.matrix = loglik.matrix,degreeProfile.matrix = connect.matrix,scale.clustering = scale.clusters,
		module.degreeStat = module.degreeStat,scale.summary.clusters = scaleCluster.modules,scale.cluster.info = scale.alpha.summary)
		
		#### plot pvalue distribution
		#X11();
		if (output.figures)
		{
			png("pvalue_distributions.png",width = 800,height = 800)
			par(mfrow = c(2,ceiling(length(hub.stat)/2)))
			for (i in 1:length(hub.stat))
			{
			 plot(sort(hub.stat[[i]][[which(names(hub.stat[[i]]) == use.pvalue)]]),main = names(hub.stat)[i],ylab = use.pvalue,xlab = "")
			 abline(h = pval,col = "red")
			}
			dev.off()
			
			#### plot gene sets by tile plot
			
			plot.df <- data.frame()
			for (i in 1:length(scalewise.hubs))
			{
			 ddf <- data.frame(gene = scalewise.hubs[[i]],scale = rep(names(scalewise.hubs)[i],length(scalewise.hubs[[i]])),is.hub = rep(TRUE,length(scalewise.hubs[[i]])))
			 plot.df <- rbind.data.frame(plot.df,ddf);
			 rm(ddf)
			}
			gene.order <- names(sort(table(do.call(c,scalewise.hubs))))
			#X11();
			png("scalewise_hub_tileplot.png",width = 350,height = 1000)
			tile.obj <- ggplot(data = plot.df,aes(x = scale,y = gene,fill = is.hub)) + theme_bw() + scale_y_discrete(limits = gene.order) + scale_x_discrete(limits = names(scalewise.hubs)) + theme_bw() +
			theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 20),axis.text.y = element_text(size = 5)) + geom_tile();
			print(tile.obj)
			dev.off()
		}
	}
	return(output)
}

get.scalwiseHubs <- function(hub.stat,pval = 0.01,use.pvalue = "perm.pvalue.BH")
{
 lapply(hub.stat,function(x,y) as.character(x[[1]][which(x[[which(names(x) == y)]] <= pval)]),y = use.pvalue)
}
