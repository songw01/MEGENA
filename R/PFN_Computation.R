#require(foreach);
#require(iterators);
#require(igraph)
#require(doParallel)

# set parallel backend
anyNA <- function(x) any(is.na(x))

###### deprecate set.parallel.backend() function in MEGENA__1.3.4-6 version. 
#set.parallel.backend <- function(num.cores = NULL)
#{
 
# max.cores <- as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))

# if (is.na(max.cores)) max.cores <- num.cores

# if (is.null(num.cores))
# {
#  num.cores <- max.cores;
# }else{
#  if (num.cores > max.cores)
  # {
   # cat("warning: num.cores exceeds maximum available in the system.")
   # num.cores <- max.cores
  # }
 # }
 
 # cl <- parallel::makeCluster(num.cores)

 # registerDoParallel(cl)
 
 # return(cl)
# }


serial.PFN <- function(sortedEdge,Ng,maxENum)
{
 iplanarityTesting <- function(epair,rows,cols,N)
{
 nrows <- c(rows,epair[1]);ncols <- c(cols,epair[2]);
 out <- planaritytest(as.integer(N),nrows,ncols)
 return(out);
}
 # initiate parameters 
 Ne <- maxENum;
 ne <- 0
 rows <- c();# Initiate row vector
 cols <- c();# initiate col vector
 weights <- c();#initiate weights vector

 # Check if the final PFG has less edges than the input PFG
 if (Ne < ne)
 {stop("The input PFG has exceeding number of edges.");}

 # Detect the starting point to start planarity checking

 ee <- 1;
 
 while (ne < Ne & ee <= nrow(sortedEdge))
 {
  if (iplanarityTesting(sortedEdge[ee,1:2],rows,cols,Ng))
  {
   rows <- c(rows,sortedEdge[ee,1])
   cols <- c(cols,sortedEdge[ee,2])
   weights <- c(weights,sortedEdge[ee,3])
   ne <- ne + 1;
  }
  ee <- ee + 1;
 }
 
 if (length(rows) == maxENum & planaritytest(as.integer(Ng),rows,cols))
 {print("PFG is complete.")}
 edgel <- cbind(rows,cols,weights);
 return(edgel);
}

# a function that purely computes PFN without any quality check on edgelist
compute.PFN.par <- function(sortedEdge,Ng,maxENum,Njob,Ncore,max.skipEdges = NULL,keep.track = TRUE,initial.links = NULL)
{

Ng <- as.integer(Ng)
  
iplanarityTesting <- function(epair,rows,cols,N)
{
 nrows <- c(rows,epair[1]);ncols <- c(cols,epair[2]);
 out <- planaritytest(as.integer(N),nrows,ncols)
 return(out);
}
  # Planarity testing function object - multiple edge
  iplanarityTestingMultiple <- function(epairs,rows,cols,N)
  {
   out <- vector("logical",length = dim(epairs)[1]);
   for (ei in 1:dim(epairs)[1])
   {
    nrows <- c(rows,epairs[ei,1]);ncols <- c(cols,epairs[ei,2]);
    out[ei] <- planaritytest(as.integer(N),nrows,ncols);
    rm(nrows,ncols);
   }
   return(out);
  }
  # Njob = number of jobs assigned per core.
  # Ncore = number of cores available.
  Ne <- maxENum;
  if (is.null(initial.links))
  {
    ne <- 0;
    rows <- c();# Initiate row vector
    cols <- c();# initiate col vector
    weights <- c();#initiate weights vector
  }else{cat("Using input links as initial links...\n");ne <- nrow(initial.links);rows <- initial.links[,1];cols <- initial.links[,2];weights <- initial.links[,3];}
  
  numedges <- nrow(sortedEdge) # total candidate edges

  # Detect the starting point to start planarity checking

  ee <- 0; # number of edges traversed

  de <- 0;

  multi <- 0.1;

  do.par <- FALSE;
  enum_qual <- Ng;

  if (is.null(max.skipEdges)) 
  {
	skip.turns <- ceiling((Ng * 10)/(Njob * Ncore))
  }else{
	skip.turns <- ceiling(max.skipEdges/(Njob * Ncore))
  }
  
  zero.hit <- 0;

  while (ne<Ne & ee<numedges & zero.hit < skip.turns)
  {
   # Set the number of edge subsets to be tested against planarity
   # evector <- (ee+1):(min(c((ee+Njobs*Ncore),numedges))); 
   # Break 
   enum_qual <- 0;  

   if (do.par)
   {
    evector <- vector("list",length = Ncore);
    nc = 1;
    while (nc <= Ncore & ee < numedges)
    {evector[[nc]] <- sortedEdge[(ee+1):min(c(numedges,ee+Njob)),];
     ee <- min(c(numedges,(ee+Njob)));nc = nc + 1;}

    evector <- evector[1:(nc-1)];    rm(nc);
    enum_all = sum(unlist(lapply(evector,function(x) {dim(x)[1]})))

    # Perform edge quality check by utilizing parallel computing.
    print(paste("Performing parallel quality checks on ",as.character(enum_all),sep = ""))

    # Subset of edges on which quality checks are performed
    outvector <- foreach(batch = 1:length(evector),.packages = "MEGENA") %dopar%
    {
     out <- iplanarityTestingMultiple(evector[[batch]],rows,cols,Ng);
     return(out);
    }
    
    enum_qual = sum(unlist(lapply(outvector,sum)));

    if (enum_qual == 0) {zero.hit = zero.hit + 1}else{zero.hit = 0;}

    print(paste("Qualified edges: ",as.character(enum_qual),sep = ""))

    # Get the list of qualified edges, and perform subsequent planarity testing.
    valid_evector <- matrix(0,nrow = 0,ncol = 3);
    for (batch in 1:length(evector))
    {
     elist <- evector[[batch]];
     elist <- elist[outvector[[batch]],];
     valid_evector <- rbind(valid_evector,elist);
     rm(elist);
    }

    ij = 0;
    while ((ij < dim(valid_evector)[1]) & (ne < Ne))
    {
     out <- iplanarityTesting(valid_evector[ij+1,],rows,cols,Ng)
     if (out)
     {rows <- c(rows,valid_evector[ij+1,1]);
     cols <- c(cols,valid_evector[ij+1,2]);
     weights <- c(weights,valid_evector[ij+1,3]);
     ne <- ne + 1;}
          
     ij = ij + 1;
    }
    rm(ij);
    print(paste(as.character(ne),"out of",as.character(Ne),"has been included."))
   
    rm(valid_evector);
   }else{
    evector <- sortedEdge[(ee+1):min(c(numedges,ee+Njob)),];
    enum_all = dim(evector)[1];
    # Perform edge-wise planarity test
    ij = 0;
    while ((ij < dim(evector)[1]) & (ne < Ne))
    {
     out <- iplanarityTesting(evector[ij+1,1:2],rows,cols,Ng)
     if (out)
     {
      rows <- c(rows,evector[ij+1,1]);
      cols <- c(cols,evector[ij+1,2]);
      weights <- c(weights,evector[ij+1,3]);
      ne <- ne + 1;
      enum_qual <- enum_qual + 1;
     }
     ij = ij + 1;
    }
    ee <- ee + dim(evector)[1];
    print(paste(as.character(ne),"out of",as.character(Ne),"has been included."))

    
    rm(ij);
   }
   
   # By comparing qualified edge numbers to the submitted edge list, determine
   # condition for parallel planarity check. 
   if (enum_qual < (enum_all/10))
   {do.par = TRUE}else{do.par = FALSE}
   
   if (keep.track) save(rows,cols,weights,file = "pfg_el.RData")
  } 
    
  if (length(rows) == Ne & planaritytest(Ng,rows,cols))
  {print("PFG is complete.")}
  edgel <- cbind(rows,cols,weights);
  return(edgel);
}

calculate.PFN <- function (edgelist, max.skipEdges = NULL,maxENum = NULL,doPar = FALSE, num.cores = NULL, keep.track = TRUE)
{
	if (is.null(num.cores)) num.cores = 1
    if (is.unsorted(rev(edgelist[[3]]))) edgelist <- edgelist[order(edgelist[[3]], decreasing = T),]
	if (is.null(max.skipEdges)) max.skipEdges = ceiling((num.cores * 1000) * 0.9999)
	
    vertex.names <- unique(c(as.character(unique(edgelist[[1]])),
        as.character(unique(edgelist[[2]]))))
    ijw <- cbind(match(edgelist[[1]], vertex.names), match(edgelist[[2]],
        vertex.names), edgelist[[3]])
    rm(edgelist)
    N <- length(vertex.names)
	
	if (is.null(maxENum)) maxENum = 3 * (N - 2)
    cat("####### PFN Calculation commences ########\n")
    if (!doPar) {
        PFN <- serial.PFN(sortedEdge = ijw, Ng = N, maxENum = maxENum)
    }
    else {
        PFN <- compute.PFN.par(sortedEdge = ijw, Ng = N, maxENum = maxENum, 
		Njob = 1000, Ncore = num.cores, max.skipEdges = max.skipEdges,
            keep.track = keep.track)
    }
    PFN <- data.frame(row = vertex.names[PFN[,1]], col = vertex.names[PFN[,2]], weight = PFN[, 3])
    return(PFN)
}


