#---------------------------------------------------------------------------------------------------------
#ModulePrincComps, which computes the first principal component of each module. 
#Also it provides the variance explained for the first 5 PCs of a module.
#It is based on the singular value decomposition.
#It takes as input datExpr  for which the rows are samples and the columns are genes.
#Here is how you would use it
#PC1=ModulePrinComps2[datExpr,color1]
#Then PC1[[1]] provides a data frame with the first PC of each module
#PC1[[2]] provides a data frame where each column contains the percentages of 
#the variance explained by the first 10 PCs by module.
#---------------------------------------------------------------------------------------------------------
# the output is a list and each element is a matrix whose number of columns equals the number of modules
# the first element of the list contains the 1st principal components of all modules, and 
# the second element of the list contains the 2nd principal components of all modules, etc...
# the last element is the expained variations of top "no.pcs" PCs in each modules
#---------------------------------------------------------------------------------------------------------

require(impute)

ModuleHubProfileSingle = function(datexpr, no_pcs=10) {

    ngenes = dim(datexpr)[1]

    print("PC by ModuleHubProfileSingle .................. ")

    corrlmatrix = cor( t(datexpr), use = "pairwise.complete.obs")
    corrlmatrix = abs(corrlmatrix)

    diag(corrlmatrix)<- 0
    kin <- apply(corrlmatrix,2,sum, na.rm=TRUE) 

    orderK  = order(-kin)

    # select no.pcs genes which are ebv
    step = as.integer(ngenes/no_pcs)
    if(step<1){step=1; }

    selIdx= seq(from=1,to=ngenes, by=step)

    # mimic values from SVD
    #
    v = t(datexpr[selIdx[1:no_pcs], ])
    d = rep(0, no_pcs)

    ret   = NULL
    ret$v = v
    ret$d = d

    return(ret)
}

ModulePrinComps = function(datexpr,modules, min_modulesize=10) {

  couleur = modules;
  datEXpr = t(datExpr)

  no.pcs=10

  allmean = mean(datexpr, na.rm=T)

  #modlevels= names(table(couleur)) #levels(factor(couleur))
  modlevels <- names(couleur)

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels

  for(i in c(1:length(modlevels)) ){
 
    #print(paste(i, modlevels[i]) )
    modulename    = modlevels[i]
    #restrict1= as.character(couleur)== modulename
	restrict1 = colnames(datexpr) %in% couleur[[i]]
    restrict1= ifelse(is.na(restrict1), FALSE, restrict1)

    if(modlevels[i]=="grey"){
       next
    }

    # in the following, rows are genes and columns are samples
    datModule=t(datexpr[, restrict1])
    dim(datModule)

    # check whether some samples have missing rate (missing value in column) >80%
    naDat      = is.na(datModule)
    nasBySample= apply(naDat, 2, sum)
    nacutoff   = 0.6*dim(datModule)[1]
    selSamples = nasBySample>=nacutoff
    if(sum(selSamples)>0){
        print("patch samples with missing rate >=0.8")
        naSampleIdx = c(1:length(nasBySample))[selSamples]
        for(x in naSampleIdx){
            #print(paste("Sample Idx=", x) )
            datModule[,x]=ifelse(is.na( datModule[,x]), allmean, datModule[,x])
        }
    }

 
    if(sum(restrict1)<min_modulesize){
       listPCs[[1] ][,i] = datModule[1,]
       next
    }

    imputed=impute.knn(as.matrix(datModule)) #$data for R version < 2.0
    datModule =imputed$data #=imputed for R < 2.9

    datModule=t(scale(t(datModule)))

    # use the hub's profuile as PC if svd doesn't converge
    #svd1=svd(datModule)
    svd1=tryCatch(svd(datModule), error=function(e) ModuleHubProfileSingle(datModule), no_pcs=no.pcs)

    mtitle=paste("PCs of ", modulename," module", sep="")

    no.samples = dim(datModule)[2]

    actualpcs=min(dim(svd1$v)[2], no.pcs)

    #cat(modulename, as.character(i), as.character(actualpcs), "\n")

    #explained variation (%)
    listPCs[[ no.pcs+1] ][,i]= (svd1$d[1:no.pcs])^2/sum(svd1$d^2)

    # this is the first principal component, the ith column of the j-th element in the list
    for (j in c(1:actualpcs) ){
         #print(j)        
         pcj=svd1$v[,j]

         # detect NAs
         jnas = is.na(pcj)
         if (sum(jnas)>0){
             break;
         }

         signhj=sign(sum(cor(pcj,  t(datModule))))

         if( !(signhj == 0 | is.na(signhj)) )  pcj=signhj* pcj
         listPCs[[j] ][,i] = pcj
    }

  }

  PC.out = listPCs
  
  PC.res <- data.frame()
  l <- length(PC.out)
  for (i in 1:(l-1))
  {
	df <- data.frame(type = rep(paste("PC",i,sep = ""),nrow(PC.out[[i]])),id = colnames(datExpr),
	as.data.frame(PC.out[[i]]))
	PC.res <- rbind.data.frame(PC.res,df)
	rownames(PC.out[[i]]) <- colnames(datExpr)
	rm(df)
  }

  PC.res <- rbind.data.frame(PC.res,
  data.frame(type = rep("variance.explained",nrow(PC.out[[l]])),id = paste("PC",1:nrow(PC.out[[l]]),sep = ""),as.data.frame(PC.out[[l]])))
  rownames(PC.out[[l]]) <- paste("PC",1:nrow(PC.out[[l]]),sep = "")

  names(PC.out) <- c(paste("PC",1:(l-1),sep = ""),"variance.explained")
  
  return(PC.res)
  
}


ModuleHubProfiles = function(datexpr, adjmatrix, couleur, min_modulesize=10, whichmodule=NULL) {

  no.pcs=10

  modlevels= names(table(couleur)) #levels(factor(couleur))

  #the first no.pcs of elements in the list are PCs and the last one is the variations
  listPCs = as.list(rep(NULL, no.pcs+1) )

  for (i in c(1:no.pcs) ){
    listPCs[[i]]     = data.frame(matrix(0,nrow=dim(datexpr)[[1]],ncol= length(modlevels))) 
    colnames(listPCs[[i]]) <- modlevels
  }

  listPCs[[ no.pcs+1] ]= data.frame(matrix(0,nrow=no.pcs,ncol= length(modlevels))) #varexplained
  colnames(listPCs[[ no.pcs+1] ]) <- modlevels


    # take the profile of the most connected gene in each module as PC
    #
    colorI = as.character(colcode.reduced)
    kin    = computeModuleLinks(adjmatrix, couleur)

    # find the hub in each module 
    #
    kinIdxed= cbind(kin, c(1:length(couleur)), couleur)
    orderK  = order(-kin)
    kinIdxed= kinIdxed[orderK, ]
    orderK  = order(kinIdxed[,3])
    kinIdxed= kinIdxed[orderK, ]
    
    hubIdx    = rep(0, length(modlevels) )
    for(z in c(1:length(modlevels)) ) {        
        isel      = modlevels[z] == kinIdxed[,3]
        ikinIdxed = kinIdxed[isel, ]

        # extract hubs' profiles
        #
        listPCs[[1] ][,z] = datexpr[,as.integer(ikinIdxed [1,2])]
        hubIdx[z] = ikinIdxed [1,2]
    }

    # extract hubs' profiles
    #
    #listPCs[[1]] = t(datexpr[,as.integer(hubIdx)])
    names(listPCs[[1]]) <- modlevels

    if(!is.null(whichmodule) ){
       wsel = modlevels==whichmodule
       widx = c(1:length(modlevels))[widx]
       return(listPCs[[1] ][,wsel])
    }

    return(listPCs)
}


#-------------------------------------------------------------------------
#Function: compute whithin-module number of connections for each gene
#Input:   
#         adjMatrix  ~ adjacency matrix 
#         colorcodeC ~ color codes (modules) for genes
#Output:  in-module number of connections of each gene
#-------------------------------------------------------------------------
computeModuleLinks = function(adjMatrix, colorcodeC, isAdjacency=TRUE, normalized=FALSE, usegreymodule=F)
{
   modnames= names( table(colorcodeC) )

   no.genes     = dim(adjMatrix)[1]
   links        = rep(0, no.genes)
   idxseq       = c(1:no.genes)   
   for (each in modnames ){
      if((usegreymodule==F) & (as.character(each)=="grey" ) ){
         next
      }

      whichmod    = each==colorcodeC
      module.size = sum(whichmod)

      if (module.size==1){
        next
      }

      modk <- apply(adjMatrix[whichmod,whichmod],2,sum, na.rm=TRUE) 
      if(!isAdjacency){
          modk <- (module.size -modk)
      }

      #normalize against the module size
      if(normalized==TRUE){
         modk = modk/module.size
      }

      #indices of genes in the current module
      idxmod = idxseq * whichmod
      
      #put the links's to the buffer      
      links[idxmod] = modk
   }
   links
}
