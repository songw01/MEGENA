############################# FET pipeline
make.Pairwise.Tables <- function(geneSets1,geneSets2,background)
{
 mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 d <- t(mem1) %*% mem2;
 b <- abs(t(mem1) %*% (mem2-1))
 c <- abs(t(mem1-1) %*% (mem2))
 a <- t(mem1-1) %*% (mem2-1);

 ij <- do.call(rbind,lapply(1:length(geneSets1),function(i,j) cbind(rep(i,length(j)),j),j = 1:length(geneSets2)))

 pairwise.tables <- lapply(1:nrow(ij),function(i,ij,a,b,c,d) as.table(matrix(c(a[ij[i,1],ij[i,2]],b[ij[i,1],ij[i,2]],c[ij[i,1],ij[i,2]],d[ij[i,1],ij[i,2]]),nrow = 2)),ij = ij,a = a,b = b,c = c,d = d)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}

do.FisherExactTest <- function(table.count,N_bg = NULL)
{
 if (is.null(N_bg)) N_bg = sum(rowSums(table.count))
 
 out <- fisher.test(x = table.count,or = 1,alternative = "greater")
 odds.ratio <- out$estimate
 p.value <- out$p.value;
 geneSet1.count <- rowSums(table.count)[2]
 geneSet2.count <- colSums(table.count)[2]
 expected.count <- geneSet1.count/N_bg * geneSet2.count
 overlap.count <- table.count[2,2];
 fold.change <- overlap.count/expected.count
 
 out <- c(N_bg,geneSet1.count,geneSet2.count,expected.count,overlap.count,fold.change,odds.ratio,p.value)
 names(out) <- c("Background","set1_size","set2_size","expected.overlap","actual.overlap","enrichment.foldchange","odds.ratio","FET_pvalue")
 return(out)
}

perform.AllPairs.FET <- function(geneSets1,geneSets2,background,adjust.FET.pvalue = T)
{
 do.multicore = F
 n.cores = NULL
 
 pairwise.tables <- make.Pairwise.Tables(geneSets1,geneSets2,background)
 
 output <- lapply(pairwise.tables,do.FisherExactTest)
 
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 return(output)
}

make.paired.Tables <- function(geneSets1,geneSets2,ij,background)
{
 
 #mem1 <- do.call(cbind,lapply(geneSets1,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))
 #mem2 <- do.call(cbind,lapply(geneSets2,function(x,y) {vec <- rep(0,length(y));vec[which(y %in% x)] <- 1;return(vec)},y = background))

 pairwise.tables <- mapply(FUN = function(s1,s2,z) {
                                         v1 <- rep(0,length(z));v1[which(z %in% s1)] <- 1; 
										 v2 <- rep(0,length(z));v2[which(z %in% s2)] <- 1;
										 # n11,n12,n21,n22
										 as.table(matrix(c(sum(abs(v1 - 1) * abs(v2 - 1)),sum(abs(v2 - 1) * v1),sum(abs(v1 - 1) * v2),sum(v1 * v2)),nrow = 2))
                                },s1 = geneSets1[ij[,1]],s2 = geneSets2[ij[,2]],MoreArgs = list(z = background),SIMPLIFY = FALSE)
 names(pairwise.tables) <- paste(names(geneSets1)[ij[,1]],"~vs~",names(geneSets2)[ij[,2]],sep = "")
 return(pairwise.tables)
}


perform.ijPairs.FET <- function(geneSets1,geneSets2,ij,background,adjust.FET.pvalue = T)
{
 pairwise.tables <- make.paired.Tables(geneSets1,geneSets2,ij,background)
 
 output <- lapply(pairwise.tables,do.FisherExactTest)
 
 # collect outputs
 output <- data.frame(set1_Name = gsub("~vs~(.*)$","",names(output)),set2_Name = gsub("^(.*)~vs~","",names(output)),as.data.frame(do.call(rbind,output)))
 int.element <- mapply(FUN = function(x,y) paste(intersect(x,y),collapse = ","),x = geneSets1[as.character(output[[1]])],y = geneSets2[as.character(output[[2]])],SIMPLIFY = F)
 int.element <- do.call(c,int.element)
 if (adjust.FET.pvalue) output <- cbind.data.frame(output,data.frame(corrected.FET.pvalue = p.adjust(output$FET_pvalue,"bonferroni"),intersects = int.element))
 
 output <- output[order(output$FET_pvalue),]
 return(output)
}

