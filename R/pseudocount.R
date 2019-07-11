# MoPAC internal pseudocount optimization function
# version: 1.0
# date: May 2018
# description: performs a single optimization for each pseudocount to minimize the variance of control genes with most frequent size only
# author:  Oscar Villarreal
# affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
# contact: oscardvillarreal AT gmail.com

get.pseudocount.simple <- function(counts,initial,groups,pseudorange) {
  dataset <- counts[counts$Category%in%c("Essential","Nonessential"),] #keep only controls
  datasetp <- unique(dataset$Gene[dataset$Category=="Essential"])
  datasetn <- unique(dataset$Gene[dataset$Category=="Nonessential"])
  # Optimize the pseudocount by minimizing the variance:
  variance <- variancep <- variancen <- numeric(length(pseudorange))
  scaling <- shifting <- weights <- list()
  for(i in 1:length(pseudorange)) {
    # Normalize data:
    logged <- cbind(dataset[,c("sgRNA","Gene","Category")],apply(dataset[,-which(colnames(dataset)%in%c("sgRNA","Gene","Category"))],2,function(x)log2(x/stats::median(x)+pseudorange[i])))
    # Calculate fold changes:
    folded <- cbind(logged[,c("sgRNA","Gene","Category")],do.call(cbind,unname(lapply(groups,function(y)logged[,setdiff(y$Sample,initial),drop=F]-logged[,intersect(y$Sample,initial),drop=T]))))
    # Average sgRNAs per gene:
    gene_reps <- do.call(rbind,lapply(split(folded[,-which(colnames(folded)%in%c("sgRNA","Gene","Category")),drop=FALSE],folded$Gene),colMeans))
    # Scale:
    meanp <- apply(gene_reps[rownames(gene_reps)%in%datasetp,,drop=F],2,mean)
    meann <- apply(gene_reps[rownames(gene_reps)%in%datasetn,,drop=F],2,mean)
    scaled <- apply(gene_reps,2,function(y){
      y <- y-mean(y[rownames(gene_reps)%in%datasetn])
      y <- -y/mean(y[rownames(gene_reps)%in%datasetp])
      return(y)
    })
    # Compute variances:
    variancep[i] <- sum(apply(scaled[rownames(gene_reps)%in%datasetp,,drop=F],2,stats::var))
    variancen[i] <- sum(apply(scaled[rownames(gene_reps)%in%datasetn,,drop=F],2,stats::var))
    variance[i] <- variancep[i] + variancen[i]
    print(paste0("Pseudocount = ",pseudorange[i],". Variance = ",variance[i]))
  }
  optimal <- which.min(variance)
  Pseudocount <- pseudorange[optimal]
  Variance <- data.frame(Pseudocount=pseudorange,Variance=variance,Variancep=variancep,Variancen=variancen)
  return(list(Pseudocount=Pseudocount,Variance=Variance,Weights=NULL,Scaling=NULL,Shifting=NULL))
}

get.pseudocount <- function(counts,initial,groups,pseudorange,shinyF=NULL) {
  # Uniformize the number of guides per control genes (make them all the same size):
  dataset <- counts[counts$Category%in%c("Essential","Nonessential"),] #keep only controls
  Nguide <- table(dataset$Gene) #guides per gene
  Ngene <- table(Nguide) #genes per guide content
  mode <- as.numeric(names(Ngene)[which.max(Ngene)]) #most frequent guide content
  less <- names(Nguide)[which(Nguide<mode)] #small genes
  more <- names(Nguide)[which(Nguide>mode)] #large genes
  tmp <- dataset
  remove <- lapply(more, function(x) { #split large control genes into equally sized groups by adding _1 _2 _3 etc)
    isGene <- tmp$Gene==x
    tmp[isGene,"Gene"] <<- paste0(tmp[isGene,"Gene"],"_",rep(1:floor(sum(isGene)/mode),each=mode))
    keep <- floor(sum(isGene)/mode)*mode
    isGene[isGene==T] <- c(rep(F,keep),rep(T,sum(isGene)-keep))
    return(which(isGene))
  })
  if(length(unlist(remove))>0) tmp <- tmp[-unlist(remove),] #remove leftovers for each gene
  tmp <- tmp[!tmp$Gene%in%less,] #remove small genes
  pctrls <- unique(tmp$Gene[tmp$Category=="Essential"])
  nctrls <- unique(tmp$Gene[tmp$Category=="Nonessential"])
  # Optimize the pseudocount by minimizing the variance:
  variance <- variancep <- variancen <- numeric(length(pseudorange))
  scaling <- shifting <- weights <- list()
  for(i in 1:length(pseudorange)) {
    # Normalize data:
    logged <- cbind(tmp[,c("sgRNA","Gene","Category")],apply(tmp[,-which(colnames(tmp)%in%c("sgRNA","Gene","Category"))],2,function(x)log2(x/stats::median(x)+pseudorange[i])))
    # Calculate fold changes:
    folded <- cbind(logged[,c("sgRNA","Gene","Category")],do.call(cbind,unname(lapply(groups,function(y)
      logged[,setdiff(y$Sample,initial),drop=F]-logged[,intersect(y$Sample,initial),drop=T]))))
    # Find variance with optimized parameters:
    optimized <- get.weights(folded,pctrls,nctrls)
    variance[i] <- optimized$Variance
    variancep[i] <- optimized$Variancep
    variancen[i] <- optimized$Variancen
    print(paste0("Pseudocount = ",pseudorange[i],". Variance = ",variance[i]))
    if(!is.null(shinyF)) shinyF(0,"Quality Control",paste0("Pseudocount = ",pseudorange[i],". Variance = ",variance[i]))
  }
  optimal <- which.min(variance)
  optimal <- data.frame(Pseudocount=pseudorange[optimal],Variance=variance[optimal])
  Variance <- data.frame(Pseudocount=pseudorange,Combined=variance,Positive=variancep,Negative=variancen)
  return(list(Optimal=optimal,Variance=Variance))
}


