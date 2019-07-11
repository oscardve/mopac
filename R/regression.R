# MoPAC internal sgRNA regression functions
# version: 1.0
# date: May 2018
# description: performs weighted average of sgRNAs sorted by rank, and bagging average of bootstrapped genes
# author:  Oscar Villarreal
# affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
# contact: oscardvillarreal AT gmail.com

use.weights <- function(folded, uniform.weights, empirical.weights) {
  folded <- folded[order(folded$Gene,folded$sgRNA),]
  Nsgrna <- table(folded$Gene) #find number of sgrnas per gene
  Ngene <- table(Nsgrna) #how many genes have 1 sgRNA, how many have 2, how many have 3, etc.
  if(length(Ngene)>1) stop("Not all genes have the same amount of guides!")
  # Sort sgRNA values:
  sorted0 <- split(as.matrix(folded[,-which(names(folded)%in%c("sgRNA","Gene","Category"))]),folded$Gene)
  sorted <- as.data.frame(do.call(rbind,lapply(sorted0,function(x) get_sorted1(x,Nsgrna[1]))))
  # Split into ranks:
  xijk <- split(sorted,rep(1:Nsgrna[1],times=(nrow(sorted)/Nsgrna[1])))
  # Bind into a 3D tensor (put samples in 1st dimension, genes in 2nd and sgRNA in 3rd):
  tensor <- array(unlist(xijk), dim = c(nrow(xijk[[1]]), ncol(xijk[[1]]), length(xijk)))
  tensor <- aperm(tensor, c(2,1,3))
  dimnames(tensor)[[1]] <- names(folded)[-which(names(folded)%in%c("sgRNA","Gene","Category"))]
  dimnames(tensor)[[2]] <- unique(folded$Gene)
  # Scoring (averaged):
  if(empirical.weights==TRUE){
    if(Nsgrna[1]>10) {
      weights <- rowMeans(apply(matrix(NA,nrow=1000,ncol=Nsgrna[1]),1,function(x){
        x[sort(sample(1:Nsgrna[1],4))] <- c(0.2792423,0.3346092,0.2440308,0.1421177)
        return(x)}),na.rm=TRUE)
    } else if(Nsgrna[1]==10) { weights <- c(0.2792423,0.2978576,0.3038050,0.2996329,0.2864142,0.2658323,0.2401169,0.2093034,0.1763553,0.1421177)
    } else if(Nsgrna[1]==9) { weights <- c(0.2792423,0.2999200,0.3051250,0.2967544,0.2785473,0.2503667,0.2172760,0.1802932,0.1421177)
    } else if(Nsgrna[1]==8) { weights <- c(0.2792423,0.3028417,0.3059574,0.2919458,0.2646420,0.2278163,0.1852644,0.1421177)
    } else if(Nsgrna[1]==7) { weights <- c(0.2792423,0.3069282,0.3055136,0.2813066,0.2413848,0.1929482,0.1421177)
    } else if(Nsgrna[1]==6) { weights <- c(0.2792423,0.3124171,0.3020701,0.2610669,0.2033731,0.1421177)
    } else if(Nsgrna[1]==5) { weights <- c(0.2792423,0.3208447,0.2892724,0.2184558,0.1421177)
    } else if(Nsgrna[1]==4) { weights <- c(0.2792423,0.3346092,0.2440308,0.1421177)
    } else if(Nsgrna[1]==3) { weights <- c(0.4746435,0.3780770,0.1472794)
    } else if(Nsgrna[1]==2) { weights <- c(0.6851414,0.3148586)
    }
    weights <- weights/sum(weights)
  } else {
    weights <- rep(1/Nsgrna[1],Nsgrna[1])
  }
  scored <- as.data.frame(apply(tensor,1,function(x) x%*%weights))
  if(nrow(scored)==1) scored <- t(scored) #for single-sample datasets
  if(identical(rownames(scored),rownames(tensor))) scored <- t(scored) #for single-gene datasets (oct23)
  colnames(scored) <- rownames(tensor)
  scored <- cbind(unique(folded[,which(colnames(folded)%in%c("Gene","Category")),drop=F]),scored)
  rownames(scored) <- colnames(tensor)
  # Genes with size larger than the mode are aggregated through bagging average:
  bootstraps <- grep("_boot",scored$Gene) #find indices of bootstrapped genes
  scored$Gene <- gsub("_boot.*","",scored$Gene) #remove the keyword "_boot"
  bootstraps <- unique(scored$Gene[bootstraps])
  bagged <- list()
  remove <- lapply(bootstraps, function(x) {
    isGene <- scored$Gene==x
    bagged[[x]] <<- cbind(scored[isGene,which(colnames(scored)%in%c("Gene","Category")),drop=F][1,,drop=F],as.data.frame(t(colMeans(scored[isGene,-grep("Gene|Category",names(scored)),drop=F]))))
    return(which(isGene))
  })
  if(length(unlist(remove))>0) {
    scored <- scored[-unlist(remove),] #remove bootstraps
    scored <- rbind(scored, do.call(rbind,bagged)) #add the averaged bootstraps
    scored <- scored[order(scored$Gene),] #reorder according to gene
  }
  if(length(grep("_boot",rownames(scored)))>0) scored <- scored[-grep("_boot",rownames(scored)),,drop=F] #do not include bootstraps in the output
  return(scored)
}

get.weights <- function(folded,pctrls,nctrls) {
  folded <- folded[order(folded$Gene,folded$sgRNA),]
  Nsgrna <- table(folded$Gene) #find number of sgrnas per gene
  Ngene <- table(Nsgrna) #how many genes have 1 sgRNA, how many have 2, how many have 3, etc.
  if(length(Ngene)>1) stop("Not all genes have the same amount of guides!")
  # Sort sgRNA values:
  sorted0 <- split(as.matrix(folded[,-which(names(folded)%in%c("sgRNA","Gene","Category"))]),folded$Gene)
  sorted <- as.data.frame(do.call(rbind,lapply(sorted0,function(x) get_sorted1(x,Nsgrna[1]))))
  # Split into ranks:
  xijk <- split(sorted,rep(1:Nsgrna[1],times=(nrow(sorted)/Nsgrna[1])))
  # Bind into a 3D tensor (put samples in 1st dimension, genes in 2nd and sgRNA in 3rd):
  tensor <- array(unlist(xijk), dim = c(nrow(xijk[[1]]), ncol(xijk[[1]]), length(xijk)))
  tensor <- aperm(tensor, c(2,1,3))
  dimnames(tensor)[[1]] <- as.list(names(folded)[-which(names(folded)%in%c("sgRNA","Gene","Category"))])
  dimnames(tensor)[[2]] <- as.list(unique(folded$Gene))
  # Regression and scoring:
  solution <- do.optimize(xp=tensor[,colnames(tensor)%in%pctrls,,drop=F], xn=tensor[,colnames(tensor)%in%nctrls,,drop=F])
  wk <- solution[["Weights"]][nrow(solution[["Weights"]]),]
  # Output bootstrapped weights:
  # utils::write.table(solution[["Weights"]],paste0("~/paper_july12/gecko_mopac_filtered/",length(wk),".txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
  # if(length(wk)<4) utils::write.table(data.frame(Rank=1:length(wk),Weight=wk),paste0("sabatini_june7/",length(wk),"-",wk[1],".txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
  Ai <- solution[["A"]][nrow(solution[["A"]]),]; Bi <- solution[["B"]][nrow(solution[["B"]]),]
  variancep <- solution[["Variancep"]][length(solution[["Variancep"]])]
  variancen <- solution[["Variancen"]][length(solution[["Variancen"]])]
  variance <- solution[["Variance"]][length(solution[["Variance"]])]
  scored <- apply(tensor,1,function(x) x%*%wk)
  scored <- as.data.frame(t(apply(scored,1,function(x) x*t(Ai) + t(Bi))))
  if(nrow(scored)==1) scored <- t(scored) #for single-column datasets
  colnames(scored) <- rownames(tensor)
  scored <- cbind(unique(folded[,c("Gene","Category")]),scored)
  rownames(scored) <- colnames(tensor)
  # Genes with size larger than the mode are aggregated through bagging average:
  bootstraps <- grep("_boot",scored$Gene) #find indices of bootstrapped genes
  scored$Gene <- gsub("_boot.*","",scored$Gene) #remove the keyword "_boot"
  bootstraps <- unique(scored$Gene[bootstraps])
  bagged <- list()
  remove <- lapply(bootstraps, function(x) {
    isGene <- scored$Gene==x
    bagged[[x]] <<- cbind(scored[isGene,c("Gene","Category")][1,],as.data.frame(t(colMeans(scored[isGene,-grep("Gene|Category",names(scored)),drop=F]))))
    return(which(isGene))
  })
  if(length(unlist(remove))>0) {
    scored <- scored[-unlist(remove),] #remove bootstraps
    scored <- rbind(scored, do.call(rbind,bagged)) #add the averaged bootstraps
    scored <- scored[order(scored$Gene),] #reorder according to gene
  }
  if(length(grep("_boot",rownames(scored)))>0) scored <- scored[-grep("_boot",rownames(scored)),,drop=F] #do not include bootstraps in the output
  return(list(Scored=scored,W=wk,A=Ai,B=Bi,Variance=variance,Variancep=variancep,Variancen=variancen))
}
