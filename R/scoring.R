# MoPAC scoring internal function
# version: 1.0
# date: May 2018
# description: internal function used by the analyze.essentiality module to perform bootstrapping of genes
# author:  Oscar Villarreal
# affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
# contact: oscardvillarreal AT gmail.com

get.scores <- function(sgRNA,uniform.weights,empirical.weights,shinyF=NULL) {
  print("Bootstrapping large genes..."); if(!is.null(shinyF)) shinyF(0.1,"Essentiality analysis","Bootstrapping large genes...")
  tmp <- sgRNA
  # Control genes with size larger than the mode are split randomly into groups of size equal to the mode, removing the left-over sgRNAs:
  if(uniform.weights==FALSE & empirical.weights==FALSE){
    Nguide <- table(sgRNA$Gene[sgRNA$Category%in%c("Essential","Nonessential")]) #guides per control
    Ngene <- table(Nguide) #genes per guide content
    mode <- as.numeric(names(Ngene)[which.max(Ngene)]) #most frequent guide content
    bigCtrl <- names(Nguide)[which(Nguide>mode)]
    remove <- lapply(bigCtrl, function(x) { # add _1 _2 _3 etc to the large control names
      isGene <- tmp$Gene==x
      tmp[isGene,"Gene"] <<- paste0(tmp[isGene,"Gene"],"_ctrl",rep(1:floor(sum(isGene)/mode),each=mode))
      keep <- floor(sum(isGene)/mode)*mode
      isGene[isGene==T] <- c(rep(F,keep),rep(T,sum(isGene)-keep))
      return(which(isGene))
    })
    if(length(unlist(remove))>0) tmp <- tmp[-unlist(remove),] #remove leftovers for each control
  }
  # Group genes by size:
  Nguide <- table(tmp$Gene) #guides per gene
  Ngene <- table(Nguide) #genes per guide content
  # if(empirical.weights==FALSE){
  mode <- as.numeric(names(Ngene)[which.max(Ngene)]) #most frequent guide content
  # } else mode <- 4
  same <- names(Nguide)[which(Nguide==mode)]; more <- names(Nguide)[which(Nguide>mode)]; less <- names(Nguide)[which(Nguide<mode)]
  inSame <- tmp$Gene%in%same; inMore <- tmp$Gene%in%more; inLess <- tmp$Gene%in%less
  tmpSame <- tmp[inSame,]
  tmpMore <- tmp[inMore,]
  # Genes with size larger than the mode are processed by generating 10 bootstraps of themselves and adding them to the genes with size equal to the mode:
  print(paste0("Size = ",mode)); if(!is.null(shinyF)) shinyF(0,"Essentiality analysis",paste0("Size = ",mode))
  if(sum(inMore)>0) {
    tmpMore <- split(tmpMore,tmpMore$Gene)
    tmpMore <- do.call(rbind,lapply(tmpMore, function(x) {
      bootstraps <- replicate(10, sample(1:nrow(x),mode))
      return(do.call(rbind, apply(bootstraps,2, function(y) x[y,])))
    }))
    tmpMore$Gene <- paste0(tmpMore$Gene,"_boot",rep(1:10,each=mode)) #label each of the bootstraps
    tmp1 <- rbind(tmpSame,tmpMore) #merge genes whose size is equal to the mode with all bootstraps
  } else tmp1 <- tmpSame
  # Process genes with same or more sgRNAs as the mode:
  if(uniform.weights==FALSE & empirical.weights==FALSE) {
    isCtrl <- tmp$Category%in%c("Essential","Nonessential") & !tmp$Gene%in%less #small controls cannot be used
    pctrls <- unique(tmp$Gene[tmp$Category=="Essential" & isCtrl]) #useful positive controls
    nctrls <- unique(tmp$Gene[tmp$Category=="Nonessential" & isCtrl]) #useful negative controls
    out <- get.weights(tmp1,pctrls,nctrls)
    print(paste0("Variance = ",out$Variance)); if(!is.null(shinyF)) shinyF(0.1,"Essentiality analysis",paste0("Variance = ",out$Variance))
    scored1 <- out$Scored; wk <- out$W; Ai <- out$A; Bi <- out$B
    Ai <- data.frame(Sample=names(Ai),Score=Ai); Bi <- data.frame(Sample=names(Bi),Score=Bi)
    # Average scores of large splitted controls:
    largectrl <- scored1[grep("_ctrl",scored1$Gene),]
    if(nrow(largectrl)>0) {
      largectrl$Gene <- substring(largectrl$Gene,1,regexpr("_ctrl",largectrl$Gene)-1)
      largectrl1 <- do.call(rbind,lapply(split(largectrl[,-c(1:2),drop=F],largectrl$Gene),colMeans))
      largectrl1 <- data.frame(Gene=rownames(largectrl1),largectrl1,check.names=FALSE)
      largectrl1 <- merge(unique(largectrl[,c("Gene","Category")]),largectrl1)
      scored1 <- rbind(scored1[-grep("_ctrl",scored1$Gene),],largectrl1)
    }
  } else scored1 <- use.weights(tmp1,uniform.weights,empirical.weights)
  # Genes with smaller size than the mode are processed by randomly subsetting the controls 3 independent times:
  if(sum(inLess)>0) {
    if(uniform.weights==FALSE & empirical.weights==FALSE) {
      print("Bootstrapping controls for small genes...")
      if(!is.null(shinyF)) shinyF(0.1,"Essentiality analysis","Bootstrapping controls for small genes...")
      tmpCtrl <- tmp[isCtrl,]
      tmpCtrlChr <- lapply(split(tmpCtrl[,c(1:3)],tmpCtrl$Gene), as.matrix) #the command do.call(rbind,) is much faster on lists of matrices than of data frames
      tmpCtrlMat <- lapply(split(tmpCtrl[,-c(1:3),drop=F],tmpCtrl$Gene), as.matrix)
    }
    scored2 <- list()
    sizes <- unique(Nguide[Nguide<mode])
    for(x in sizes[sizes>1]) {
      print(paste0("Size = ",x))
      if(!is.null(shinyF)) shinyF(0,"Essentiality analysis",paste0("Size = ",x))
      isGene <- tmp$Gene%in%names(Nguide)[which(Nguide==x)]
      tmpLess <- tmp[isGene,]
      if(uniform.weights==FALSE & empirical.weights==FALSE) {
        # Generate multiple models using resampled controls (from all controls the same amount of sgRNAs):
        tmpCtrlChr1 <- data.frame(do.call(rbind,lapply(tmpCtrlChr, function(y) { return(y[sample(1:mode,x),]) })),stringsAsFactors=FALSE,check.names=FALSE)
        tmpCtrlMat1 <- data.frame(do.call(rbind,lapply(tmpCtrlMat, function(y) { return(y[sample(1:mode,x),,drop=F]) })),stringsAsFactors=FALSE,check.names=FALSE)
        tmp2 <- rbind(tmpLess, cbind(tmpCtrlChr1,tmpCtrlMat1))
        scored2[[x]] <- get.weights(tmp2,pctrls,nctrls)$Scored
        for(i in 2:10) {
            tmpCtrlMat1 <- data.frame(do.call(rbind,lapply(tmpCtrlMat, function(y) { return(y[sample(1:mode,x),,drop=F]) })),stringsAsFactors=FALSE,check.names=FALSE)
            tmp2 <- rbind(tmpLess, cbind(tmpCtrlChr1,tmpCtrlMat1))
            scored2[[x]][,-c(1,2)] <- scored2[[x]][,-c(1,2)] + get.weights(tmp2,pctrls,nctrls)$Scored[,-c(1,2)]
        }
        scored2[[x]][,-c(1,2)] <- scored2[[x]][,-c(1,2)]/10
      } else scored2[[x]] <- use.weights(tmpLess,uniform.weights,empirical.weights)
    }
    scored2 <- do.call(rbind,scored2); scored2 <- scored2[scored2$Gene%in%less,]
    scored1 <- rbind(scored1,scored2); scored1 <- scored1[order(scored1$Gene),]
  }
  if(uniform.weights==FALSE & empirical.weights==FALSE) {
    return(list(Scores=scored1,A=Ai,B=Bi,W=wk))
  } else return(list(Scores=scored1))
}
