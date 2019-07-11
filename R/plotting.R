# MoPAC internal plotting functions
# version: 1.0
# date: May 2018
# description: plotting functions used to generate reports
# author:  Oscar Villarreal
# affiliation: University of Texas MD Anderson Cancer Center. Laboratory of Dr. Han Xu
# contact: oscardvillarreal AT gmail.com

### FASTQ ----------------------------------------------------------------------------------------------------------------

plot.mapped1 <- function(anno){
  p1 <- ggplot2::ggplot(anno, ggplot2::aes_string("File","`Mapped%`")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5)) +
    ggplot2::xlab("") +
    ggplot2::ylab("% mapped")
  return(p1)
}

plot.mapped2 <- function(anno){
  anno$Unmapped <- anno$Total_reads-anno$Mapped
  my_data <- reshape2::melt(anno[,c("File","Mapped","Unmapped")], id.vars="File", variable.name="Mapping", value.name="Depth")
  p2 <- ggplot2::ggplot(my_data, ggplot2::aes_string(x="File", y="Depth", fill="Mapping")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(
      legend.title=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),legend.position="top") +
    ggplot2::xlab("") +
    ggplot2::ylab("No. of reads in the fastq file")
  return(p2)
}

plot.zeros0 <- function(anno){
  p1 <- ggplot2::ggplot(anno,ggplot2::aes_string("File","Zeros")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(
      legend.title=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5)) +
    ggplot2::xlab("File") +
    ggplot2::ylab("Zeros")
  return(p1)
}

plot.gini <- function(anno){
  p1 <- ggplot2::ggplot(anno,ggplot2::aes_string("File","Gini")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5)) +
    ggplot2::xlab("File") +
    ggplot2::ylab("Gini")
  return(p1)
}

### Read counts ----------------------------------------------------------------------------------------------------------

plot.classification <- function(id,counts) {
  id1 <- id[id$sgRNA%in%rownames(counts),c("sgRNA","Gene","Category")]
  id2 <- id1[!id1$Category=="Nontargeting",]
  tmp1 <- data.frame(table(id1$Category))
  tmp2 <- data.frame(table(unique(id2[,c("Gene","Category")])$Category))
  names(tmp1)[1] <- names(tmp2)[1] <- "Category"
  p1 <- ggplot2::ggplot(tmp1, ggplot2::aes_string(x="factor(1)", y="Freq", fill="Category")) +
    ggplot2::ggtitle("sgRNA level") +
    ggplot2::geom_bar(stat="identity") +
    ggplot2:: coord_polar(theta = "y") +
    ggplot2::theme(
      axis.text.y=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      legend.title=ggplot2::element_blank(),
      text=ggplot2::element_text(face="bold",color="black",size=12)) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_fill_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",Nontargeting="blue"))
  p2 <- ggplot2::ggplot(tmp2, ggplot2::aes_string(x="factor(1)", y="Freq", fill="Category")) +
    ggplot2::ggtitle("Gene level") +
    ggplot2::geom_bar(stat="identity") +
    ggplot2:: coord_polar(theta="y") +
    ggplot2::theme(
      axis.text.y=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      legend.title=ggplot2::element_blank(),
      text=ggplot2::element_text(face="bold",color="black",size=12)) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_fill_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4"))
  tg1 = gridExtra::tableGrob(tmp1,rows=NULL)
  tg2 = gridExtra::tableGrob(tmp2,rows=NULL)
  return(gridExtra::grid.arrange(p1,p2,tg1,tg2,nrow=2))
}

plot.sizes <- function(library, counts_mat){
  lib1 <- library[library$sgRNA%in%rownames(counts_mat),c("sgRNA","Gene","Category")]
  lib2 <- lib1[!lib1$Category=="Nontargeting",]
  tmp1 <- split(lib2,lib2$Category)
  tmp2 <- do.call(rbind,lapply(1:length(tmp1),function(x){
    out <- data.frame(table(table(tmp1[[x]]$Gene)))
    names(out) <- c("Content","Frequency")
    out$Category <- names(tmp1)[x]
    return(out)
  }))
  tmp2$Category <- factor(tmp2$Category, levels=c("Essential","SEED","Nonessential","Nontargeting","Other"))
  tmp2$Content <- as.character(tmp2$Content)
  p1 <- ggplot2::ggplot(tmp2, ggplot2::aes_string(x="Content", y="Frequency", fill="Category")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
      legend.position="top",
      legend.title=ggplot2::element_blank(),
      text=ggplot2::element_text(face="bold",color="black",size=15)) +
    ggplot2::xlab("sgRNAs per gene") +
    ggplot2::ylab("No.of genes") +
    ggplot2::scale_fill_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
  return(p1)
}

plot.depths <- function(counts, anno){
  # tmp2 <- subset(counts, !duplicated(counts$sgRNA)) #if some sgRNAs are repeated, keep only the first occurence
  if(!"Category"%in%colnames(counts)) counts$Category <- "Other"
  tmp3 <- do.call(rbind,lapply(split(counts[,anno$Sample],counts$Category),colSums))
  tmp3 <- reshape2::melt(tmp3,value.name="Depth",varnames=c("Category","Sample"))
  tmp3$Category <- factor(tmp3$Category, levels=c("Essential","SEED","Nonessential","Nontargeting","Other"))
  p1 <- ggplot2::ggplot(tmp3, ggplot2::aes_string(x="Sample", y="Depth", fill="Category")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
      legend.title=ggplot2::element_blank(),
      axis.title.x=ggplot2::element_blank(),
      legend.position="top",text=ggplot2::element_text(face="bold",color="black",size=15)) +
    ggplot2::ylab("No. of counts") +
    ggplot2::scale_fill_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
  return(p1)
}

plot.zeros <- function(counts, anno){
  if(!"Category"%in%colnames(counts)) counts$Category <- "Other"
  # tmp2 <- subset(counts, !duplicated(counts$sgRNA)) #if some sgRNAs are repeated, keep only the first occurence
  tmp3 <- do.call(rbind,lapply(split(counts[,anno$Sample],counts$Category), function(x) apply(x,2,function(y) sum(y==0))))
  tmp3 <- reshape2::melt(tmp3,value.name="Zeroes",varnames=c("Category","Sample"))
  tmp3$Category <- factor(tmp3$Category, levels=c("Essential","SEED","Nonessential","Nontargeting","Other"))
  p1 <- ggplot2::ggplot(tmp3, ggplot2::aes_string(x="Sample", y="Zeroes", fill="Category")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(
      axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
      legend.title=ggplot2::element_blank(),
      axis.title.x=ggplot2::element_blank(),
      legend.position="top",text=ggplot2::element_text(face="bold",color="black",size=15)) +
    ggplot2::ylab("No. of zeros") +
    ggplot2::scale_fill_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
  return(p1)
}

# plot.cumulative <- function(counts,anno){
#   tmp0 <- data.frame(sgRNA=rownames(counts),counts,check.names=F)
#   tmp1 <- reshape2::melt(tmp0, id.vars="sgRNA",value.name="Reads",variable.name="Sample")
#   tmp1 <- merge(tmp1,anno)
#   p1 <- ggplot2::ggplot(tmp1, ggplot2::aes(log10(Reads+1),col=Sample)) +
#     ggplot2::stat_ecdf(size=1) +
#     ggplot2::ylab("CDF")
#   return(p1)
# }

### Normalized counts ---------------------------------------------------------------------------------------------------

plot.pseudocount <- function(pseudo) {
  variance <- pseudo$Variance
  optimal <- pseudo$Optimal
  melted <- reshape2::melt(variance,id.vars = "Pseudocount",value.name = "Variance",variable.name = "Category")
  melted <- rbind(melted,data.frame(Pseudocount=optimal$Pseudocount,Variance=optimal$Variance,Category="Combined"))
  ggplot2::ggplot(melted, ggplot2::aes_string("Pseudocount","Variance")) +
    ggplot2::facet_wrap(~Category,nrow=1,scales="free",strip.position = "top") +
    ggplot2::geom_line(size=2) +
    ggplot2::geom_point(size=5) +
    ggplot2::geom_point(data=melted[melted$Pseudocount==optimal$Pseudocount,],size=5,ggplot2::aes(col="Selected")) +
    ggplot2::theme(
      legend.position=c(0.15,0.9),
      text=ggplot2::element_text(face="bold",color="black",size=15),
      strip.text = ggplot2::element_text(face="bold",color="white",size=15),
      strip.background=ggplot2::element_rect(fill="black"),
      legend.title=ggplot2::element_blank(),
      panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
      panel.background=ggplot2::element_blank())
}

plot.weights <- function(weights) {
  weights$Rank <- as.factor(weights$Rank)
  ggplot2::ggplot(weights, ggplot2::aes_string("Rank","Weight")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(
      text=ggplot2::element_text(face="bold",color="black",size=15),
      title=ggplot2::element_text(face="bold",color="black",size=15),
      strip.text = ggplot2::element_text(face="bold",color="white",size=15),
      legend.text=ggplot2::element_text(size=13),
      panel.background=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_blank(),
      panel.grid.minor=ggplot2::element_blank(),
      panel.border=ggplot2::element_rect(fill=NA,size=1),
      legend.position="bottom",
      legend.justification=c(0,0),
      legend.title=ggplot2::element_blank(),
      strip.background=ggplot2::element_rect(fill="black"),
      axis.ticks=ggplot2::element_line(size=1,color="black"),
      axis.text = ggplot2::element_text(size=15,face="bold",color="black"),
      axis.text.x=ggplot2::element_text(size=15),
      axis.title=ggplot2::element_text(size=15,face="bold",color="black"),
      axis.ticks.length=grid::unit(0.15,"cm"))
}

plot.scaling <- function(scaling) {
  ggplot2::ggplot(scaling, ggplot2::aes_string("Sample","Scaling")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(
      text=ggplot2::element_text(face="bold",color="black",size=15),
      title=ggplot2::element_text(face="bold",color="black",size=15),
      strip.text = ggplot2::element_text(face="bold",color="white",size=15),
      legend.text=ggplot2::element_text(size=13),
      panel.background=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_blank(),
      panel.grid.minor=ggplot2::element_blank(),
      panel.border=ggplot2::element_rect(fill=NA,size=1),
      legend.position="bottom",
      legend.justification=c(0,0),
      legend.title=ggplot2::element_blank(),
      strip.background=ggplot2::element_rect(fill="black"),
      axis.ticks=ggplot2::element_line(size=1,color="black"),
      axis.text = ggplot2::element_text(size=15,face="bold",color="black"),
      axis.text.x=ggplot2::element_text(size=15,angle=90,hjust=1,vjust=0.5),
      axis.title=ggplot2::element_text(size=15,face="bold",color="black"),
      axis.ticks.length=grid::unit(0.15,"cm"))
}

plot.shifting <- function(shifting) {
  ggplot2::ggplot(shifting,ggplot2::aes_string("Sample","Shifting")) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(
      text=ggplot2::element_text(face="bold",color="black",size=15),
      title=ggplot2::element_text(face="bold",color="black",size=15),
      strip.text = ggplot2::element_text(face="bold",color="white",size=15),
      legend.text=ggplot2::element_text(size=13),
      panel.background=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_blank(),
      panel.grid.minor=ggplot2::element_blank(),
      panel.border=ggplot2::element_rect(fill=NA,size=1),
      legend.position="bottom",
      legend.justification=c(0,0),
      legend.title=ggplot2::element_blank(),
      strip.background=ggplot2::element_rect(fill="black"),
      axis.ticks=ggplot2::element_line(size=1,color="black"),
      axis.text = ggplot2::element_text(size=15,face="bold",color="black"),
      axis.text.x=ggplot2::element_text(size=15,angle=90,hjust=1,vjust=0.5),
      axis.title=ggplot2::element_text(size=15,face="bold",color="black"),
      axis.ticks.length=grid::unit(0.15,"cm"))
}

### Fold changes ---------------------------------------------------------------------------------------------------

plot.correlation <- function(counts){
  tmp1 <- stats::cor(counts,method="pearson")
  return(gplots::heatmap.2(tmp1, cexRow=1, cexCol=1, margins=c(15,15), #cexRow and cexCol change the font size
            key = FALSE, dendrogram = "row", lhei = c(0.1,1), lwid=c(0.1,1),
            trace="none", density.info="none", #show trace lines?
            cellnote=round(tmp1,2), notecol="white", notecex=1.2,
            col=grDevices::colorRampPalette(c("black","blue"))))
}

plot.pca <- function(pca) {
  p1 <- ggplot2::ggplot(pca, ggplot2::aes_string("Z1","Z2",label="Sample")) +
    ggrepel::geom_text_repel() +
    ggplot2::geom_point() +
    ggplot2::ggtitle("2nd vs 1st") +
    ggplot2::xlab("First Principal Component for Samples") +
    ggplot2::ylab("Second Principal Component")
  p2 <- ggplot2::ggplot(pca, ggplot2::aes_string("Z1","Z3",label="Sample")) +
    ggrepel::geom_text_repel() +
    ggplot2::geom_point() +
    ggplot2::ggtitle("3rd vs 1st") +
    ggplot2::xlab("First Principal Component for Samples") +
    ggplot2::ylab("Third Principal Component")
  p3 <- gridExtra::grid.arrange(p1,p2,nrow=1)
  return(p3)
}

plot.pve <- function(pve_sample) {
  pve <- data.frame(PC=seq(1:length(pve_sample)),pve_sample=pve_sample)
  cpve <- data.frame(PC=seq(1:length(pve_sample)),Cpve_sample=cumsum(pve_sample))
  p1 <- ggplot2::ggplot(pve, ggplot2::aes_string("PC","pve_sample")) +
    ggplot2::geom_point(col="darkblue",shape="O",size=5) +
    ggplot2::geom_line(col="blue") +
    ggplot2::xlab("Principal Component") +
    ggplot2::ylab("Proportion of Variance Explained")
  p2 <- ggplot2::ggplot(cpve, ggplot2::aes_string("PC","Cpve_sample")) +
    ggplot2::geom_line(col="blue") +
    ggplot2::geom_point(col="darkblue",shape="O",size=5) +
    ggplot2::xlab("Principal Component") +
    ggplot2::ylab("Cumulative PVE")
  p3 <- gridExtra::grid.arrange(p1,p2,nrow=1)
  return(p3)
}

plot.distributions <- function(logged,anno,formula){
  my_data <- reshape2::melt(logged,id.vars=c("sgRNA","Gene","Category"),value.name="Log2",variable.name="Sample")
  my_data <- merge(my_data,anno)
  p1 <- ggplot2::ggplot(my_data[my_data$Category!="Other",],ggplot2::aes_string("Log2",color="Category")) +
    ggplot2::theme(
      legend.position="top",
      legend.title=ggplot2::element_blank(),
      panel.background=ggplot2::element_blank(),
      panel.grid.major=ggplot2::element_blank(),
      panel.grid.minor=ggplot2::element_blank(),
      axis.line=ggplot2::element_line(colour="black",size=1),
      panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
      text=ggplot2::element_text(face="bold",color="black",size=15),
      axis.ticks=ggplot2::element_line(colour="black",size=1),
      strip.background=ggplot2::element_rect(fill="#0033cc"),
      strip.text = ggplot2::element_text(face="bold",color="white")
    ) +
    ggplot2::geom_density(ggplot2::aes_string(fill="Category"),stat="density",alpha=0.2,lwd=0) +
    ggplot2::geom_line(ggplot2::aes_string(fill="Category"),stat="density",size=3) +
    ggplot2::geom_line(data=my_data[my_data$Category=="Other",],stat="density",linetype="dashed",size=2) +
    ggplot2::ylab("Density") +
    ggplot2::facet_grid(stats::as.formula(formula)) +
    ggplot2::scale_fill_discrete(guide=FALSE) +
    ggplot2::scale_color_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
  print(p1)
}

plot.reproducibility <- function(dataset,title) {
  if(ncol(dataset)<3) return(NULL)
  colors <- c(Other=1,Essential=2,Nonessential=3,Nontargeting=4)
  others <- dataset$Category=="Other"
  graphics::pairs(dataset[,-1],main=title,
    lower.panel = function(x,y) {
      graphics::par(usr=c(0,1,0,1))
      graphics::text(0.5,0.5,round(stats::cor(x,y),digits=3),cex=2)
    },
    upper.panel = function(x,y) {
      graphics::points(x[others],y[others],pch=19,col=colors[dataset$Category[others]])
      graphics::points(x[!others],y[!others],pch=19,col=colors[dataset$Category[!others]])
      graphics::abline(a=0,b=1,col="green")
    }
  )
}

### Replicate-averaged fold changes -------------------------------------------------------------------------------

plot.folds <- function(folded, formula, level){
  if(level=="sgRNA"){
    tmp1 <- reshape2::melt(folded,id.vars=c("sgRNA","Gene","Category"),value.name="Fold",variable.name="Sample")
  } else {
    tmp1 <- folded[!folded$Category=="Nontargeting",]
    tmp1 <- reshape2::melt(tmp1,id.vars=c("Gene","Category"),value.name="Fold",variable.name="Sample")
  }
  # tmp1 <- merge(tmp1,anno1)
  p1 <- ggplot2::ggplot(tmp1,ggplot2::aes_string("Category","Fold",col="Category")) +
    ggplot2::geom_boxplot(lwd=2,fatten=1/4) +
    ggplot2::stat_summary(fun.y=mean,geom="point",shape=18,size=2) +
    ggplot2::facet_wrap(stats::as.formula(formula)) +
    # ggplot2::facet_grid(stats::as.formula(formula)) +
    ggplot2::theme(
      legend.position="bottom",
      panel.background=ggplot2::element_blank(),
      panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
      legend.title=ggplot2::element_blank(),
      axis.ticks.x=ggplot2::element_blank(),
      strip.background=ggplot2::element_rect(fill="#0033cc"),
      strip.text = ggplot2::element_text(face="bold",color="white"),
      axis.text.x=ggplot2::element_blank()) +
    ggplot2::scale_color_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
  return(p1)
}





