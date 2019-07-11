options(shiny.maxRequestSize=1000*1024^2) #allow large file upload

shinyServer(function(input, output, session){

  ready <- reactiveValues(fastq=0, library=0) # to control the appearance/dissapearance of widgets
  output$fastq <- reactive({ return(ready$fastq) }) # an output to detect the value of the flag in the UI
  outputOptions(output, name="fastq", suspendWhenHidden=FALSE) # a condition to use in the conditional panels of the UI

  ### PLOTTING FUNCTIONS ### ---------------------------------------------------

  plot.mapped1 <- function(anno){
    p1 <- ggplot2::ggplot(anno, ggplot2::aes_string("File","`Mapped%`")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
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
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
      ggplot2::xlab("") +
      ggplot2::ylab("No. of reads in the fastq file")
    return(p2)
  }

  plot.gini <- function(anno){
    colnames(anno)[1] <- "Sample"
    p1 <- ggplot2::ggplot(anno,ggplot2::aes_string("Sample","Gini")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
                     text=ggplot2::element_text(face="bold",color="black",size=15),
                     axis.title.x=ggplot2::element_blank()) +
      ggplot2::ylab("Gini index") # ggplot2::xlab("File") +
    return(p1)
  }

  plot.classification <- function(id,counts) {
    if(!"Category"%in%colnames(id)) id$Category <- "Other"
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
    return(gridExtra::grid.arrange(p1,tg1,p2,tg2,nrow=2))
  }

  plot.sizes <- function(library, counts_mat){
    if(!"Category"%in%colnames(library)) library$Category <- "Other"
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

  plot.depths <- function(counts){
    # tmp1 <- data.frame(sgRNA=rownames(counts),counts,check.names=F)
    # tmp1 <- merge(tmp1[tmp1$sgRNA%in%id$sgRNA,],id,all.x=T)
    if(!"Category"%in%colnames(counts)) counts$Category <- "Other"
    # tmp2 <- subset(counts, !duplicated(counts$sgRNA)) #if some sgRNAs are repeated, keep only the first occurence
    tmp3 <- do.call(rbind,lapply(split(counts[,-which(colnames(counts)%in%c("sgRNA","Gene","Category")),drop=F],counts$Category),colSums))
    tmp3 <- reshape2::melt(tmp3,value.name="Depth",varnames=c("Category","Sample"))
    tmp3$Category <- factor(tmp3$Category, levels=c("Essential","SEED","Nonessential","Nontargeting","Other"))
    p1 <- ggplot2::ggplot(tmp3, ggplot2::aes_string(x="Sample", y="Depth", fill="Category")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        legend.title=ggplot2::element_blank(),
        axis.title.x=ggplot2::element_blank(),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
      ggplot2::ylab("No. of counts") +
      ggplot2::scale_fill_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
    return(p1)
  }

  plot.zeros <- function(counts){
    # tmp1 <- data.frame(sgRNA=rownames(counts),counts,check.names=F)
    # tmp1 <- merge(tmp1[tmp1$sgRNA%in%id$sgRNA,],id,all.x=T)
    if(!"Category"%in%colnames(counts)) counts$Category <- "Other"
    # tmp2 <- subset(counts, !duplicated(counts$sgRNA)) #if some sgRNAs are repeated, keep only the first occurence
    tmp3 <- do.call(rbind,lapply(split(counts[,-which(colnames(counts)%in%c("sgRNA","Gene","Category")),drop=F],counts$Category), function(x) apply(x,2,function(y) sum(y==0))))
    tmp3 <- reshape2::melt(tmp3,value.name="Zeroes",varnames=c("Category","Sample"))
    tmp3$Category <- factor(tmp3$Category, levels=c("Essential","SEED","Nonessential","Nontargeting","Other"))
    p1 <- ggplot2::ggplot(tmp3, ggplot2::aes_string(x="Sample", y="Zeroes", fill="Category")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(
        axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
        legend.title=ggplot2::element_blank(),
        axis.title.x=ggplot2::element_blank(),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15)) +
      ggplot2::ylab("No. of zeros") +
      ggplot2::scale_fill_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
    return(p1)
  }

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
        # legend.position="top",
        legend.position=c(0.15,0.9),
        text=ggplot2::element_text(face="bold",color="black",size=15),
        strip.text = ggplot2::element_text(face="bold",color="white",size=15),
        strip.background=ggplot2::element_rect(fill="black"),
        legend.title=ggplot2::element_blank(),
        panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
        panel.background=ggplot2::element_blank())
  }

  plot.correlation <- function(counts,method){
    tmp1 <- stats::cor(counts,method=method)
    return(gplots::heatmap.2(tmp1, cexRow=1, cexCol=1, margins=c(15,15), #cexRow and cexCol change the font size
                             key = FALSE, dendrogram = "row", lhei = c(0.1,1), lwid=c(0.1,1),
                             trace="none", density.info="none", #show trace lines?
                             cellnote=round(tmp1,2), notecol="white", notecex=1.2,
                             col=grDevices::colorRampPalette(c("black","blue"))))
  }

  plot.pca <- function(pca,Y) {
    p1 <- ggplot2::ggplot(pca, ggplot2::aes_string("Z1",Y,label="Sample")) +
      ggrepel::geom_text_repel() +
      ggplot2::geom_point() +
      ggplot2::ggtitle("2nd vs 1st") +
      ggplot2::xlab("First Principal Component for Samples") +
      ggplot2::ylab("Second Principal Component") +
      ggplot2::theme_linedraw()
    p2 <- ggplot2::ggplot(pca, ggplot2::aes_string("Z1",Y,label="Sample")) +
      ggrepel::geom_text_repel() +
      ggplot2::geom_point() +
      ggplot2::ggtitle("3rd vs 1st") +
      ggplot2::xlab("First Principal Component for Samples") +
      ggplot2::ylab("Third Principal Component") +
      ggplot2::theme_linedraw()
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
      ggplot2::ylab("Proportion of Variance Explained") +
      ggplot2::theme_linedraw()
    p2 <- ggplot2::ggplot(cpve, ggplot2::aes_string("PC","Cpve_sample")) +
      ggplot2::geom_line(col="blue") +
      ggplot2::geom_point(col="darkblue",shape="O",size=5) +
      ggplot2::xlab("Principal Component") +
      ggplot2::ylab("Cumulative PVE") +
      ggplot2::theme_linedraw()
    p3 <- gridExtra::grid.arrange(p1,p2,nrow=1)
    return(p3)
  }

  # plot.distributions <- function(logged,anno,formula){
  #   if(!"Category"%in%colnames(logged)) logged$Category <- "Other"
  #   my_data <- reshape2::melt(logged,id.vars=c("sgRNA","Gene","Category"),value.name="Log2",variable.name="Sample")
  #   my_data <- merge(my_data,anno)
  #   if(input$lc.style=="Density plot") {
  #     p1 <- ggplot2::ggplot(my_data[my_data$Category!="Other",],ggplot2::aes_string("Log2",color="Category")) +
  #       ggplot2::theme(
  #         title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
  #         legend.text = ggplot2::element_text(size=15,color="black"),
  #         legend.position="top",
  #         legend.title=ggplot2::element_blank(),
  #         panel.background=ggplot2::element_blank(),
  #         panel.grid.major=ggplot2::element_blank(),
  #         panel.grid.minor=ggplot2::element_blank(),
  #         axis.line=ggplot2::element_line(colour="black",size=1),
  #         panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
  #         text=ggplot2::element_text(face="bold",color="black",size=15),
  #         axis.ticks=ggplot2::element_line(colour="black",size=1),
  #         strip.background=ggplot2::element_rect(fill="#0033cc"),
  #         strip.text = ggplot2::element_text(face="bold",color="white")
  #       ) +
  #       ggplot2::geom_density(ggplot2::aes_string(fill="Category"),stat="density",alpha=0.2,lwd=0) +
  #       ggplot2::geom_line(ggplot2::aes_string(fill="Category"),stat="density",size=3) +
  #       ggplot2::geom_line(data=my_data[my_data$Category=="Other",],stat="density",linetype="dashed",size=2) +
  #       ggplot2::ylab("Density") +
  #       ggplot2::facet_grid(stats::as.formula(formula)) +
  #       ggplot2::scale_fill_discrete(guide=FALSE) +
  #       ggplot2::scale_color_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
  #   } else {
  #     p1 <- ggplot2::ggplot(my_data,ggplot2::aes_string("Category","Log2",col="Category")) +
  #       ggplot2::xlab("") +
  #       ggplot2::geom_boxplot(lwd=2,fatten=1/4) +
  #       ggplot2::stat_summary(fun.y=mean,geom="point",shape=18,size=2) +
  #       ggplot2::facet_grid(stats::as.formula(formula)) +
  #       ggplot2::theme(
  #         title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
  #         legend.text = ggplot2::element_text(size=15,color="black"),
  #         legend.position="top",
  #         text=ggplot2::element_text(face="bold",color="black",size=15),
  #         panel.background=ggplot2::element_blank(),
  #         panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
  #         legend.title=ggplot2::element_blank(),
  #         axis.ticks.x=ggplot2::element_blank(),
  #         strip.background=ggplot2::element_rect(fill="#0033cc"),
  #         strip.text = ggplot2::element_text(face="bold",color="white"),
  #         axis.text.x=ggplot2::element_blank()) +
  #       ggplot2::scale_color_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
  #   }
  #   print(p1)
  # }

  plot.folds <- function(folded, loglevel) {
    if(input$lfc.style=="Box plot") {
      p1 <- ggplot2::ggplot(folded,ggplot2::aes_string("Category","Value",col="Category")) +
        ggplot2::xlab("") +
        ggplot2::geom_boxplot(lwd=2,fatten=1/4) +
        ggplot2::stat_summary(fun.y=mean,geom="point",shape=18,size=2) +
        ggplot2::facet_wrap(stats::as.formula("~Condition+Replicate"),nrow=1) +
        ggplot2::theme(
          title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
          legend.text = ggplot2::element_text(size=15,color="black"),
          legend.position="top",
          text=ggplot2::element_text(face="bold",color="black",size=15),
          panel.background=ggplot2::element_blank(),
          panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
          legend.title=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          strip.background=ggplot2::element_rect(fill="#0033cc"),
          strip.text = ggplot2::element_text(face="bold",color="white"),
          axis.text.x=ggplot2::element_blank())
      ggplot2::scale_color_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
      if(loglevel=="Fold"){ p1 <- p1 + ggplot2::ylab("Log-fold-change") } else p1 <- p1 + ggplot2::ylab("Log-read-count")
    } else {
      p1 <- ggplot2::ggplot(folded[folded$Category!="Other",],ggplot2::aes_string("Value",color="Category")) +
        ggplot2::theme(
          title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
          legend.text = ggplot2::element_text(size=15,color="black"),
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
        ggplot2::geom_line(data=folded[folded$Category=="Other",],stat="density",linetype="dashed",size=2) +
        ggplot2::facet_wrap(stats::as.formula("~Condition+Replicate"),nrow=1,strip.position="top") + ggplot2::coord_flip() +
        ggplot2::scale_fill_discrete(guide=FALSE)
      ggplot2::scale_color_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
      if(loglevel=="Fold"){ p1 <- p1 + ggplot2::xlab("Log-fold-change") } else p1 <- p1 + ggplot2::xlab("Log-read-count")
    }
    # if("sgRNA"%in%colnames(folded)) { p1 <- p1 + ggplot2::ggtitle("sgRNA log-fold-changes") } else p1 <- p1 + ggplot2::ggtitle("Gene log-fold-changes")
    return(p1)
  }

  plot.ssmd <- function(folds) {
    SSMD1 <- apply(folds[,-which(colnames(folds)%in%c("sgRNA","Gene","Category"))],2,function(y){
      muP <- mean(y[folds$Category=="Essential"])
      muN <- mean(y[folds$Category=="Nonessential"])
      varP <- var(y[folds$Category=="Essential"])
      varN <- var(y[folds$Category=="Nonessential"])
      return((muN-muP)/sqrt(varP+varN))
    })
    SSMD2 <- data.frame(Sample=names(SSMD1),SSMD=SSMD1)
    ggplot2::ggplot(SSMD2,ggplot2::aes_string("Sample","SSMD")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=0.5),
                     text=ggplot2::element_text(face="bold",color="black",size=15),
                     axis.title.x=ggplot2::element_blank()) +
      ggplot2::ylab("SSMD") # ggplot2::xlab("File") +
  }

  plot.scores <- function(scored,anno,formula="~Condition+Replicate"){
    # print(file=stderr(),head(scored$scored_reps))
    # print(file=stderr(),anno)

    if(!"Category"%in%colnames(scored)) scored$Category <- "Other"
    my_data <- reshape2::melt(scored[scored$Category!="Nontargeting",],id.vars=c("Gene","Category"),value.name="Score",variable.name="Sample")
    my_data <- merge(my_data,anno)
    if(input$ea.style=="Density plot") {
      p1 <- ggplot2::ggplot(my_data[my_data$Category!="Other",],ggplot2::aes_string("Score",color="Category")) +
        ggplot2::theme(
          title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
          legend.text = ggplot2::element_text(size=15,color="black"),
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
        ggplot2::facet_wrap(stats::as.formula("~Condition+Replicate"),nrow=1,strip.position="top") + ggplot2::coord_flip() +
        # ggplot2::facet_grid(stats::as.formula(formula)) +
        ggplot2::scale_fill_discrete(guide=FALSE) +
        ggplot2::scale_color_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
    } else {
      p1 <- ggplot2::ggplot(my_data,ggplot2::aes_string("Category","Score",col="Category")) +
        ggplot2::xlab("") +
        ggplot2::geom_boxplot(lwd=2,fatten=1/4) +
        ggplot2::stat_summary(fun.y=mean,geom="point",shape=18,size=2) +
        ggplot2::facet_wrap(stats::as.formula("~Condition+Replicate"),nrow=1) +
        # ggplot2::facet_grid(stats::as.formula(formula)) +
        ggplot2::theme(
          title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
          legend.text = ggplot2::element_text(size=15,color="black"),
          legend.position="top",
          text=ggplot2::element_text(face="bold",color="black",size=15),
          panel.background=ggplot2::element_blank(),
          panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
          legend.title=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(),
          strip.background=ggplot2::element_rect(fill="#0033cc"),
          strip.text = ggplot2::element_text(face="bold",color="white"),
          axis.text.x=ggplot2::element_blank()) +
        ggplot2::scale_color_manual(name="Category",values=c(Other="gray50",Essential="#F8766D",Nonessential="#00BFC4",SEED="red",Nontargeting="blue"))
    }
    print(p1)
  }

  plot.reproducibility <- function(dataset,title) {
    if(!"Category"%in%colnames(dataset)) dataset$Category <- "Other"
    if(ncol(dataset)<3) return(NULL)
    colors <- c(Other=1,Essential=2,Nonessential=3,Nontargeting=4)
    others <- dataset$Category=="Other"
    graphics::pairs(dataset[,-which(colnames(dataset)=="Category"),drop=FALSE],#main=title,
                    lower.panel = function(x,y) {
                      graphics::par(usr=c(0,1,0,1))
                      graphics::text(0.5,0.5,round(stats::cor(x,y),digits=3),cex=2)
                    },
                    upper.panel = function(x,y) {
                      graphics::abline(a=0,b=1,col="gray",lwd=5)
                      graphics::points(x[others],y[others],pch=19,col=colors[dataset$Category[others]])
                      graphics::points(x[!others],y[!others],pch=19,col=colors[dataset$Category[!others]])
                    },
                    oma=c(3,3,9,3)
    )
    par(xpd = TRUE)
    graphics::legend("topright",col=1:length(colors),legend=names(colors),pch=16,bty='n',ncol=length(colors),pt.cex=3,cex=1.5)
  }

  ### INPUT ### -----------------------------------------------------------

  volumes <- shinyFiles::getVolumes()
  paths <- reactiveValues(out=character(0),fastq=character(0))

  # Output file:
  # A) in the "FASTQ Alignment" tab:
  observeEvent(input$dir_out,{
    shinyFiles::shinyDirChoose(input,'dir_out',roots=volumes,session=session,restrictions=system.file(package='base'))
    paths$out <- shinyFiles::parseDirPath(volumes,input$dir_out)
  })
  output$dir_out_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<b>Please select a working directory:</b><br/>")))
    } else return(HTML(paste0("<b>Working directory selected:</b><br/>")))
  })
  output$dir_out_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    return(paths$out)
  })
  observeEvent(input$dir_fastq,{
    shinyFiles::shinyDirChoose(input,'dir_fastq',roots=volumes,session=session,restrictions=system.file(package='base'))
    paths$fastq <- shinyFiles::parseDirPath(volumes,input$dir_fastq)
  })
  output$dir_fastq_text1 <- renderUI({
    if(identical(paths$fastq,character(0))) {
      return(HTML(paste0("<b>Please select the directory containing FASTQ files:</b><br/>")))
    } else return(HTML(paste0("<b>FASTQ directory selected:</b><br/>")))
  })
  output$dir_fastq_text2 <- renderText({
    req(paths$fastq!=0)
    return(paths$fastq)
  })
  output$q1 <- downloadHandler(
    filename = "sample_counts.txt",
    content = function(file) {
      write.table(MoPAC::dang_cck81,file,quote=F,row.names=F,col.names=T,sep='\t')
    }
  )
  output$q3 <- downloadHandler(
    filename = "sample_library.txt",
    content = function(file) {
      write.table(MoPAC::sgRNA_library[,c("sgRNA","Gene")],file,quote=F,row.names=F,col.names=T,sep='\t')
    }
  )
  # B) in the "Quality control" tab:
  observeEvent(input$dir_qc,{
    shinyFiles::shinyDirChoose(input,'dir_qc',roots=volumes,session=session,restrictions=system.file(package='base'))
    paths$out <- shinyFiles::parseDirPath(volumes,input$dir_qc)
  })
  output$dir_qc_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<b>Please select a working directory:</b><br/>")))
    } else return(HTML(paste0("<b>Working directory selected:</b><br/>")))
  })
  output$dir_qc_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    return(paths$out)
  })
  output$counts_qc <- renderUI({
    if(ready$fastq==0) {
      return(NULL)
    } else return(HTML(paste0("<b>Table of read counts:</b><br/>",ready$fastq)))
  })
  output$library_qc <- renderUI({
    if(ready$library==0) {
      return(NULL)
    } else return(HTML(paste0("<b>sgRNA library file:</b><br/>",ready$library)))
  })
  # C) in the "Essentiality analysis" tab:
  output$dir_ea_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<b>Please run the Quality Control module first.</b><br/>")))
    } else return(HTML(paste0("<b>Working directory selected:</b><br/>")))
  })
  output$dir_ea_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    return(paths$out)
  })
  # D) in the "Differential essentiality analysis" tab:
  output$dir_dea_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<b>Please run the Quality Control module first.</b><br/>")))
    } else return(HTML(paste0("<b>working directory selected:</b><br/>")))
  })
  output$dir_dea_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    return(paths$out)
  })
  # E) in the "Network analysis" tab:
  output$dir_net_text1 <- renderUI({
    if(identical(paths$out,character(0))) {
      return(HTML(paste0("<b>Please run the Quality Control module first.</b><br/>")))
    } else return(HTML(paste0("<b>Working directory selected:</b><br/>")))
  })
  output$dir_net_text2 <- renderText({
    req(!identical(paths$out,character(0)))
    diff_out <- basename(list.files(path=paths$out,pattern="_differential",full.names=TRUE))
    comparisons <- substring(diff_out,1,regexpr("_differential",diff_out)-1)
    updateSelectInput(session, "comparisons",
                      label = "Select condition comparison:",
                      choices = comparisons,
                      selected = comparisons[1]
    )
    return(paths$out)
  })

  ### FASTQ ### -----------------------------------------------------------

  # Read counts from FASTQ files:
  counts <- eventReactive(input$button_run, {
    withProgress(message='Input',value=0,detail="Initializing...",{
      essential=unlist(strsplit(input$essentials,"\n"))
      nonessential=unlist(strsplit(input$nonessentials,"\n"))
      nontargeting=unlist(strsplit(input$nontargeting,"\n"))
      if(identical(essential,character(0))) essential <- NULL
      if(identical(nonessential,character(0))) nonessential <- NULL
      if(identical(nontargeting,character(0))) nontargeting <- NULL

      validate(need(paths$out!=0,"Please enter an output directory."))
      validate(need(paths$fastq!=0,"Please enter a directory containing FASTQ files."))

      fastq <- MoPAC::read.fastq(id=paths$out,
                               fastq=paths$fastq,
                               library=input$file.library$datapath,
                               spacer_start=input$spacer.range[1]-1,
                               spacer_length=input$spacer.range[2]-(input$spacer.range[1]-1),
                               reverse=input$spacer.rev,
                               complement=input$spacer.comp,
                               essential=essential,
                               nonessential=nonessential,
                               nontargeting=nontargeting,
                               shinyF=incProgress)
      incProgress(0,"Input","Finalizing...")

      counts <- fastq$Counts
      annotation <- fastq$Annotation
      library <- fastq$Library
      counts_mat <- as.matrix(counts[,sapply(counts,is.numeric)])
      rownames(counts_mat) <- counts$sgRNA
    })
    showNotification("Please fill out the sample information in the 'File annotation' tab.",action=NULL,duration=20,closeButton=TRUE,type="message")

    return(list(Counts=counts,Library=library,Annotation=annotation,Counts_mat=counts_mat))
  })

  output$mapped2 <- renderPlot({
    req(counts()$Annotation)
    plot.mapped2(counts()$Annotation)
  })
  output$mapped1 <- renderPlot({
    req(counts()$Annotation)
    plot.mapped1(counts()$Annotation)
  })

  output$anno_table <- rhandsontable::renderRHandsontable({
    rtable <- counts()$Annotation
    rownames(rtable) <- NULL
    rhandsontable::rhandsontable(rtable) # converts the R dataframe to rhandsontable object
  })

  output$text_anno1 <- renderUI({
    if(!input$save_anno) {
      return(HTML(paste0("Please fill out the table below. Then click on the red button in order to merge the FASTQ files belonging to the same condition and replicate. <br/>")))
    } else if(!is.null(counts1()$Counts)) return(HTML(paste0("Merged table of counts saved at <br/>",ready$fastq)))
  })

  output$data_input <- DT::renderDataTable({ #DT is necessary to make the row names appear
    if(input$fastq.level=="FASTQ files") {
      req(counts())
      return(counts()$Counts)
    } else if(input$fastq.level=="Conditions") {
      req(counts1())
      return(counts1()$Counts)
    } else return(NULL)
  }, options = list(pageLength = 30), rownames=FALSE)

  ### COUNTS ### -----------------------------------------------------------

  counts1 <- eventReactive(c(input$save_anno,input$button_qc),{
    # Either the FASTQ module must have been run, or the RUN button pressed:
    if(input$save_anno==0 & input$button_qc==0) return(NULL)

    withProgress(message='Read counts',value=0,detail="Initializing...",{
      validate(need(paths$out!=0,"Please enter an output directory first."))

      essential=unlist(strsplit(input$essentials1,"\n"))
      nonessential=unlist(strsplit(input$nonessentials1,"\n"))
      nontargeting=unlist(strsplit(input$nontargeting1,"\n"))
      if(identical(essential,character(0))) essential <- NULL
      if(identical(nonessential,character(0))) nonessential <- NULL
      if(identical(nontargeting,character(0))) nontargeting <- NULL

      if(identical(paths$fastq,character(0))) {
        merged <- MoPAC::read.counts(id=paths$out,
                                     library=input$file.library2$datapath,
                                     counts=input$file.counts$datapath,
                                     essential=essential,
                                     nonessential=nonessential,
                                     nontargeting=nontargeting,
                                     shinyF=incProgress)
      } else {
        validate(need(!is.null(input$anno_table),"Please run the FASTQ alignment module first."))
        rtable <- rhandsontable::hot_to_r(input$anno_table)
        rtable1 <- rtable[!is.na(rtable[,1]),]
        rtable1$Condition <- stringr::str_replace_all(rtable1$Condition,"[^[:alnum:]]","") #remove non-alphanumeric characters in the sample names
        validate(need(sum(rtable1=="")==0,"Please fill out the FASTQ file annotation first."))
        if(sum(rtable1=="")>0)
          showNotification("Please fill out the FASTQ file annotation first.",action=NULL,duration=10,closeButton=TRUE,type="message")
        write.csv(rtable1,file=paste0(paths$out,"/fastq_annotation.csv"),row.names=FALSE,quote=FALSE)
        merged <- MoPAC::read.counts(id=paths$out,
                                     library=counts()$Library,
                                     counts=counts()$Counts,
                                     annotation=rtable1,
                                     essential=essential,
                                     nonessential=nonessential,
                                     nontargeting=nontargeting,
                                     shinyF=incProgress)
        ready$fastq <- paste0(paths$out,"/counts.txt")
        ready$library <- paste0(paths$out,"/library.txt")
      }

      counts <- merged$Counts
      annotation <- merged$Annotation
      library <- merged$Library
      counts_mat <- as.matrix(counts[,sapply(counts,is.numeric)])
      rownames(counts_mat) <- counts$sgRNA
    })
    return(list(Counts=counts,Library=library,Annotation=annotation,Counts_mat=counts_mat))
  })

  output$classification <- renderPlot({
    req(counts1())
    plot.classification(counts1()$Library,counts1()$Counts_mat)
  })
  output$sizes <- renderPlot({
    req(counts1())
    plot.sizes(counts1()$Library,counts1()$Counts_mat)
  })
  output$gini <- renderPlot({
    req(counts1()$Annotation)
    plot.gini(counts1()$Annotation)
  })
  output$depths <- renderPlot({
    req(counts1())
    plot.depths(counts1()$Counts)
  })
  output$zeros <- renderPlot({
    req(counts1())
    plot.zeros(counts1()$Counts)
  })

  ### Quality Control ### -----------------------------------------------

  qc <- eventReactive(input$save_anno1, {
    withProgress(message='Quality Control',value=0,detail="Initializing...",{
      validate(need(paths$out!=0,"Please enter an output directory."))

      if(input$pseudo_default==TRUE) {
        pseudocount <- 0.15
      } else if(input$pseudo_optimize==TRUE) {
        pseudocount <- NULL
      } else pseudocount <- input$pseudo.count

      rtable <- rhandsontable::hot_to_r(input$anno1_table)
      rtable1 <- rtable[!is.na(rtable[,1]),]
      write.csv(rtable1,file=paste0(paths$out,"/annotation.csv"),row.names=FALSE,quote=FALSE)

      qc <- MoPAC::quality.control(id=paths$out,
                                   counts=counts1()$Counts,
                                   annotation=rtable1,
                                   pseudocount=pseudocount,
                                   report=input$qc.report,
                                   shinyF=incProgress)

      incProgress(0.1,"Quality Control","Finalizing...")

      sgRNA_reps_mat <- as.matrix(qc$sgRNA_lfc_reps[,sapply(qc$sgRNA_lfc_reps,is.numeric),drop=FALSE])
      rownames(sgRNA_reps_mat) <- qc$sgRNA_lfc_reps$sgRNA
      if(ncol(sgRNA_reps_mat)>2){
        pca1_sample <- stats::prcomp(t(sgRNA_reps_mat),scale=F)
        pca1 <- data.frame(Z1=pca1_sample$x[,1],Z2=pca1_sample$x[,2],Z3=pca1_sample$x[,3],Sample=colnames(sgRNA_reps_mat))
        pve1_sample <- 100*pca1_sample$sdev^2/sum(pca1_sample$sdev^2)
      } else pca1 <- pve1_sample <- NULL

      gene_reps_mat <- as.matrix(qc$gene_lfc_reps[,sapply(qc$gene_lfc_reps,is.numeric),drop=FALSE])
      rownames(gene_reps_mat) <- qc$gene_lfc_reps$Gene
      if(ncol(sgRNA_reps_mat)>2){
        pca2_sample <- stats::prcomp(t(sgRNA_reps_mat),scale=F)
        pca2 <- data.frame(Z1=pca2_sample$x[,1],Z2=pca2_sample$x[,2],Z3=pca2_sample$x[,3],Sample=colnames(sgRNA_reps_mat))
        pve2_sample <- 100*pca2_sample$sdev^2/sum(pca2_sample$sdev^2)
      } else pca2 <- pve2_sample <- NULL

      updateSelectInput(session, "use_reps1",
                        label = "Choose samples to proceed with the analysis:",
                        choices = qc$Annotation$Sample,
                        selected = 1)
      conditions <- unique(qc$Annotation1$Condition)
      updateSelectInput(session, "lfc.condition",
                        label = "Condition:",
                        choices = c("All conditions",conditions),
                        selected = "All conditions"
      )
      updateSelectInput(session, "qc.condition",
                        label = "Visualize replicates of condition:",
                        choices = conditions,
                        selected = conditions[1]
      )
    })
    return(list(QC=qc,sgRNA_reps_mat=sgRNA_reps_mat,gene_reps_mat=gene_reps_mat,PCA1=pca1,PVE1=pve1_sample,PCA2=pca2,PVE2=pve2_sample))
  })

  output$text_anno2 <- renderUI({
    if(!input$save_anno1) {
      return(HTML(paste0("Please fill out the table below and then click on the red button when finished in order to compute the sgRNA fold changes. For each condition, please write its corresponding 'Day 0' or plasmid in the 'Control' column. For each 'Day 0' or plasmid condition, leave the 'Control' value blank. <br/>")))
    } else if(!is.null(qc()$QC)) return(HTML(paste0("sgRNA log-fold-changes saved at ",paths$out,"/sgRNA_lfc.txt")))
  })

  output$anno1_table <- rhandsontable::renderRHandsontable({
    req(counts1()$Annotation)
    rtable <- counts1()$Annotation
    rownames(rtable) <- NULL
    rhandsontable::rhandsontable(rtable) # converts the R dataframe to rhandsontable object
  })

  output$variance <- renderPlot({
    req(qc()$QC$Pseudo)
    if(nrow(qc()$QC$Pseudo$Variance)==1) return(NULL)
    plot.pseudocount(qc()$QC$Pseudo)
  })

  observeEvent(input$log.lfc,{
    if(input$log.lfc=="Log-fold-change"){
      conditions <- unique(qc()$QC$Annotation1$Condition)
    } else conditions <- unique(qc()$QC$Annotation$Condition)
    updateSelectInput(session, "lfc.condition",
                      label = "Condition:",
                      choices = c("All conditions",conditions),
                      selected = "All conditions")
  })

  output$folds <- renderPlot({
    req(qc(),input$lfc.condition,input$folds.level)
    annotation2 <- qc()$QC$Annotation

    if(input$lfc.condition=="All conditions") {
      annotation2 <- annotation2[,c("Sample","Condition","Replicate")]
    } else annotation2 <- annotation2[annotation2$Condition==input$lfc.condition,c("Sample","Condition","Replicate")]

    # Dataset selection:
    loglevel = "Fold"
    if(input$folds.level=="sgRNA"){
      if(input$log.lfc=="Log-fold-change") {
        dataset <- qc()$QC$sgRNA_lfc_reps
      } else {
        dataset <- qc()$QC$sgRNA_log_reps
        loglevel = "Log"
      }
      if(!"Category"%in%colnames(dataset)) dataset$Category <- "All"
      dataset <- reshape2::melt(dataset,id.vars=c("sgRNA","Gene","Category"),value.name="Value",variable.name="Sample")
    } else {
      if(input$log.lfc=="Log-fold-change") {
        dataset <- qc()$QC$gene_lfc_reps
      } else {
        dataset <- qc()$QC$gene_log_reps
        loglevel = "Log"
      }
      if(!"Category"%in%colnames(dataset)) dataset$Category <- "All"
      dataset <- reshape2::melt(dataset,id.vars=c("Gene","Category"),value.name="Value",variable.name="Sample")
    }
    dataset <- merge(dataset,annotation2)
    plot.folds(dataset, loglevel=loglevel)
  }, height = function() { session$clientData$output_folds_width * 0.5 })

  output$ssmd <- renderPlot({
    req(qc())
    if(!"Category"%in%colnames(qc()$QC$sgRNA_lfc_reps)) return(NULL)
    plot.ssmd(qc()$QC$sgRNA_lfc_reps)
  })

  observeEvent(input$corr_pca,{
    if(input$corr_pca=="Correlation") {
      updateRadioButtons(session,"corr_pca_method",label="Method:",choices=c("pearson","spearman"),selected="pearson")
    } else updateRadioButtons(session,"corr_pca_method",label="Components:",choices=c("2nd vs 1st","3rd vs 1st"),selected="2nd vs 1st")
  })

  output$correlation <- renderPlot({
    req(input$corr_pca=="Correlation")
    req(input$corr_pca_method%in%c("pearson","spearman"))
    validate(need(ncol(qc()$sgRNA_reps_mat)>1,"At least two samples are required in order to compute correlations."))
    if(input$corr.level=="sgRNA") {
      dataset <- qc()$sgRNA_reps_mat
    } else dataset <- qc()$gene_reps_mat
    plot.correlation(dataset,method=input$corr_pca_method)
  }, height = function() { session$clientData$output_correlation_width })

  output$pca1 <- renderPlot({
    req(input$corr_pca=="Principal component analysis")
    validate(need(ncol(qc()$sgRNA_reps_mat)>2,"More samples are required in order to compute PCA"))
    if(input$corr_pca_method=="2nd vs 1st")  {
      Y <- "Z2"
      Ylab <- "Second principal component"
    } else if(input$corr_pca_method=="3rd vs 1st")  {
      Y <- "Z3"
      Ylab <- "Third principal component"
    } else return(NULL)
    if(input$corr.level=="sgRNA") {
      dataset <- qc()$PCA1
    } else dataset <- qc()$PCA2
    ggplot2::ggplot(dataset, ggplot2::aes_string("Z1",Y,label="Sample")) +
      ggrepel::geom_text_repel() +
      ggplot2::geom_point() +
      ggplot2::ggtitle(input$corr_pca_method) +
      ggplot2::xlab("First principal component") +
      ggplot2::ylab(Ylab) +
      ggplot2::theme(
        title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
        legend.text = ggplot2::element_text(size=15,color="black"),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15),
        panel.background=ggplot2::element_blank(),
        panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
        legend.title=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank(),
        strip.background=ggplot2::element_rect(fill="#0033cc"),
        strip.text = ggplot2::element_text(face="bold",color="white"),
        axis.text.x=ggplot2::element_blank())
  },height = function() { session$clientData$output_pca1_width })

  output$pve <- renderPlot({
    req(input$corr_pca=="Principal component analysis")
    validate(need(ncol(qc()$sgRNA_reps_mat)>2,"More samples are required in order to compute PCA"))
    if(input$corr.level=="sgRNA") {
      pve_sample <- qc()$PVE1
    } else pve_sample <- qc()$PVE2
    pve <- data.frame(PC=seq(1:length(pve_sample)),pve_sample=pve_sample)
    cpve <- data.frame(PC=seq(1:length(pve_sample)),Cpve_sample=cumsum(pve_sample))
    p1 <- ggplot2::ggplot(pve, ggplot2::aes(PC,pve_sample)) +
      ggplot2::geom_point(col="darkblue",shape="O",size=5) +
      ggplot2::geom_line(col="blue") +
      ggplot2::xlab("Principal Component") +
      ggplot2::ylab("Proportion of Variance Explained") +
      ggplot2::theme(
        title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
        legend.text = ggplot2::element_text(size=15,color="black"),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15),
        panel.background=ggplot2::element_blank(),
        panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
        legend.title=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank(),
        strip.background=ggplot2::element_rect(fill="#0033cc"),
        strip.text = ggplot2::element_text(face="bold",color="white"),
        axis.text.x=ggplot2::element_blank())
    p2 <- ggplot2::ggplot(cpve, ggplot2::aes(PC,Cpve_sample)) +
      ggplot2::geom_line(col="blue") +
      ggplot2::geom_point(col="darkblue",shape="O",size=5) +
      ggplot2::xlab("Principal Component") +
      ggplot2::ylab("Cumulative PVE") +
      ggplot2::theme(
        title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
        legend.text = ggplot2::element_text(size=15,color="black"),
        legend.position="top",
        text=ggplot2::element_text(face="bold",color="black",size=15),
        panel.background=ggplot2::element_blank(),
        panel.border=ggplot2::element_rect(colour="black",fill=NA,size=1),
        legend.title=ggplot2::element_blank(),
        axis.ticks.x=ggplot2::element_blank(),
        strip.background=ggplot2::element_rect(fill="#0033cc"),
        strip.text = ggplot2::element_text(face="bold",color="white"),
        axis.text.x=ggplot2::element_blank())
    gridExtra::grid.arrange(p1,p2,nrow=1)
  })

  observeEvent(input$reproducibility.stage,{
    if(input$reproducibility.stage=="Log-fold-change"){
      conditions <- unique(qc()$QC$Annotation1$Condition)
    } else conditions <- unique(qc()$QC$Annotation$Condition)
    updateSelectInput(session, "qc.condition",
                      label = "Visualize replicates of condition:",
                      choices = conditions,
                      selected = conditions[1])
  })

  output$reproducibility <- renderPlot({
    req(input$qc.condition)
    validate(need(ncol(qc()$sgRNA_reps_mat)>2,"More samples are required in order to show reproducibility."))
    if(input$reproducibility.stage=="Log-fold-change") {
      annotation2 <- qc()$QC$Annotation1[qc()$QC$Annotation1$Condition==input$qc.condition,]
      grp <- split(annotation2$Sample,annotation2[c("Condition")],drop=TRUE)
    } else {
      annotation2 <- qc()$QC$Annotation[qc()$QC$Annotation$Condition==input$qc.condition,]
      grp <- split(annotation2$Sample,annotation2[c("Condition")],drop=TRUE)
    }
    if(input$reproducibility.level=="sgRNA"){
      if(input$reproducibility.stage=="Log-fold-change") {
        cols <- which(colnames(qc()$QC$sgRNA_lfc_reps)%in%c("Category",as.character(grp[[1]])))
        dataset <- qc()$QC$sgRNA_lfc_reps[,cols,drop=F]
      } else {
        cols <- which(colnames(qc()$QC$sgRNA_log_reps)%in%c("Category",as.character(grp[[1]])))
        dataset <- qc()$QC$sgRNA_log_reps[,cols,drop=F]
      }
    } else {
      if(input$reproducibility.stage=="Log-fold-change") {
        cols <- which(colnames(qc()$QC$gene_lfc_reps)%in%c("Category",as.character(grp[[1]])))
        dataset <- qc()$QC$gene_lfc_reps[,cols,drop=F]
      } else {
        cols <- which(colnames(qc()$QC$gene_log_reps)%in%c("Category",as.character(grp[[1]])))
        dataset <- qc()$QC$gene_log_reps[,cols,drop=F]
      }
    }
    plot.reproducibility(dataset,names(grp))
  },height = function() { session$clientData$output_reproducibility_width })

  output$data_qc <- DT::renderDataTable({ #DT is necessary to make the row names appear
    req(qc())
    if(input$qc.level=="sgRNA") {
      if(input$qc.stage=="Log-fold-change") {
        if(input$qc.averaged==FALSE) return(qc()$QC$sgRNA_lfc_reps) else return(qc()$QC$sgRNA_lfc)
      } else {
        if(input$qc.averaged==FALSE) return(qc()$QC$sgRNA_log_reps) else return(qc()$QC$sgRNA_log)
      }
    } else {
      if(input$qc.stage=="Log-fold-change") {
        if(input$qc.averaged==FALSE) return(qc()$QC$gene_lfc_reps) else return(qc()$QC$gene_lfc)
      } else {
        if(input$qc.averaged==FALSE) return(qc()$QC$gene_log_reps) else return(qc()$QC$gene_log)
      }
    }
  }, options = list(pageLength = 30), rownames=FALSE)

  ### 2-Tail RRA ### -----------------------------------------------------------

  ranges <- reactiveValues(x=NULL,y=NULL,xd=NULL,yd=NULL,xr=NULL,yr=NULL)

  rra <- eventReactive(input$button_ea,{
    req(qc())

    withProgress(message='2-tail RRA',value=0,detail="Initializing...",{
      validate(need(paths$out!=0,"Please enter an output directory."))

      rra <- MoPAC::RRA.2tail(id=paths$out,
                              conditions=NULL,
                              pvalue=input$pvalue,
                              fraction=input$fraction,
                              shinyF=incProgress)
      conditions1 <- names(rra$p.value)
      updateSelectInput(session, "rra.condition",
                        label = "Visualize essentiality p values of condition:",
                        choices = conditions1,
                        selected = conditions1[1])
      updateSelectizeInput(session=session,inputId='rragene',choices=rra$p.value[[1]]$Gene,selected=NULL,server=TRUE)
    })
    return(rra)
  })

  plot.rra <- function() {
    pvalue <- rra()$p.value[[input$rra.condition]]

    pvalue$p <- -log10(pvalue$p)
    if(!"Category"%in%colnames(pvalue)) pvalue$Category <- "Other"
    pvalue$Category[pvalue$Category=="Essential"] <- "Nonfiltered-essential"
    pvalue$Category[pvalue$Category=="Nonessential"] <- "Nonfiltered-nonessential"
    pvalue$Category[pvalue$Gene%in%rra()$Essential] <- "Filtered-essential"
    pvalue$Category[pvalue$Gene%in%rra()$Nonessential] <- "Filtered-nonessential"
    colnames(pvalue)[which(colnames(pvalue)==input$rra.condition)] <- "X"

    # pvalue$Category <- factor(pvalue$Category,levels=c("Filtered-essential","Nonfiltered-essential","Filtered-nonessential","Nonfiltered-nonessential","Nontargeting","Other"))
    ggplot2::ggplot(pvalue,ggplot2::aes_string("X","p",col="Category")) +
      ggplot2::ggtitle(paste0(input$rra.condition," (select area)")) +
      ggplot2::geom_point(data=pvalue[pvalue$Category=="Other",],size=3,col="black",pch=1) +
      ggplot2::geom_point(data=pvalue[pvalue$Category!="Other",],size=3,pch=16) +
      ggplot2::geom_point(data=pvalue[pvalue$Gene%in%input$rragene,],size=4,pch=21,col="black",fill="blue") +
      ggplot2::geom_point(data=pvalue[pvalue$Gene==flags$click1,],size=5,pch=21,col="black",fill="goldenrod") +
      ggplot2::scale_color_manual(name="Category",values=c(Other="black",`Filtered-essential`="red",`Filtered-nonessential`="green",`Nonfiltered-essential`="darkred",`Nonfiltered-nonessential`="darkgreen",Nontargeting="darkgreen")) +
      ggplot2::xlab("Log-Fold-Change") + ggplot2::ylab("-log10(P-value)") +
      ggplot2::guides(col=ggplot2::guide_legend(nrow=2,byrow=FALSE)) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size=15,color="black"),
            axis.text = ggplot2::element_text(size=15,color="black"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
  }

  observeEvent(c(rra(),input$rra.condition,input$rragene,input$clickrrazoom$x), {
    output$rra <- renderPlot({
      req(rra(),input$rra.condition)
      isolate(plot.rra())
    },height = function() {session$clientData$output_rra_width })
  })

  # After hovering on the plot:
  output$plotrrainfo <- renderUI({
    req(rra(),input$rra.condition)
    pvalue <- rra()$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    hover <- input$rra_hover
    dataset <- pvalue
    point <- nearPoints(dataset,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=colnames(pvalue)[which(colnames(pvalue)==input$rra.condition)],yvar="p")
    if(nrow(point)==0) return(NULL)
    info <- dataset[dataset$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220

  # After clicking on the plot:
  flags <- reactiveValues(click1 = 0, click2 = 0) # to control the appearance/dissapearance of widgets

  observeEvent(input$clickrrazoom$x, {
    pvalue <- rra()$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    if(!"Category"%in%colnames(pvalue)) pvalue$Category <- "Other"
    pvalue$Category[pvalue$Category=="Essential"] <- "Unfiltered-essential"
    pvalue$Category[pvalue$Category=="Nonessential"] <- "Unfiltered-nonessential"
    pvalue$Category[pvalue$Gene%in%rra()$Essential] <- "Filtered-essential"
    pvalue$Category[pvalue$Gene%in%rra()$Nonessential] <- "Filtered-nonessential"

    point <- nearPoints(pvalue,input$clickrrazoom,threshold=10,maxpoints=1,xvar=colnames(pvalue)[which(colnames(pvalue)==input$rra.condition)],yvar="p")
    flags$click1 <- point$Gene
    output$inforra <- renderUI({
      req(rra(),flags$click1)
      gene <- flags$click1 #the sample on which the user clicked
      info <- paste0("<font color=#0000FF><b>Gene name: ",gene,"</b></font>")
      info <- paste(info,paste0("Category: ",point$Category),sep="<br/>")
      for(condition in names(rra()$p.value)){
        info <- paste(info,paste0("<br/><font color=#0000FF><b>Condition: ",condition,"</b></font>"),sep="<br/>")
        info <- paste(info,paste0("Fold-change = ",rra()$p.value[[condition]][rra()$p.value[[condition]]$Gene==flags$click1,3]),sep="<br/>")
        info <- paste(info,paste0("P-value = ",rra()$p.value[[condition]][rra()$p.value[[condition]]$Gene==flags$click1,"p"]),sep="<br/>")
        info <- paste(info,paste0("FDR = ",rra()$p.value[[condition]][rra()$p.value[[condition]]$Gene==flags$click1,"FDR"]),sep="<br/>")
        info <- paste(info,paste0("Depletion rank = ",rra()$p.value[[condition]][rra()$p.value[[condition]]$Gene==flags$click1,"depleted"]),sep="<br/>")
        info <- paste(info,paste0("Enrichment rank = ",rra()$p.value[[condition]][rra()$p.value[[condition]]$Gene==flags$click1,"enriched"]),sep="<br/>")
      }
      return(HTML(info))
    })
  })

  # Zoom in on Gene distributions:
  output$rrazoom <- renderPlot({
    req(rra(),input$rra.condition)
    pvalue <- rra()$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    if(!"Category"%in%colnames(pvalue)) pvalue$Category <- "Other"
    pvalue$Category[pvalue$Category=="Essential"] <- "Nonfiltered-essential"
    pvalue$Category[pvalue$Category=="Nonessential"] <- "Nonfiltered-nonessential"
    pvalue$Category[pvalue$Gene%in%rra()$Essential] <- "Filtered-essential"
    pvalue$Category[pvalue$Gene%in%rra()$Nonessential] <- "Filtered-nonessential"
    colnames(pvalue)[which(colnames(pvalue)==input$rra.condition)] <- "X"
    # pvalue$Category <- factor(pvalue$Category,levels=c("Filtered-essential","Nonfiltered-essential","Filtered-nonessential","Nonfiltered-nonessential","Nontargeting","Other"))
    p <- ggplot2::ggplot(pvalue,ggplot2::aes_string("X","p",col="Category")) +
      ggplot2::ggtitle(paste0(input$rra.condition," (click on a point)")) +
      ggplot2::scale_color_manual(name="Category",values=c(Other="black",`Filtered-essential`="red",`Filtered-nonessential`="green",`Nonfiltered-essential`="darkred",`Nonfiltered-nonessential`="darkgreen",Nontargeting="darkgreen")) +
      ggplot2::xlab("Log-Fold-Change") + ggplot2::ylab("-log10(P-value)") +
      ggplot2::guides(col=ggplot2::guide_legend(nrow=2,byrow=FALSE)) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size=15,color="black"),
            axis.text = ggplot2::element_text(size=15,color="black"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    if(!is.null(input$plotrra_brush)) {
      brush <- input$plotrra_brush
      pvalue1 <- pvalue[pvalue$X>brush$xmin & pvalue$X<brush$xmax & pvalue$p>brush$ymin & pvalue$p<brush$ymax,]
      p <- p + ggplot2::geom_point(data=pvalue1[pvalue1$Category=="Other",],size=3,col="black",pch=1) +
        ggplot2::geom_point(data=pvalue1[pvalue1$Category!="Other",],size=3,pch=16) +
        ggplot2::geom_point(data=pvalue1[pvalue1$Gene%in%input$rragene,],size=4,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=pvalue1[pvalue1$Gene==flags$click1,],size=5,pch=21,col="black",fill="goldenrod")
    } else {
      p <- p + ggplot2::geom_point(data=pvalue[pvalue$Category=="Other",],size=3,col="black",pch=1) +
        ggplot2::geom_point(data=pvalue[pvalue$Category!="Other",],size=3,pch=16) +
        ggplot2::geom_point(data=pvalue[pvalue$Gene%in%input$rragene,],size=4,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=pvalue[pvalue$Gene==flags$click1,],size=5,pch=21,col="black",fill="goldenrod")
    }
    return(p)
  },height = function() { session$clientData$output_rrazoom_width })

  # After hovering on the plot:
  output$plotrrazoominfo <- renderUI({
    req(rra(),input$rra.condition)
    pvalue <- rra()$p.value[[input$rra.condition]]
    pvalue$p <- -log10(pvalue$p)
    hover <- input$rrazoom_hover
    dataset <- pvalue
    point <- nearPoints(dataset,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=colnames(dataset)[which(colnames(pvalue)==input$rra.condition)],yvar="p")
    if(nrow(point)==0) return(NULL)
    info <- dataset[dataset$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220


  ### ESSENTIALITY ### ------------------------------------------------------------

  scores <- eventReactive(rra(),{
    req(rra())
    withProgress(message='Essentiality analysis',value=0,detail="Initializing...",{

      scores <- MoPAC::analyze.essentiality(id=paths$out,
                                            empirical.weights=ifelse(input$empirical.weight=="Pre-trained",T,F),
                                            read.filtered=input$read.filtered,
                                            filtered.reps=input$use_reps1,
                                            shinyF=incProgress)

      # rra1 <- lapply(rra()$p.value, function(x) {
      #   if(!"Category"%in%colnames(x)) {
      #     condition <- colnames(x)[2]
      #     x <- data.frame(Gene=x$Gene,Category="Other",p=x$p,stringsAsFactors=F)
      #   } else {
      #     condition <- colnames(x)[3]
      #     x <- x[,c("Gene","Category","p")]
      #   }
      #   x <- merge(x,scores()$scored[,c("Gene",condition)])
      # })

      conditions <- colnames(scores$scored)[-which(colnames(scores$scored)%in%c("Gene","Category"))]
      updateSelectInput(session, "ea.condition",
                        label = "Condition:",
                        choices = c("All conditions",conditions),
                        selected = "All conditions"
      )
      updateSelectInput(session, "condition1",
                        label = "First condition:",
                        choices = conditions,
                        selected = conditions[1]
      )
      updateSelectInput(session, "condition2",
                        label = "Second condition:",
                        choices = conditions,
                        selected = conditions[2]
      )
      updateSelectizeInput(session=session,inputId='dagene',choices=scores$scored$Gene,selected=NULL,server=TRUE)
    })
    return(scores)
  })

  output$scores <- renderPlot({
    req(scores(),input$ea.condition)
    if(input$ea.condition=="All conditions") {
      annotation2 <- qc()$QC$Annotation1
    } else annotation2 <- qc()$QC$Annotation1[qc()$QC$Annotation1$Condition==input$ea.condition,]
    plot.scores(scores()$scored_reps, annotation2)
  }, height = function() { session$clientData$output_scores_width * 0.5 })

  output$data_ea <- DT::renderDataTable({ #DT is necessary to make the row names appear
    req(scores())
    if(input$ea.averaged==FALSE) {
      return(scores()$scored_reps)
    } else return(scores()$scored)
  }, options = list(pageLength = 30), rownames=FALSE)

  ### Differential essentiality ### ---------------------------------------------

  diff <- eventReactive(input$button_da,{
    req(scores())
    withProgress(message='Differential analysis',value=0,detail="Initializing...",{
      if(input$normalization=="Controls"){
        read.unenriched <- FALSE
      } else read.unenriched <- TRUE
      diff <- MoPAC::analyze.differential(id=paths$out,
                                           condition1=input$condition1,
                                           condition2=input$condition2,
                                           read.unenriched=read.unenriched,
                                           report=input$da.report)$Gene
      if(!"Category"%in%colnames(diff))
        diff <- data.frame(Gene=diff$Gene,Category="Other",diff[,-1],check.names=FALSE,stringsAsFactors=FALSE)
      if(!"P(replicates)"%in%colnames(diff))
        diff$`P(replicates)` <- 2*pnorm(-abs(diff$Z))
      # diff$logP <- -log10(diff$`P(sgRNAs)`)
    })
    diff_out <- basename(list.files(path=paths$out,pattern="_differential",full.names=TRUE))
    comparisons <- substring(diff_out,1,regexpr("_differential",diff_out)-1)
    updateSelectInput(session, "comparisons",
                      label = "Condition comparisons:",
                      choices = comparisons,
                      selected = comparisons[1]
    )
    return(diff)
  })

  # Plotting distributions:
  output$plot1 <- renderPlot({
    req(scores(),diff())
    ggplot2::ggplot(diff(),ggplot2::aes_string(colnames(diff())[3],colnames(diff())[4],col="Category")) +
      ggplot2::ggtitle("Essentiality (select area)") +
      ggplot2::geom_abline(slope=1,size=2,col="darkgray") +
      ggplot2::geom_point(data=diff()[diff()$Category=="Other",],size=3,col="black",pch=1) +
      ggplot2::geom_point(data=diff()[diff()$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
      ggplot2::geom_point(data=diff()[diff()$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
      ggplot2::geom_point(data=diff()[diff()$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod") +
      ggplot2::scale_color_manual(name="Category",values=c(Other="black",`Essential`="red",`Nonessential`="green",Nontargeting="darkgreen")) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size=15,color="black"),
            axis.text = ggplot2::element_text(size=15,color="black"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
  },height = function() {session$clientData$output_plot1_width })

  # After hovering on the plot:
  output$plot1info <- renderUI({
    req(diff())
    hover <- input$plot1_hover
    point <- nearPoints(diff(),hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=colnames(diff())[3],yvar=colnames(diff())[4])
    if(nrow(point)==0) return(NULL)
    info <- diff()[diff()$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
                    # paste0("<b>Category: </b>",info$Category),
                    # paste0("<b>Z score: </b>",info$Z),
                    # paste0("<b>P value: </b>",info$P),
                    # paste0("<b>FDR: </b>",info$FDR),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220


  observeEvent(input$click2$x, {
    req(diff())
    point <- nearPoints(diff(),input$click2,threshold=5,maxpoints=1,addDist=TRUE,xvar=colnames(diff())[3],yvar=colnames(diff())[4])
    flags$click2 <- point$Gene
    output$infoscore <- renderUI({
      req(diff(),flags$click2)
      gene <- flags$click2 #the sample on which the user clicked
      info <- paste0("<br/><font color=#0000FF><b>Gene name: ",gene,"</b></font>")
      info <- paste(info,paste0("Category: ",point$Category),sep="<br/>")
      info <- paste(info,paste0("Z score: ",formatC(point$Z,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("P value: ",formatC(point$P,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("FDR: ",formatC(point$FDR,format="e",digits=2)),sep="<br/>")
      return(HTML(info))
    })
  })

  # Zoom in on Gene distributions:
  output$plot2 <- renderPlot({
    req(scores(),diff())
    p <- ggplot2::ggplot(diff(),ggplot2::aes_string(colnames(diff())[3],colnames(diff())[4],col="Category")) +
      ggplot2::ggtitle("Essentiality (click on a point)") +
      ggplot2::geom_abline(slope=1,size=2,col="darkgray") +
      ggplot2::scale_color_manual(name="Category",values=c(Other="black",`Essential`="red",`Nonessential`="green",Nontargeting="darkgreen")) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size=15,color="black"),
            axis.text = ggplot2::element_text(size=15,color="black"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    if(!is.null(input$plot1_brush)) {
      brush <- input$plot1_brush
      dataset1 <- diff()[diff()[,3]>brush$xmin & diff()[,3]<brush$xmax & diff()[,4]>brush$ymin & diff()[,4]<brush$ymax,]
      p <- p + ggplot2::geom_point(data=dataset1[dataset1$Category=="Other",],size=3,col="black",pch=1) +
        ggplot2::geom_point(data=dataset1[dataset1$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
        ggplot2::geom_point(data=dataset1[dataset1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=dataset1[dataset1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    } else {
      p <- p + ggplot2::geom_point(data=diff()[diff()$Category=="Other",],size=3,col="black",pch=1) +
        ggplot2::geom_point(data=diff()[diff()$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
        ggplot2::geom_point(data=diff()[diff()$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=diff()[diff()$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    }
    return(p)
  },height = function() { session$clientData$output_plot2_width })

  # After hovering on the plot:
  output$plot2info <- renderUI({
    req(diff())
    hover <- input$plot2_hover
    point <- nearPoints(diff(),hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=colnames(diff())[3],yvar=colnames(diff())[4])
    if(nrow(point)==0) return(NULL)
    info <- diff()[diff()$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
                    # paste0("<b>Category: </b>",info$Category),
                    # paste0("<b>Z score: </b>",info$Z),
                    # paste0("<b>P value: </b>",info$P),
                    # paste0("<b>FDR: </b>",info$FDR),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220

  output$weights <- renderPlot({
    req(scores())
    if(is.null(scores()$weights)) return(NULL)
    weights <- scores()$weights
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
        axis.text.x=ggplot2::element_text(size=15),#,angle=90,hjust=1,vjust=0.5),
        axis.title=ggplot2::element_text(size=15,face="bold",color="black"),
        axis.ticks.length=grid::unit(0.15,"cm"))
  })

  output$scaling <- renderPlot({
    req(scores())
    if(is.null(scores()$scaling)) return(NULL)
    ggplot2::ggplot(scores()$scaling, ggplot2::aes_string("Sample","Score")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::ylab("Scaling factor") +
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
  })

  output$shifting <- renderPlot({
    req(scores())
    if(is.null(scores()$shifting)) return(NULL)
    ggplot2::ggplot(scores()$shifting,ggplot2::aes_string("Sample","Score")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::ylab("Shifting factor") +
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
  })

  ### Volcano plots ###

  output$plot3 <- renderPlot({
    req(scores(),diff())
    tmp <- diff()
    if(input$p.volcano=="sgRNAs") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    ggplot2::ggplot(tmp,ggplot2::aes(Z,logP,col=Category)) +
      ggplot2::ggtitle("Differential essentiality (select area)") +
      ggplot2::geom_point(data=tmp[tmp$Category=="Other",],size=3,col="black",pch=1) +
      ggplot2::geom_point(data=tmp[tmp$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
      ggplot2::geom_point(data=tmp[tmp$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
      ggplot2::geom_point(data=tmp[tmp$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod") +
      ggplot2::xlab("Z score") + ggplot2::ylab("-log10(P-value)") +
      ggplot2::scale_color_manual(name="Category",values=c(Other="black",`Essential`="red",`Nonessential`="green",Nontargeting="darkgreen")) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size=15,color="black"),
            axis.text = ggplot2::element_text(size=15,color="black"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
  },height = function() {session$clientData$output_plot3_width })

  # After hovering on the plot:
  output$plot3info <- renderUI({
    req(diff())
    hover <- input$plot3_hover
    tmp <- diff()
    if(input$p.volcano=="sgRNAs") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    point <- nearPoints(tmp,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar="Z",yvar="logP")
    if(nrow(point)==0) return(NULL)
    info <- tmp[tmp$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
                    # paste0("<b>Category: </b>",info$Category),
                    # paste0("<b>Z score: </b>",info$Z),
                    # paste0("<b>P value: </b>",info$P),
                    # paste0("<b>FDR: </b>",info$FDR),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220

  # Zoom in on Gene distributions:
  output$plot4 <- renderPlot({
    req(scores(),diff())
    tmp <- diff()
    if(input$p.volcano=="sgRNAs") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    p <- ggplot2::ggplot(tmp,ggplot2::aes(Z,logP,col=Category)) +
      ggplot2::ggtitle("Differential essentiality (click on a point)") +
      ggplot2::xlab("Z score") + ggplot2::ylab("-log10(P-value)") +
      ggplot2::scale_color_manual(name="Category",values=c(Other="black",`Essential`="red",`Nonessential`="green",Nontargeting="darkgreen")) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size=15,color="black"),
            axis.text = ggplot2::element_text(size=15,color="black"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    if(!is.null(input$plot3_brush)) {
      brush <- input$plot3_brush
      dataset1 <- tmp[tmp$Z>brush$xmin & tmp$Z<brush$xmax & tmp$logP>brush$ymin & tmp$logP<brush$ymax,]
      p <- p + ggplot2::geom_point(data=dataset1[dataset1$Category=="Other",],size=3,col="black",pch=1) +
        ggplot2::geom_point(data=dataset1[dataset1$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
        ggplot2::geom_point(data=dataset1[dataset1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=dataset1[dataset1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    } else {
      p <- p + ggplot2::geom_point(data=tmp[tmp$Category=="Other",],size=3,col="black",pch=1) +
        ggplot2::geom_point(data=tmp[tmp$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
        ggplot2::geom_point(data=tmp[tmp$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=tmp[tmp$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    }
    return(p)
  },height = function() { session$clientData$output_plot3_width })

  # After hovering on the plot:
  output$plot4info <- renderUI({
    req(diff())
    hover <- input$plot4_hover
    tmp <- diff()
    if(input$p.volcano=="sgRNAs") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    point <- nearPoints(tmp,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar="Z",yvar="logP")
    if(nrow(point)==0) return(NULL)
    info <- tmp[tmp$Gene==point$Gene,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),sep="<br/>")
                    # paste0("<b>Category: </b>",info$Category),
                    # paste0("<b>Z score: </b>",info$Z),
                    # paste0("<b>P value: </b>",info$P),
                    # paste0("<b>FDR: </b>",info$FDR),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220

  observeEvent(input$click4$x, {
    req(diff())
    tmp <- diff()
    if(input$p.volcano=="sgRNAs") {
      tmp$logP <- -log10(tmp$`P(sgRNAs)`)
    } else tmp$logP <- -log10(tmp$`P(replicates)`)
    point <- nearPoints(tmp,input$click4,threshold=5,maxpoints=1,addDist=TRUE,xvar="Z",yvar="logP")
    point1 <- tmp[tmp$Gene==point$Gene,]
    flags$click2 <- point$Gene
    output$infoscore <- renderUI({
      req(tmp,flags$click2)
      gene <- flags$click2 #the sample on which the user clicked
      info <- paste0("<br/><font color=#0000FF><b>Gene name: ",gene,"</b></font>")
      info <- paste(info,paste0("Category: ",point1$Category),sep="<br/>")
      info <- paste(info,paste0("Z score: ",formatC(point1$Z,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("<br/>",colnames(tmp)[3]," score: ",formatC(as.numeric(point1[3]),format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("P(-): ",formatC(point1$`P(-)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("FDR(-): ",formatC(point1$`FDR(-)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("<br/>",colnames(tmp)[4]," score: ",formatC(as.numeric(point1[4]),format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("P(+): ",formatC(point1$`P(+)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("FDR(+): ",formatC(point1$`FDR(+)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("<br/>Replicate P value: ",formatC(point1$`P(replicates)`,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("<br/>sgRNA P value: ",formatC(point1$`P(sgRNAs)`,format="e",digits=2)),sep="<br/>")
      return(HTML(info))
    })
  })

  ### LOG-FOLD-CHANGES

  output$plot5 <- renderPlot({
    req(scores(),diff())
    lfc1 <- qc()$QC$sgRNA_lfc
    if(!"Category"%in%colnames(lfc1)) lfc1$Category <- "Other"
    ggplot2::ggplot(lfc1,ggplot2::aes_string(input$condition1,input$condition2,col="Category")) +
      ggplot2::ggtitle("sgRNA log-fold-changes (select area)") +
      ggplot2::geom_abline(slope=1,size=2,col="darkgray") +
      ggplot2::geom_point(data=lfc1[lfc1$Category=="Other",],size=3,col="black",pch=1) +
      ggplot2::geom_point(data=lfc1[lfc1$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
      ggplot2::geom_point(data=lfc1[lfc1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
      ggplot2::geom_point(data=lfc1[lfc1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod") +
      ggplot2::scale_color_manual(name="Category",values=c(Other="black",`Essential`="red",`Nonessential`="green",Nontargeting="darkgreen")) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size=15,color="black"),
            axis.text = ggplot2::element_text(size=15,color="black"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
  },height = function() {session$clientData$output_plot5_width })

  # After hovering on the plot:
  output$plot5info <- renderUI({
    req(scores(),diff())
    hover <- input$plot5_hover
    point <- nearPoints(qc()$QC$sgRNA_lfc,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=input$condition1,yvar=input$condition2)
    if(nrow(point)==0) return(NULL)
    info <- qc()$QC$sgRNA_lfc[qc()$QC$sgRNA_lfc$sgRNA==point$sgRNA,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),
                    paste0("<b>sgRNA: </b>",info$sgRNA),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220

  output$plot6 <- renderPlot({
    req(scores(),diff())
    lfc1 <- qc()$QC$sgRNA_lfc
    if(!"Category"%in%colnames(lfc1)) lfc1$Category <- "Other"
    p <- ggplot2::ggplot(lfc1,ggplot2::aes_string(input$condition1,input$condition2,col="Category")) +
      ggplot2::ggtitle("sgRNA log-fold-changes (click on a point)") +
      ggplot2::geom_abline(slope=1,size=2,col="darkgray") +
      ggplot2::scale_color_manual(name="Category",values=c(Other="black",`Essential`="red",`Nonessential`="green",Nontargeting="darkgreen")) +
      ggplot2::theme(title = ggplot2::element_text(face="bold",color="deepskyblue3",size=15),
            legend.title = ggplot2::element_blank(),
            legend.position = "bottom",
            legend.text = ggplot2::element_text(size=15,color="black"),
            axis.text = ggplot2::element_text(size=15,color="black"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(colour="black",fill=NA,size=1))
    if(!is.null(input$plot5_brush)) {
      brush <- input$plot5_brush
      dataset1 <- lfc1[lfc1[,input$condition1]>brush$xmin & lfc1[,input$condition1]<brush$xmax & lfc1[,input$condition2]>brush$ymin & lfc1[,input$condition2]<brush$ymax,]
      p <- p + ggplot2::geom_point(data=dataset1[dataset1$Category=="Other",],size=3,col="black",pch=1) +
        ggplot2::geom_point(data=dataset1[dataset1$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
        ggplot2::geom_point(data=dataset1[dataset1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=dataset1[dataset1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    } else {
      p <- p + ggplot2::geom_point(data=lfc1[lfc1$Category=="Other",],size=3,col="black",pch=1) +
        ggplot2::geom_point(data=lfc1[lfc1$Category%in%c("Essential","Nonessential"),],size=3,pch=1) +
        ggplot2::geom_point(data=lfc1[lfc1$Gene%in%input$dagene,],size=6,pch=21,col="black",fill="blue") +
        ggplot2::geom_point(data=lfc1[lfc1$Gene==flags$click2,],size=6,pch=21,col="black",fill="goldenrod")
    }
    return(p)
  },height = function() {session$clientData$output_plot6_width })

  # After hovering on the plot:
  output$plot6info <- renderUI({
    req(scores(),diff())
    hover <- input$plot6_hover
    point <- nearPoints(qc()$QC$sgRNA_lfc,hover,threshold=5,maxpoints=1,addDist=TRUE,xvar=input$condition1,yvar=input$condition2)
    if(nrow(point)==0) return(NULL)
    info <- qc()$QC$sgRNA_lfc[qc()$QC$sgRNA_lfc$sgRNA==point$sgRNA,]
    # calculate point position INSIDE the image as percent of total dimensions, from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    # create style property fot tooltip, background color is set so tooltip is a bit transparent, z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ","left:",left_px+2,"px; top:",top_px+2,"px;")
    output <- paste(paste0("<b>Gene: </b>",info$Gene),
                    paste0("<b>sgRNA: </b>",info$sgRNA),sep="<br/>")
    wellPanel(style=style, HTML(output))
  }) #reference: https://gitlab.com/snippets/16220

  observeEvent(input$click6$x, {
    req(scores(),diff())
    point <- nearPoints(qc()$QC$sgRNA_lfc,input$click6,threshold=5,maxpoints=1,addDist=TRUE,xvar=input$condition1,yvar=input$condition2)
    flags$click2 <- point$Gene
    point1 <- diff()[diff()$Gene==point$Gene,]
    output$infoscore <- renderUI({
      req(diff(),flags$click2)
      gene <- flags$click2 #the sample on which the user clicked
      info <- paste0("<br/><font color=#0000FF><b>Gene name: ",gene,"</b></font>")
      if("Category"%in%colnames(diff())) info <- paste(info,paste0("Category: ",point1$Category),sep="<br/>")
      info <- paste(info,paste0("Z score: ",formatC(point1$Z,format="e",digits=2)),sep="<br/>")
      info <- paste(info,paste0("Replicate P value: ",formatC(point1$`P(replicates)`,format="e",digits=2)),sep="<br/>")
      return(HTML(info))
    })
  })

  ### TABLE OF SCORES ###

  output$table_da <- DT::renderDataTable({ #DT is necessary to make the row names appear
    req(scores(),diff())
    if(input$data_da=="Normalized pair"){
      return(diff())
      # return(diff()[,-which(colnames(diff())%in%"logP")])
    } else if(input$data_da=="Pre-normalized replicates") {
      return(scores()$scored_reps)
    } else return(scores()$scored)
  }, options = list(pageLength = 30), rownames=FALSE)

  ### NETWORK ### -------------------------------------------------------

  # This section contains internal functions of the STRINGdb R package used for network analysis.
  # Source: http://www.string-db.org
  #
  #  copyright:   Andrea Franceschini (Swiss Institute of Bioinformatics) andrea.franceschini@isb-sib.ch
  #

  stringdb_merge.with.order <- function(x,y, ..., sort=T) {
    #  copyright:   Andrea Franceschini (Swiss Institute of Bioinformatics) andrea.franceschini@isb-sib.ch
    # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
    add.id.column.to.data <- function(DATA) {
      data.frame(DATA, id... = seq_len(nrow(DATA)))
    }
    # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
    order.by.id...and.remove.it <- function(DATA) {
      # gets in a data.frame with the "id..." column.  Orders by it and returns it
      if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
      ss_r <- order(DATA$id...)
      ss_c <- colnames(DATA) != "id..."
      DATA[ss_r, ss_c]
    }
    if(sort==F){ return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
    } else {return(merge(x=x,y=y,..., sort = sort))}
  }

  stringdb_multi_map_df <- function (dfToMap, dfMap, strColsFrom, strColFromDfMap, strColToDfMap, caseSensitive = FALSE) {
    #  copyright:   Andrea Franceschini (Swiss Institute of Bioinformatics) andrea.franceschini@isb-sib.ch
    delColDf <- function (df, colName) {
      if (colName %in% names(df))
        return(df[, -which(names(df) %in% c(colName))])
      else return(df)
    }
    tempMatr = matrix(NA, length(strColsFrom), nrow(dfToMap))
    for (i in 1:length(strColsFrom)) {
      if (!caseSensitive) {
        tempMatr[i, ] = as.vector(dfToMap[, strColsFrom[i]])
        dfToMap[, strColsFrom[i]] = toupper(iconv(dfToMap[,strColsFrom[i]],"WINDOWS-1252","UTF-8"))
      }
    }
    if (!caseSensitive) {
      dfMap[,strColFromDfMap] = toupper(iconv(dfMap[,strColFromDfMap],"WINDOWS-1252","UTF-8"))
    }
    dfMap2 = unique(subset(dfMap, select = c(strColFromDfMap,strColToDfMap)))
    df2 = stringdb_merge.with.order(dfToMap, dfMap2, by.x = strColsFrom[1],by.y = strColFromDfMap, all.x = TRUE, sort = FALSE)
    if (length(strColsFrom) > 1) {
      for (i in 2:length(strColsFrom)) {
        dfna = delColDf(subset(df2, is.na(as.vector(df2[,strColToDfMap]))), strColToDfMap)
        dfgood = subset(df2, !is.na(as.vector(df2[, strColToDfMap])))
        df3 = stringdb_merge.with.order(dfna, dfMap2, by.x = strColsFrom[i],by.y = strColFromDfMap, all.x = TRUE, sort = FALSE)
        df2 = rbind(dfgood, df3)
      }
    }
    for (i in 1:length(strColsFrom)) {
      if (!caseSensitive && length(tempMatr[i, ]) == length(df2[,strColsFrom[i]]))
        df2[, strColsFrom[i]] = tempMatr[i, ]
    }
    return(df2)
  }

  # observeEvent(input$button_diff, {
  #   req(paths$out)
  #   diff_out <- basename(list.files(path=paths$out,pattern="_differential",full.names=TRUE))
  #   comparisons <- substring(diff_out,1,regexpr("_differential",diff_out)-1)
  #   updateSelectInput(session, "comparisons",
  #                     label = "Condition comparisons:",
  #                     choices = comparisons,
  #                     selected = comparisons[1]
  #   )
  # })

  net_data_in <- eventReactive(input$comparisons,{
    if(input$comparisons=="") return(NULL)
    conditions <- input$comparisons
    condition1 <- substring(conditions,1,regexpr("_",conditions)-1)
    condition2 <- substring(conditions,regexpr("_",conditions)+1)
    updateSelectInput(session, "comparisons_essential",
                      label = "Select differentially essential condition:",
                      choices = c(condition1,condition2),
                      selected = condition1
    )
    net_data_in <- utils::read.table(paste0(paths$out,"/",input$comparisons,"_differential.txt"),header=TRUE,sep="",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,row.names=NULL)

    if(!"Category"%in%colnames(net_data_in))
      net_data_in <- data.frame(Gene=net_data_in$Gene,Category="Other",net_data_in[,-1],check.names=FALSE,stringsAsFactors=FALSE)

    # if(!"P"%in%colnames(net_data_in)) {
    #   net_data_in$P <- 1
    #   net_data_in$FDR <- 1
    # }
    if(colnames(net_data_in)[3]==input$comparisons_essential) {
      net_data_in <- net_data_in[order(net_data_in$Z,decreasing=FALSE),]
      num <- sum(net_data_in$Z<-1*input$comparisons_sd)
    } else if(colnames(net_data_in)[4]==input$comparisons_essential) {
      net_data_in <- net_data_in[order(net_data_in$Z,decreasing=TRUE),]
      num <- sum(net_data_in$Z>input$comparisons_sd)
    }
    net_data_in$Rank <- 1:nrow(net_data_in)
    return(net_data_in)
  })

  output$comparisons_number <- renderUI({
    req(net_data_in())
    text <- "Number of genes = "
    num <- 0
    if(colnames(net_data_in())[3]==input$comparisons_essential) {
      num <- sum(net_data_in()$Z<(-1*input$comparisons_sd))
    } else if(colnames(net_data_in())[4]==input$comparisons_essential) {
      num <- sum(net_data_in()$Z>input$comparisons_sd)
    }
    text <- paste0(text,num)
    if(num>400) text <- paste0(text,"<br/>Only the top 400 genes will be used.")
    return(HTML(text))
  })

  network <- eventReactive(input$button_net,{
    #  copyright:   Andrea Franceschini (Swiss Institute of Bioinformatics) andrea.franceschini@isb-sib.ch
    if(input$list_source=="Custom list") {
      gene_list <- unlist(strsplit(input$gene.list,"\n"))
      if(identical(gene_list,character(0))) return(NULL)
      net_data_in <- data.frame(Gene=unlist(strsplit(input$gene.list,"\n")),stringsAsFactors=FALSE)
    } else {
      if(input$comparisons=="") return(NULL)
      if(colnames(net_data_in())[3]==input$comparisons_essential) {
        net_data_in <- net_data_in()[net_data_in()$Z<(-1*input$comparisons_sd),]
      } else if(colnames(net_data_in())[4]==input$comparisons_essential) {
        net_data_in <- net_data_in()[net_data_in()$Z>input$comparisons_sd,]
      }
    }
    if(nrow(net_data_in)>400) net_data_in <- net_data_in[1:400,,drop=FALSE]
    withProgress(message="Network analysis",value=0.1,detail="Initializing...",{
      # Search for gene aliases and assign ids:
      incProgress(0.3,"Network analysis","Searching for aliases...")
      net_data = stringdb_multi_map_df(net_data_in,MoPAC::stringdb_9606$Aliases,"Gene","alias","STRING_id")
      naDf = subset(net_data, is.na(STRING_id))
      if(nrow(naDf)>0) showNotification(paste("Warning: we couldn't map to STRING ",as.integer((nrow(naDf)/nrow(net_data))*100),"% of your identifiers",sep=""),action=NULL,duration=10,closeButton=TRUE,type="warning")
      net_data = subset(net_data, !is.na(STRING_id))
      net_data <- net_data[!duplicated(net_data$STRING_id),] #I added this line because viznetwork does not allow duplicated nodes
      # Search for protein information:
      incProgress(0.2,"Network analysis","Searching for protein information...")
      net_data = merge(net_data, MoPAC::stringdb_9606$Proteins,by.x="STRING_id",by.y="protein_external_id",all.x=TRUE,sort=FALSE)
      # Search for interactions of top genes:
      incProgress(0.2,"Network analysis","Searching for protein interactions...")
      topID <- net_data$STRING_id
      graph0 <- igraph::graph.data.frame(MoPAC::stringdb_9606$PPI[MoPAC::stringdb_9606$PPI$combined_score >= input$comparisons_tresh,], FALSE)
      graph1 <- igraph::induced.subgraph(graph0, which(igraph::V(graph0)$name %in% topID))
      # Gene clustering: #uncomment to add clustering
      # if (input$net_algorithm == "fastgreedy") fgreedy <- igraph::fastgreedy.community(graph1,merges=TRUE,modularity=TRUE)
      # if (input$net_algorithm == "walktrap") fgreedy <- igraph::walktrap.community(graph1,merges=TRUE,modularity=TRUE)
      # if (input$net_algorithm == "spinglass") fgreedy <- igraph::spinglass.community(graph1,merges=TRUE,modularity=TRUE)
      # if (input$net_algorithm == "edge.betweenness") fgreedy <- igraph::edge.betweenness.community(graph1,merges=TRUE,modularity=TRUE)
      # memb = igraph::membership(fgreedy)
      # clusters = NULL
      # for(i in 1:max(memb)) clusters[[i]] = names(igraph::membership(fgreedy)[igraph::membership(fgreedy) == i])
      # clusters <- do.call(rbind,lapply(1:length(clusters),function(x){
      #   if(length(clusters[[x]])>1) group <- x else group <- "Unlinked"
      #   return(data.frame(id=clusters[[x]],group=group))
      # }))
      # updateSelectInput(session, "net_cluster",
      #                   label="Network cluster to analyze:",
      #                   choices=c("All",unique(clusters$group)),
      #                   selected="All")
      # Creat network:
      incProgress(0.2,"Network analysis","Creating network...")
      network <- visNetwork::toVisNetworkData(graph1)
      labels <- net_data[net_data$STRING_id%in%topID,]
      labels <- merge(labels,net_data_in,sort=FALSE)
      colnames(labels)[which(colnames(labels)=="STRING_id")] <- "id"
      colnames(labels)[which(colnames(labels)=="Gene")] <- "label"
      # labels <- merge(labels,clusters) #uncomment to add clustering
      if(input$list_source!="Custom list")
        labels$title <- paste0(
          "<p>Gene = ",labels$label,
          # "<p>Rank = ",labels$Rank,
          "<br>Z score = ",signif(labels$Z,4),
           # "<br>Reproducibility FDR = ",signif(labels$FDR,4),
           "<br>Protein size = ",labels$protein_size,
           # "<br>Cluster = ",labels$group,
           "</p>")
      # labels$color.background = "blue"
      # labels$shadow.size = 1-labels$FDR
      # network$edges$title <- paste0("<p>Interaction score = ",network$edges$combined_score,"</p>")
      network$nodes <- labels
    })
    # network$nodes <- network$nodes[order(network$nodes$label),]
    return(list(topID=topID,network=network))
  })

  output$network1 <- visNetwork::renderVisNetwork({
    req(network())
    visNetwork::visNetwork(nodes=network()$network$nodes,edges=network()$network$edges,height="600px") %>%
      visNetwork::visOptions(highlightNearest=TRUE,nodesIdSelection=TRUE) %>%
      # visNetwork::visOptions(selectedBy="group",highlightNearest=TRUE,nodesIdSelection=TRUE) %>% #uncomment to add clustering
      visNetwork::visEdges(arrows=list(to=list(enabled=FALSE,scaleFactor=2)),color=list(color="black",highlight="red")) %>%
      visNetwork::visLayout(randomSeed=123) %>%
      visNetwork::visEvents(selectNode = "function(nodes) {Shiny.onInputChange('current_node_id', nodes);;}") %>%
      visNetwork::visPhysics(solver='repulsion')
        # stabilizationIterationsDone = "function(params) {network.stopSimulation();;}")
      # visNetwork::visPhysics(enabled=FALSE,stabilization=list(iterations=100))
    # visInteraction(hover = TRUE) %>%
  })
  # visNodes(color=list(border="black",highlight=list(background="blue",border="black"))) %>%

  pathways <- eventReactive(input$button_paths,{
    #  copyright:   Andrea Franceschini (Swiss Institute of Bioinformatics) andrea.franceschini@isb-sib.ch
    delColDf <- function (df, colName) {
      #  copyright:   Andrea Franceschini (Swiss Institute of Bioinformatics) andrea.franceschini@isb-sib.ch
      if (colName %in% names(df))
        return(df[, -which(names(df) %in% c(colName))])
      else return(df)
    }
    renameColDf <- function (df, colOldName, colNewName) {
      #  copyright:   Andrea Franceschini (Swiss Institute of Bioinformatics) andrea.franceschini@isb-sib.ch
      if (!(colOldName %in% names(df)))
        print(paste("ERROR: We cannot find ", colOldName, " in the data frame.", sep = ""))
      names(df)[which(names(df) == colOldName)] <- colNewName
      return(df)
    }
    req(network())
    withProgress(message="Pathway enrichment",value=0.2,detail="Initializing...",{
      if(input$net_database=="Kyoto Encyclopedia of Genes and Genomes") {
        category <- "KEGG"
      } else category <- input$net_database
      minScore <- input$net_minscore; iea <- input$net_iea
      methodMT="fdr";backgroundV=NULL
      ann = MoPAC::stringdb_9606$Annotation
      temp_category = category
      ann = subset(ann, category==temp_category)
      if(!is.null(minScore) & category!="Tissue" & category!="Disease") minScore = NULL
      if(!is.null(minScore)) {
        if(minScore<0 || minScore>5) cat("\nWARNING: minScore must be from 0 to 5 \n")
        ann = subset(ann, type >= minScore)
      }
      if(!is.null(backgroundV)) ann = subset(ann, STRING_id %in% backgroundV)
      if(!iea) ann = subset(ann, type!="IEA")
      if(nrow(ann)==0) {
        cat("\nWARNING: No annotations are present for this species with the following input parameters.\n              Please try to change species and/or parameters or try to contact our support.\n")
        return(NULL)
      }
      incProgress(0.3,"Pathway enrichment","Searching for enrichment...")
      # if(input$net_cluster=="All") { #uncomment to add clustering
        topID <- network()$topID
      # } else topID <- network()$network$nodes$id[network()$network$nodes$group==input$net_cluster] uncomment to add clustering
      if(input$net_method=="STRINGdb"){
        annHits = subset(ann, STRING_id %in% topID)
        dfcount = suppressWarnings(sqldf::sqldf("select term_id, count(STRING_id) as proteinsCount from ann group by term_id",stringsAsFactors = FALSE))
        dfcountHits = suppressWarnings(sqldf::sqldf("select term_id, count(STRING_id) as hits from annHits group by term_id",stringsAsFactors = FALSE))
        dfcountMerged = merge(dfcount,dfcountHits,by.x="term_id",by.y="term_id",all.x=TRUE)
        dfcountMerged = subset(dfcountMerged, proteinsCount <= 1500)
        dfcountMerged2 = data.frame(dfcountMerged, n = nrow(MoPAC::stringdb_9606$Proteins) - dfcountMerged$proteinsCount, k = length(unique(annHits$STRING_id)))
        dfcountMerged3 = data.frame(dfcountMerged2, pvalue = stats::phyper(dfcountMerged2$hits - 1, dfcountMerged2$proteinsCount, dfcountMerged2$n, dfcountMerged2$k, FALSE))
        dfcountMerged4 = dfcountMerged3
        if(!is.null(methodMT)) dfcountMerged4 = data.frame(dfcountMerged3, pvalue_fdr = stats::p.adjust(dfcountMerged3$pvalue, method = methodMT, n = nrow(subset(dfcountMerged3, !is.na(pvalue)))))
        incProgress(0.3,"Pathway enrichment","Attaching pathway description...")
        dfcountMerged4 = subset(dfcountMerged4, !is.na(pvalue))
        dfcountMerged4 = delColDf(dfcountMerged4, "n")
        dfcountMerged4 = delColDf(dfcountMerged4, "k")
        annDesc = MoPAC::stringdb_9606$Description
        dfcountMerged5 = plyr::arrange(merge(dfcountMerged4,annDesc,by.x="term_id",by.y="term_id",all.x=TRUE),pvalue)
        result = renameColDf(dfcountMerged5,"proteinsCount", "proteins")
      } else if(input$net_method=="NEAT (slower)") {
        incProgress(0.3,"Pathway enrichment","Running NEAT...")
        ann1 <- merge(ann,MoPAC::stringdb_9606$Description)
        ann1 <- split(ann1$STRING_id,ann1$term_description)
        graph0 <- igraph::graph.data.frame(MoPAC::stringdb_9606$PPI[MoPAC::stringdb_9606$PPI$combined_score >= input$comparisons_tresh,], FALSE)
        graph <- igraph::get.adjacency(igraph::induced.subgraph(graph0, which(igraph::V(graph0)$name %in% ann$STRING_id)))
        result = neat::neat(alist=list(Top=topID),blist=ann1,network=graph,nettype='undirected',nodes=rownames(graph),alpha=0.01)
        result <- result[order(result$pvalue),]
        rownames(result) <- 1:nrow(result)
      }
      incProgress(0.2,"Pathway enrichment","Printing...")
    })
    return(result)
  })

  output$table_net <- DT::renderDataTable({
    req(network(),pathways())
    return(pathways())
  },options=list(pageLength=30),rownames=FALSE)

  output$infonet <- renderUI({
    req(network())
    if("STRINGdb" %in% rownames(installed.packages()) == TRUE) {
      require(STRINGdb)
      string_db <- STRINGdb::STRINGdb$new(version="10",species=9606)
      link <- string_db$get_link(network()$topID)
      url <- a("Click to visualize in the STRING website", href=link, target="_blank")
      tagList(url)
    } else {
      HTML("STRINGdb installation not found.")
    }
  })

  output$infonode <- renderUI({
    req(network())
    id <- input$current_node_id$nodes
    req(id)
    info <- paste0("<br/><font color=#0000FF><b>",
                   network()$network$nodes$label[network()$network$nodes$id==id]," information:</b><br/>",
                   network()$network$nodes$annotation[network()$network$nodes$id==id],
                   "</font>")
    # info <- paste(info,paste0("FDR: ",point$FDR),sep="<br/>")
    return(HTML(info))
  })


}) #print(file=stderr(),labels)



