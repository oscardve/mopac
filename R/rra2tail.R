#' @title Gene essentiality 2-tail robust rank aggregation
#'
#' @description Performs 2-tail a-RRA analysis of gene essentiality. Within MoPAC's framework, this is necessary if there are no control genes available for the normalization of two conditions when computing the differential essentiality scores, or if the control genes should be filtered when computing the essentiality scores.
#'
#' @param id an optional identifier for the analysis. All input will be read from and all output will be stored in a directory with this name. If \code{NULL}, two data frames must be supplied containing the sgRNA-level and gene-level replicate-averaged log-fold-changes. Defaults to \code{NULL}.
#' @param sgRNA an optional data frame containing the replicate-averaged log-fold-changes for each sgRNA. If \code{NULL}, the corresponding file in the identifier directory will be read.
#' @param conditions a vector of condition names to be analyzed. All conditions must be column names in the replicate-averaged log-fold-change datasets. Defaults to all conditions in the dataset.
#' @param pvalue numeric. P-value threshold for depletion/enrichment. Defaults to 0.05.
#' @param fraction numeric. Mininmum fraction of conditions in which a gene must satisfy the p-value threshold. Defaults to 0.8 (i.e. the gene should be significant in at least 80\% of all conditions).
#' @param shinyF logical. An internal parameter used by the graphical user interface. If using command-line leave as \code{NULL}.
#'
#' @details
#' A modified version of the alpha-robust-rank-aggregation algorithm from MAGeCK-RRA is implemented whereby the skewness of the sgRNA distributions is simultaneously assessed for both negative and positive selection.
#' A gene with a significant p value may thus be either significantly depleted or significantly enriched in the screen, and it is therefore necessary to consider the sign of its log-fold-change as well in order to distinguish between negative and positive selection.
#'
#' If \code{id} is not \code{NULL}, the following files are produced in the directory with name equal to the identifier:
#' \itemize{
#'  \item{\code{RRA.txt} files, one for each condition.}
#'  \item{\code{essential.txt}, \code{nonessential.txt}, \code{nondifferential.txt}, \code{enriched.txt} and \code{unenriched.txt}, the filtered lists of essential, nonessential, nondifferential, enriched and unenriched genes.}
#' }
#'
#' @return A list containing the following components:
#' \item{p.value}{List of tables of p values and FDR for all genes of each condition, indicating the ranking of depletion/enrichment.}
#' \item{Essential}{List of filtered essential genes (i.e. those whose 2-tail p value is smaller than the threshold in the minimal fraction of conditions).}
#' \item{Nonessential}{List of filtered non-essential genes (i.e. those whose 2-tail p value is larger than the threshold in the minimal fraction of conditions).}
#' \item{Nondifferential}{List of genes which were not found to be significantly differentially essential.}
#' \item{Enriched}{List of genes which were found to be significantly enriched.}
#' \item{Unenriched}{List of genes which were not found to be significantly enriched in any of the conditions analyzed.}
#'
#' @author Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
#' @keywords RRA essentiality
#'
#' @references Li, W. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol. 15, 554 (2014).
#'
#' @examples
#'
#' ## Example I: Two-tail RRA using an identifier to produce output:
#'
#' ## 1. Prepare read counts information:
#' counts <- read.counts(id="DANG_CCK81", counts=MoPAC::dang_cck81, library=MoPAC::sgRNA_library)
#'
#' ## 2. Manually fill out either of the two annotation files produced by read.counts.
#' ## In this example we will fill it out using R:
#' annotation <- counts$Annotation
#' annotation$Condition <- c("Plasmid0","Plasmid1",rep("DANG",4),rep("CCK81",3))
#' annotation$Replicate <- c("A","A","A","B","C","D","A","B","C")
#' annotation$Control <- c("","",rep("Plasmid0",4),rep("Plasmid1",3))
#' utils::write.csv(annotation,"DANG_CCK81/annotation.csv",quote=FALSE,row.names=FALSE)
#'
#' ## 3. Get pre-processed fold changes:
#' qc <- quality.control(id="DANG_CCK81")
#'
#' ## 4. Finally get the two-tail RRA essentiality analysis:
#' significance <- RRA.2tail(id="DANG_CCK81")
#'
#' ## Example II: Two-tail RRA without using an identifier:
#'
#' ## 1. Prepare read counts information:
#' counts <- read.counts(counts=MoPAC::dang_cck81, library=MoPAC::sgRNA_library)
#'
#' ## 2. Prepare sample annotation as a data frame:
#' annotation <- counts$Annotation
#' annotation$Condition <- c("Plasmid0","Plasmid1",rep("DANG",4),rep("CCK81",3))
#' annotation$Replicate <- c("A","A","A","B","C","D","A","B","C")
#' annotation$Control <- c("","",rep("Plasmid0",4),rep("Plasmid1",3))
#'
#' ## 3. Get pre-processed fold changes:
#' qc <- quality.control(counts=counts$Counts, annotation=annotation)
#'
#' ## 4. Finally get the two-tail RRA essentiality analysis:
#' significance <- RRA.2tail(sgRNA=qc$sgRNA_lfc)
#'
#' ##--------------------------------------------------------------------------------------
#'
#' @export
#'

RRA.2tail <- function(id=NULL,sgRNA=NULL,conditions=NULL,pvalue=0.05,fraction=0.8,shinyF=NULL) {
  set.seed(123)
  print("Importing data...")
  # Import sgRNA-level replicate-averaged log-fold-changes ----------------------------------------------------
  if(is.null(id) & is.null(sgRNA)) { # insufficient information provided
    stop("Please specify a working directory or provide the replicate-averaged sgRNA log-fold-changes.")
  } else if(!is.null(id) & is.null(sgRNA)) { # is a table of read sgRNA_reps available in the working directory?
    if(file.exists(paste0(id,"/sgRNA_lfc.txt"))) {
      sgRNA <- utils::read.table(paste0(id,"/sgRNA_lfc.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Please provide the replicate-averaged sgRNA log-fold-changes generated by the quality control module.")
  } else if(is.data.frame(sgRNA)) {
    sgRNA <- sgRNA
    sgRNA[,sapply(sgRNA,is.factor)] <- sapply(sgRNA[,sapply(sgRNA,is.factor)],as.character)
  } else if(is.character(sgRNA)) { # was a file path provided?
    if(file.exists(sgRNA)) {
      sgRNA <- utils::read.table(sgRNA,header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("sgRNA replicate-averaged log-fold-change file not found in the specified address.")
  } else stop("Please provide a valid table of replicate-averaged sgRNA log-fold-changes.")
  if(!"sgRNA"%in%colnames(sgRNA) | !"Gene"%in%colnames(sgRNA))
    stop("Please run the quality control module first.")
  # Check conditions:
  if(!is.null(conditions)) {
    if(sum(conditions%in%colnames(sgRNA))<length(conditions))
      stop("Selected conditions not found in the column names of replicate-averaged log-fold-changes.")
  } else conditions <- colnames(sgRNA)[-which(colnames(sgRNA)%in%c("sgRNA","Gene","Category"))]
  # Two-tail RRA --------------------------------------------------------------------------------------------------------
  # Remove genes with a single sgRNA:
  NsgRNA <- table(sgRNA$Gene)
  sgRNA <- sgRNA[sgRNA$Gene%in%names(NsgRNA)[NsgRNA>1],]
  # Remove nontargeting sgRNAs:
  if("Category"%in%colnames(sgRNA))
    sgRNA <- sgRNA[sgRNA$Category!="Nontargeting",]
  # Gene-level log-fold-change from average:
  gene <- do.call(rbind,lapply(split(sgRNA[,-which(colnames(sgRNA)%in%c("sgRNA","Gene","Category")),drop=FALSE],sgRNA$Gene),colMeans))
  gene <- data.frame(Gene=rownames(gene),gene,check.names=FALSE,stringsAsFactors=FALSE)
  # Two-tail p values and FDR of essentiality:
  p <- list()
  for(C in conditions) {
    print(paste0("Computing essentiality P-value on ",C,"..."))
    if(!is.null(shinyF)) shinyF(1.0/length(conditions),message='2-tail RRA',detail=paste0("Computing essentiality P-value on ",C))
    p[[C]] <- merge(gene[,c("Gene",C)],getPValues(sgRNA$Gene,sgRNA$sgRNA,sgRNA[,C],0.5))
    p[[C]] <- p[[C]][order(p[[C]]$p,p[[C]]$rho),]
    p[[C]]$FDR <- stats::p.adjust(p[[C]]$p,method="fdr")
    p[[C]]$depleted[p[[C]][,C]<0] <- 1:sum(p[[C]][,C]<0)
    p[[C]]$enriched[p[[C]][,C]>0] <- 1:sum(p[[C]][,C]>0)
    if(!is.null(id)) utils::write.table(p[[C]],paste0(id,"/",C,".RRA.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
  }
  # Control gene and non-differential gene filtering ---------------------------------------------------------------------
  # Find significant genes in each condition:
  # a) Small p value and negative fold change = "low" (essential).
  # b) Large p value = "mid" (nonessential).
  # c) Small p value and positive fold change = "high" (enriched).
  rownames(gene) <- gene$Gene
  low <- mid <- high <- gene[,conditions,drop=FALSE]*0
  for(C in colnames(low)) {
    low[,C] <- as.numeric(rownames(low) %in% p[[C]]$Gene[p[[C]][,C]<0 & p[[C]]$p<pvalue])
    mid[,C] <- as.numeric(rownames(mid) %in% p[[C]]$Gene[p[[C]]$p>pvalue])
    high[,C] <- as.numeric(rownames(high) %in% p[[C]]$Gene[p[[C]][,C]>0 & p[[C]]$p<pvalue])
  }
  # Find genes which are not enriched in any condition:
  high$Any <- apply(high[,,drop=F],1,function(x)sum(x)>0)
  unenriched <- rownames(high)[high$Any==FALSE]
  # Find genes which are significant in more than 80% of the conditions:
  low$Good <- apply(low[,,drop=F],1,function(x)sum(x)>fraction*ncol(low))
  mid$Good <- apply(mid[,,drop=F],1,function(x)sum(x)>fraction*ncol(mid))
  high$Good <- apply(high[,,drop=F],1,function(x)sum(x)>fraction*ncol(high))
  essential <- rownames(low)[low$Good==TRUE]
  nonessential <- rownames(mid)[mid$Good==TRUE]
  enriched <- rownames(high)[high$Good==TRUE]
  nondifferential <- setdiff(union(essential,nonessential),enriched) #generate a new list independent of the original one
  # Filter the original lists of controls according to the significant genes:
  if("Category"%in%colnames(sgRNA)) {
    essential <- intersect(sgRNA$Gene[sgRNA$Category=="Essential"],essential) #filter the original list
    nonessential <- intersect(sgRNA$Gene[sgRNA$Category=="Nonessential"],nonessential) #filter the original list
  }
  # Output -------------------------------------------------------------------------------------------------------------------
  if(!is.null(id)) {
    utils::write.table(data.frame(Essential=essential),paste0(id,"/essential.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(data.frame(Nonessential=nonessential),paste0(id,"/nonessential.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(data.frame(Nondifferential=nondifferential),paste0(id,"/nondifferential.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(data.frame(Enriched=enriched),paste0(id,"/enriched.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(data.frame(Unenriched=unenriched),paste0(id,"/unenriched.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
  }
  print("Done.")
  return(list(p.value=p,Essential=essential,Nonessential=nonessential,Nondifferential=nondifferential,Enriched=enriched,Unenriched=unenriched))
}
