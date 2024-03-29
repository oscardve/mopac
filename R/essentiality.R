#' @title Gene essentiality scoring.
#'
#' @description Generates gene essentiality scores for each condition in the sample annotation file.
#'
#' @param id an identifier for the analysis. All input will be read from and all output will be stored in a directory with this name. If \code{NULL}, two data frames must be supplied containing the sgRNA-level log-fold-changes per replicate and annotation. Defaults to \code{NULL}.
#' @param sgRNA_reps an optional data frame containing the sgRNA log-fold-changes for each replicate. If \code{NULL}, the log-fold-change file generated previously will be read.
#' @param annotation an optional data frame containing the sample annotation information. If \code{NULL}, the annotation file generated previously will be read.
#' @param uniform.weights logical. Should uniform weight values be used for sgRNA weighted averaging? Defaults to \code{FALSE}.
#' @param empirical.weights logical. Should empirical weight values be used for sgRNA weighted averaging? Defaults to \code{FALSE}.
#' @param read.filtered logical. Should the filtered list of control genes generated by the 2-tail RRA algorithm be used to optimize the regression parameters? Requires an identifier to be set, otherwise use \code{filtered.essential} and \code{filtered.nonessential} instead. Defaults to \code{FALSE}.
#' @param filtered.essential an optional vector of names for filtered essential genes. Takes precedence over \code{read.filtered}. Defaults to \code{NULL}.
#' @param filtered.nonessential an optional vector of names for filtered non-essential genes. Takes precedence over \code{read.filtered}. Defaults to \code{NULL}.
#' @param filtered.reps an optional vector of names for filtered sample names. Defaults to \code{NULL}.
#' @param shinyF logical. An internal parameter used by the graphical user interface. If using command-line leave as \code{NULL}.
#'
#' @details
#' The gene-scoring parameters are first optimized for those genes whose size is equal to the most frequent gene size.
#' Genes with a larger size are scored by bootstrap-aggregating them 10 times so that the size of each bootstrap is equal to the most frequent size, and the same parameters are used.
#' Genes with smaller size are scored by generating 10 sets of scoring parameters from bootstrapping the control genes for each gene size, and afterwards bootstrap-aggregating the scores from each set.
#' The posibility of employing filtered sets of control genes requires running the 2tail RRA module beforehand. \cr
#'
#' If \code{id} is not \code{NULL}, the following files are produced in the directory with name equal to the identifier:
#' \itemize{
#'  \item{\code{scored_reps.txt} and \code{scored.txt}, containing the gene-level scores for each condition before and after averaging the replicates.}
#'  \item{\code{scored_weights.txt}, \code{scored_scaling.txt} and \code{scored_shifting.txt}, the optimized final regression parameters employing the selected control genes whose size is equal to the most frequent gene size in the dataset. Only produced if \code{empirical.weights} is set to \code{TRUE}}
#' }
#'
#' @return A list containing the following components:
#' \item{scored_reps}{gene scores for each replicate.}
#' \item{scored}{replicate-averaged gene scores.}
#' \item{weights, scaling and shifting}{the optimized final regression parameters employing the selected control genes whose size is equal to the most frequent gene size in the dataset. If not applicable, returns \code{NULL}.}
#'
#' @author Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
#' @keywords essentiality
#'
#' @examples
#'
#' ## Example I. Generating output files through an analysis identifier:
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
#' qc <- quality.control(id="DANG_CCK81", report=FALSE)
#'
#' ## 4. Optionally do two-tail RRA essentiality analysis for control gene filtering:
#' significance <- RRA.2tail(id="DANG_CCK81")
#'
#' ## 5. Gene essentiality analysis:
#' scores <- analyze.essentiality(id="DANG_CCK81",empirical.weights=TRUE)
#'
#' ## Example II. Returning a data frame without generating output:
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
#' ## 4. Optionally do two-tail RRA essentiality analysis:
#' significance <- RRA.2tail(sgRNA=qc$sgRNA_lfc)
#'
#' ## 5. Gene essentiality analysis:
#' essentiality <- analyze.essentiality(sgRNA_reps=qc$sgRNA_lfc_reps,
#'                                      empirical.weights=TRUE, annotation=annotation)
#'
#' ##--------------------------------------------------------------------------------------
#'
#' @export
#'

analyze.essentiality <- function(id=NULL,sgRNA_reps=NULL,annotation=NULL,uniform.weights=FALSE,empirical.weights=FALSE,read.filtered=FALSE,filtered.essential=NULL,filtered.nonessential=NULL,filtered.reps=NULL,shinyF=NULL) {
  # id=NULL;sgRNA_reps=NULL;annotation=NULL;empirical.weights=FALSE;read.filtered=FALSE;filtered.essential=NULL;filtered.nonessential=NULL;filtered.reps=NULL;shinyF=NULL
  # sgRNA_reps <- read.table("~/paper_feb20/cbp30_mopac_selftrained/sgRNA_lfc_reps.txt",
  #                          header=TRUE,sep="",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,row.names=NULL)
  # annotation <- read.csv("~/paper_feb20/cbp30_mopac_selftrained/annotation.csv")
  # filtered.essential <- read.table("~/paper_feb20/cbp30_mopac_selftrained/essential.txt",header=T,stringsAsFactors=F)[,1]
  # filtered.nonessential <- read.table("~/paper_feb20/cbp30_mopac_selftrained/nonessential.txt",header=T,stringsAsFactors=F)[,1]

  set.seed(123)
  # Read sgRNA log-fold-changes  -------------------------------------------------------------------------------
  print("Importing data..."); if(!is.null(shinyF)) shinyF(0.1,"Essentiality analysis","Importing data...")
  if(is.null(id) & is.null(sgRNA_reps)) { # insufficient information provided
    stop("Please specify a working directory or provide the sgRNA log-fold-changes and sample annotation.")
  } else if(!is.null(id) & is.null(sgRNA_reps)) { # is a table of read sgRNA_reps available in the working directory?
    if(file.exists(paste0(id,"/sgRNA_lfc_reps.txt"))) {
      sgRNA <- utils::read.table(paste0(id,"/sgRNA_lfc_reps.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
      # } else if(file.exists(paste0(id,"/fastq_sgRNA_reps.txt"))) {
      #   sgRNA <- utils::read.table(paste0(id,"/fastq_sgRNA_reps.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Please provide the sgRNA log-fold-changes and sample annotation generated by the quality control module.")
  } else if(is.data.frame(sgRNA_reps)) { # was a data frame provided as a table of read sgRNA_reps?
    sgRNA <- sgRNA_reps
    sgRNA[,sapply(sgRNA,is.factor)] <- sapply(sgRNA[,sapply(sgRNA,is.factor)],as.character)
    # colnames(sgRNA) <- stringr::str_replace_all(colnames(sgRNA),"[^[:alnum:]]","") #remove non-alphanumeric characters in the sample names
  } else if(is.character(sgRNA_reps)) { # was a file path provided as a table of read sgRNA_reps?
    if(file.exists(sgRNA_reps)) {
      sgRNA <- utils::read.table(sgRNA_reps,header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("sgRNA log-fold-change file not found in the specified address.")
    # colnames(sgRNA) <- stringr::str_replace_all(colnames(sgRNA),"[^[:alnum:]]","") #remove non-alphanumeric characters in the sample names
  } else stop("Please provide a valid table of sgRNA log-fold-changes.")
  if(!"sgRNA"%in%colnames(sgRNA) | !"Gene"%in%colnames(sgRNA))
    stop("Please run the quality control module first.")
  if(uniform.weights==FALSE & empirical.weights==FALSE & !"Category"%in%colnames(sgRNA))
    stop("Please set empirical.weights to TRUE in the absence of control genes.")
  # Read annotation -------------------------------------------------------------------------
  if(is.null(id) & is.null(annotation)) { # insufficient information provided
    stop("Please specify a working directory or provide a data frame of annotation.")
  } else if(!is.null(id) & is.null(annotation)) { # is a table of annotation available in the working directory?
    if(file.exists(paste0(id,"/annotation.csv"))) {
      annotation0 <- utils::read.csv(paste0(id,"/annotation.csv"),header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
      if(sum(is.na(annotation0[,-which(colnames(annotation0)%in%"Control")]))>0)
        annotation0 <- openxlsx::read.xlsx(paste0(id,"/annotation.xlsx"),check.names=FALSE,na.strings="")
      if(sum(is.na(annotation0[,-which(colnames(annotation0)%in%"Control")]))>0)
        stop("Please fill out the sample file annotation before proceeding.")
    } else stop("Please provide a table of annotation or run the input module first.")
  } else if(is.data.frame(annotation)) { # was a data frame provided as a table of read counts?
    annotation0 <- annotation
    annotation0[,sapply(annotation0,is.factor)] <- sapply(annotation0[,sapply(annotation0,is.factor)],as.character)
  }  else if(is.character(annotation)) {
    if(file.exists(annotation)) {
      if(tools::file_ext(annotation)=="xlsx") {
        annotation0 <- openxlsx::read.xlsx(annotation,check.names=FALSE,na.strings="")
      } else if(tools::file_ext(annotation)=="csv") {
        annotation0 <- utils::read.csv(annotation,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
      } else annotation0 <- utils::read.table(annotation,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
    } else stop("Annotation file not found in the specified address.")
    if(sum(is.na(annotation0[,-which(colnames(annotation0)%in%"Control")]))>0)
      stop("Please fill out the sample file annotation before proceeding.")
  } else stop("The annotation object provided could not be read.")
  if(sum(!c("Sample","Condition","Replicate","Control")%in%colnames(annotation0))>0)
    stop("Please make sure the annotation contains the columns: Sample, Condition, Replicate and Control")
  # Filter controls only if the RRA module requests it (the files such as "essential.txt" were already filtered if that is the case):
  if("Category"%in%colnames(sgRNA)) {
    if(!is.null(filtered.essential)) {
      sgRNA$Category[sgRNA$Category=="Essential" & !sgRNA$Gene%in%filtered.essential] <- "Other"
    } else if(read.filtered==TRUE) {
      if(!file.exists(paste0(id,"/essential.txt"))) stop(paste0("File ",id,"/essential.txt not found."))
      filtered.essential <- utils::read.table(paste0(id,"/essential.txt"),header=TRUE,sep="",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,row.names=NULL)$Essential
      sgRNA$Category[sgRNA$Category=="Essential" & !sgRNA$Gene%in%filtered.essential] <- "Other"
    }
    if(!is.null(filtered.nonessential)) {
      sgRNA$Category[sgRNA$Category=="Nonessential" & !sgRNA$Gene%in%filtered.nonessential] <- "Other"
    } else if(read.filtered==TRUE) {
      if(!file.exists(paste0(id,"/nonessential.txt"))) stop(paste0("File ",id,"/nonessential.txt not found."))
      filtered.nonessential <- utils::read.table(paste0(id,"/nonessential.txt"),header=TRUE,sep="",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,row.names=NULL)$Nonessential
      sgRNA$Category[sgRNA$Category=="Nonessential" & !sgRNA$Gene%in%filtered.nonessential] <- "Other"
    }
  }
  # Filter replicates:
  if(!is.null(filtered.reps)) {
    sgRNA <- sgRNA[,which(colnames(sgRNA)%in%c("sgRNA","Gene","Category",filtered.reps))]
    annotation0 <- annotation0[annotation0$Sample%in%filtered.reps,]
    if(nrow(annotation0)==0) stop("Replicates chosen not found.")
  }
  # Get scores:
  if(uniform.weights==TRUE) {
    scored_reps <- do.call(rbind,lapply(split(sgRNA[,-which(colnames(sgRNA)%in%c("sgRNA","Gene","Category")),drop=FALSE],sgRNA$Gene),colMeans))
    scored_reps <- data.frame(Gene=rownames(scored_reps),scored_reps,check.names=FALSE,stringsAsFactors=FALSE)
    if("Category"%in%colnames(sgRNA))
      scored_reps <- merge(unique(sgRNA[,which(colnames(sgRNA)%in%c("Gene","Category")),drop=FALSE]),scored_reps)
  } else scoring <- get.scores(sgRNA,uniform.weights,empirical.weights,shinyF=shinyF)
  scored_reps <- scoring$Scores
  if(uniform.weights==FALSE & empirical.weights==FALSE) {
    weights <- data.frame(Rank=1:length(scoring$W),Weight=scoring$W)
    scaling <- scoring$A
    shifting <- scoring$B
  } else weights <- scaling <- shifting <- NULL
  # Replicate average, standard deviation and p value:
  annotation1 <- annotation0[annotation0$Sample%in%colnames(scored_reps),]
  if(nrow(annotation1)==0)
    stop("Mismatch found between the annotation samples and the sgRNA column names.")
  groups <- split(as.character(annotation1$Sample),annotation1[c("Condition")],drop=T)
  scored <- cbind(scored_reps[,which(colnames(scored_reps)%in%c("Gene","Category")),drop=F],
                  do.call(cbind,lapply(groups,function(y)rowMeans(scored_reps[,y,drop=F]))))
  # Save output:
  if(!is.null(id)) {
    utils::write.table(scored_reps,paste0(id,"/scored_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(scored,paste0(id,"/scored.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    if(uniform.weights==FALSE & empirical.weights==FALSE) {
      utils::write.table(weights,paste0(id,"/scored_weights.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
      utils::write.table(scaling,paste0(id,"/scored_scaling.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
      utils::write.table(shifting,paste0(id,"/scored_shifting.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    }
  }
  print("Done.")
  return(list(scored_reps=scored_reps,scored=scored,weights=weights,scaling=scaling,shifting=shifting))
}
