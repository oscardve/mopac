#' @title Pre-processing and quality assessment of CRISPR screen.
#'
#' @description Pre-processes and assesses the quality of a CRISPR screen by computing log-fold-changes with optimized pseudo-count, and generates a pdf report.
#'
#' @param id an optional identifier for the analysis. All input and output will be stored in a directory with this name. No output files will be produced if not supplied.
#' @param counts an optional data frame or file path containing the read count information, which may be in either integer or normalized form. If \code{NULL}, the counts file in the id directory will be read. Defaults to \code{NULL}.
#' @param annotation an optional data frame or file path containing the sample annotation information. If \code{NULL}, the annotation file in the id directory will be read. Defaults to \code{NULL}.
#' @param pseudocount numeric. An optional pseudo-count for log2 transformation. If \code{NULL}, default optimization is performed over the values 0.01, 0.05, 0.1, 0.15, 0.2, 0.25 and 0.3. Defaults to 0.15.
#' @param report logical. Should a pdf report of quality assessment be generated? Requires the id to not be \code{NULL}, and a pandoc installation. Defaults to \code{FALSE}.
#' @param shinyF logical. An internal parameter used by the graphical user interface. If using command-line mode leave as \code{NULL}.
#'
#' @details
#' Log-fold-changes at the sgRNA level are computed before and after averaging the replicates. If an identifier is supplied, the following files are produced in the directory with name equal to the identifier:
#' \itemize{
#'  \item{\code{sgRNA_log2.txt}, the log2-transformed read counts normalized by sequencing depth.}
#'  \item{\code{sgRNA_lfc_reps.txt} and \code{sgRNA_lfc.txt}, the log-fold-changes before and after averaging the replicates.}
#'  \item{\code{qc_variance.txt}, \code{qc_weights.txt}, \code{qc_scaling.txt} and \code{qc_shifting.txt}, the optimized preliminary regression parameters employing the unfiltered list of control genes whose size is equal to the most frequent gene size in the dataset.}
#'  \item{\code{qc.pdf}, a report for the visualization of the quality assessment.}
#' }
#'
#' @return A list containing the following components:
#' \item{Counts}{raw integer counts.}
#' \item{Pseudo}{optimal pseudocount and pseudocount versus variance function.}
#' \item{Annotation}{a data frame with sample annotation information, which is merged for fastq files belonging to the same sample.}
#' \item{Annotation1}{a data frame with sample annotation information, not including the control (i.e. day-0) conditions.}
#' \item{sgRNA_log_reps}{the log-counts obtained through the optimal or user-defined pseudocount.}
#' \item{sgRNA_lfc_reps}{log-fold-changes for each sgRNA of each replicate.}
#' \item{sgRNA_lfc}{log-fold-changes for the repliate-averaged sgRNAs.}
#'
#' @author Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
#' @keywords quality control
#'
#' @examples
#'
#' ## Example I. Pre-processing read counts without producing output:
#'
#' ## 1. Prepare counts as a data frame:
#' counts <- read.counts(counts = MoPAC::dang_cck81, library=MoPAC::sgRNA_library)
#'
#' ## 2. Prepare sample annotation as a data frame:
#' annotation <- counts$Annotation
#' annotation$Condition <- c("Plasmid0","Plasmid1",rep("DANG",4),rep("CCK81",3))
#' annotation$Replicate <- c("A","A","A","B","C","D","A","B","C")
#' annotation$Control <- c("","",rep("Plasmid0",4),rep("Plasmid1",3))
#'
#' ## 3. Get sgRNA-level and gene-level fold changes:
#' qc <- quality.control(counts=counts$Counts, annotation=annotation)
#'
#'
#' ##--------------------------------------------------------------------------------------
#'
#' @export
#'

quality.control <- function(id=NULL,counts=NULL,annotation=NULL,pseudocount=0.15,report=FALSE,shinyF=NULL) {
  #######################################################################################################################################
  ########################### INPUT #####################################################################################################
  #######################################################################################################################################
  if(!is.null(id)) { if(!dir.exists(id)) stop("The working directory selected was not found") }
  # Read counts -------------------------------------------------------------------------
  print("Importing data..."); if(!is.null(shinyF)) shinyF(0.1,"Quality Control","Importing data...")
  if(is.null(id) & is.null(counts)) { # insufficient information provided
    stop("Please specify a working directory or provide a table of read counts.")
  } else if(!is.null(id) & is.null(counts)) { # is a table of read counts available in the working directory?
    if(file.exists(paste0(id,"/counts.txt"))) {
      counts0 <- utils::read.table(paste0(id,"/counts.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Please provide a table of read counts or run the FASTQ alignment module first.")
  } else if(is.data.frame(counts)) { # was a data frame provided as a table of read counts?
    counts0 <- counts
    counts0[,sapply(counts0,is.factor)] <- sapply(counts0[,sapply(counts0,is.factor)],as.character)
  } else if(is.character(counts)) { # was a file path provided as a table of read counts?
    if(file.exists(counts)) {
      if(tools::file_ext(counts)=="xlsx") {
        counts0 <- openxlsx::read.xlsx(counts)
      } else if(tools::file_ext(counts)=="csv") {
        counts0 <- utils::read.csv(counts,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"")
      } else counts0 <- utils::read.table(counts,header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Counts file not found in the specified address.")
  } else stop("Please provide a valid table of read counts.")
  colnames(counts0)[grep("sgRNA|sgrna|SGRNA",colnames(counts0))] <- "sgRNA"
  colnames(counts0)[grep("Gene|gene|GENE",colnames(counts0))] <- "Gene"
  colnames(counts0)[grep("Category|category|CATEGORY",colnames(counts0))] <- "Category"
  if(!"sgRNA"%in%colnames(counts0) | !"Gene"%in%colnames(counts0))
    stop("Please run the input module first.")
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
  # Find groups with common reference sample:
  annotation0$Control[is.na(annotation0$Control)] <- annotation0$Condition[is.na(annotation0$Control)]
  annotation0$Control[annotation0$Control==""] <- annotation0$Condition[annotation0$Control==""]
  if(sum(annotation0$Control%in%annotation0$Condition)!=nrow(annotation0))
    stop("Please make sure that the control conditions exist in the 'Condition' column.")
  initial <- annotation0$Sample[annotation0$Condition%in%unique(annotation0$Control)]
  annotation1 <- annotation0[-which(annotation0$Sample%in%initial),] #remove day 0
  if(nrow(annotation1)==0) {
    print("No future timepoints found. Returning table of read counts...")
    if(!is.null(shinyF))
      shinyF(0.1,"Quality Control","No future timepoints found. Returning table of read counts...")
    return(counts)
  }
  #######################################################################################################################################
  ########################### PROCESSING ################################################################################################
  #######################################################################################################################################
  print("Pre-processing data...")
  if(!is.null(shinyF))
    shinyF(0.1,"Quality Control","Pre-processing data...")
  # Optimize pseudocount and normalize data if the data is not normalized --------------------------------------
  groups <- split(annotation0[,"Sample",drop=F],annotation0$Control)
  isNormed <- sum(counts0[,which(sapply(counts0,is.numeric))[1]]%%1!=0)>0 #check if data is normalized or not
  if(!isNormed) {
    if(is.null(pseudocount)) {
      if(!"Category"%in%colnames(counts0))
        stop("Cannot optimize the pseudocount without gene category information.")
      print("Optimizing parameters")
      if(!is.null(shinyF))
        shinyF(0.1,"Quality Control","Optimizing parameters...")
      pseudorange <- c(0.01,seq(0.05,0.3,0.05)) # if the pseudocount is null, the pseudocount will be optimized within this range
      pseudo <- get.pseudocount(counts0,initial,groups,pseudorange,shinyF)
      pseudocount <- pseudo$Optimal$Pseudocount
    } else pseudo <- NULL
    print(paste0("Using a pseudocount of ",pseudocount))
    sgRNA_log_reps <- cbind(counts0[,which(colnames(counts0)%in%c("sgRNA","Gene","Category"))],
                    apply(counts0[,-which(colnames(counts0)%in%c("sgRNA","Gene","Category"))],2,function(x)log2(x/stats::median(x)+pseudocount)))
  } else {
    sgRNA_log_reps <- counts0 #if the data is already normalized, make no change
    pseudo <- NULL
  }
  sgRNA_lfc_reps <- cbind(sgRNA_log_reps[,which(colnames(sgRNA_log_reps)%in%c("sgRNA","Gene","Category"))],do.call(cbind,unname(lapply(groups,function(y)
    sgRNA_log_reps[,setdiff(y$Sample,initial),drop=F]-sgRNA_log_reps[,intersect(y$Sample,initial),drop=T]))))
  # Gene level -------------------------------------------------------------------------------------------------
  # Remove genes with a single sgRNA:
  NsgRNA <- table(sgRNA_log_reps$Gene)
  gene_log_reps <- sgRNA_log_reps[sgRNA_log_reps$Gene%in%names(NsgRNA)[NsgRNA>1],]
  # Remove nontargeting sgRNAs:
  gene_log_reps <- gene_log_reps[gene_log_reps$Category!="Nontargeting",]
  # Average sgRNA log2 reads per gene:
  gene_log_reps <- do.call(rbind,lapply(split(sgRNA_log_reps[,-which(colnames(sgRNA_log_reps)%in%c("sgRNA","Gene","Category")),drop=FALSE],sgRNA_log_reps$Gene),colMeans))
  gene_log_reps <- data.frame(Gene=rownames(gene_log_reps),gene_log_reps,check.names=FALSE,stringsAsFactors=FALSE)
  if("Category"%in%colnames(sgRNA_log_reps))
    gene_log_reps <- merge(unique(sgRNA_log_reps[,which(colnames(sgRNA_log_reps)%in%c("Gene","Category")),drop=FALSE]),gene_log_reps)
  # Gene log-fold-change:
  gene_lfc_reps <- cbind(gene_log_reps[,which(colnames(gene_log_reps)%in%c("Gene","Category"))],do.call(cbind,unname(lapply(groups,function(y)
    gene_log_reps[,setdiff(y$Sample,initial),drop=F]-gene_log_reps[,intersect(y$Sample,initial),drop=T]))))
  # Replicate average -------------------------------------------------------------------------------------------
  # Log reads:
  groups <- split(as.character(annotation0$Sample),annotation0[,c("Condition"),drop=F],drop=T)
  sgRNA_log <- do.call(cbind,lapply(groups,function(y)rowMeans(sgRNA_log_reps[,y,drop=F])))
  sgRNA_log <- cbind(sgRNA_log_reps[,which(colnames(sgRNA_log_reps)%in%c("sgRNA","Gene","Category"))],sgRNA_log)
  gene_log <- do.call(cbind,lapply(groups,function(y)rowMeans(gene_log_reps[,y,drop=F])))
  gene_log <- cbind(gene_log_reps[,which(colnames(gene_log_reps)%in%c("Gene","Category")),drop=FALSE],gene_log)
  # Fold changes:
  groups <- split(as.character(annotation1$Sample),annotation1[,c("Condition"),drop=F],drop=T)
  sgRNA_lfc <- do.call(cbind,lapply(groups,function(y)rowMeans(sgRNA_lfc_reps[,y,drop=F])))
  sgRNA_lfc <- cbind(sgRNA_lfc_reps[,which(colnames(sgRNA_lfc_reps)%in%c("sgRNA","Gene","Category"))],sgRNA_lfc)
  gene_lfc <- do.call(cbind,lapply(groups,function(y)rowMeans(gene_lfc_reps[,y,drop=F])))
  gene_lfc <- cbind(gene_lfc_reps[,which(colnames(gene_lfc_reps)%in%c("Gene","Category"))],gene_lfc)
  #######################################################################################################################################
  ########################### OUTPUT ####################################################################################################
  #######################################################################################################################################
  rowSort <- order(sgRNA_log_reps$Gene,sgRNA_log_reps$sgRNA) #careful: some sgRNAs are repeated for multiple genes!
  infoCols <- colnames(sgRNA_log_reps)[colnames(sgRNA_log_reps)%in%c("sgRNA","Gene","Category")]
  sgRNA_log_reps <- sgRNA_log_reps[rowSort,c(infoCols,annotation0$Sample[order(annotation0$Condition)])]
  sgRNA_lfc_reps <- sgRNA_lfc_reps[rowSort,c(infoCols,annotation1$Sample[order(annotation1$Condition)])]
  sgRNA_log <- sgRNA_log[rowSort,c(infoCols,gtools::mixedsort(colnames(sgRNA_log)[-which(colnames(sgRNA_log)%in%c("sgRNA","Gene","Category"))]))]
  sgRNA_lfc <- sgRNA_lfc[rowSort,c(infoCols,gtools::mixedsort(colnames(sgRNA_lfc)[-which(colnames(sgRNA_lfc)%in%c("sgRNA","Gene","Category"))]))]
  infoCols <- colnames(gene_lfc_reps)[colnames(gene_lfc_reps)%in%c("Gene","Category")]
  gene_log_reps <- gene_log_reps[,c(infoCols,gtools::mixedsort(colnames(gene_log_reps)[-which(colnames(gene_log_reps)%in%c("Gene","Category"))]))]
  gene_lfc_reps <- gene_lfc_reps[,c(infoCols,gtools::mixedsort(colnames(gene_lfc_reps)[-which(colnames(gene_lfc_reps)%in%c("Gene","Category"))]))]
  gene_log <- gene_log[,c(infoCols,gtools::mixedsort(colnames(gene_log)[-which(colnames(gene_log)%in%c("Gene","Category"))]))]
  gene_lfc <- gene_lfc[,c(infoCols,gtools::mixedsort(colnames(gene_lfc)[-which(colnames(gene_lfc)%in%c("Gene","Category"))]))]
  if(!is.null(id)) {
    print("Generating output...")
    if(!is.null(shinyF))
      shinyF(0.1,"Quality Control","Generating output...")
    utils::write.csv(annotation0,paste0(id,"/annotation.csv"),quote=FALSE,row.names=FALSE)
    utils::write.table(counts0,paste0(id,"/counts.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(sgRNA_log_reps,paste0(id,"/sgRNA_log_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(sgRNA_lfc_reps,paste0(id,"/sgRNA_lfc_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(sgRNA_log,paste0(id,"/sgRNA_log.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(sgRNA_lfc,paste0(id,"/sgRNA_lfc.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(gene_log_reps,paste0(id,"/gene_log_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(gene_lfc_reps,paste0(id,"/gene_lfc_reps.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(gene_log,paste0(id,"/gene_log.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    utils::write.table(gene_lfc,paste0(id,"/gene_lfc.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    if(!isNormed & !is.null(pseudo)) {
      utils::write.table(pseudo$Variance,paste0(id,"/qc_variance.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
      utils::write.table(pseudo$Optimal,paste0(id,"/qc_pseudocount.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
    }
    # Generate PDF file of quality control --------------------------------------------------------------------
    if(report==TRUE) {
      if(Sys.which("pandoc")=="")
        stop("MoPAC requires pandoc to be installed in order to generate reports. Please (1) run MoPAC within Rstudio, which already comes with pandoc installed or (2) check your pandoc installation.")
      print("Generating quality control report...")
      if(!is.null(shinyF))
        shinyF(0.1,"Quality Control","Generating quality control report...")
      if(file.exists(paste0(id,"/fastq.txt"))){
        rmarkdown::render(system.file("rmd/qcfastq.Rmd",package="MoPAC"),output_file="qc.pdf",output_dir=id,#intermediates_dir=id,knit_root_dir=id,
                          params = list( #paste0(getwd(),"/",id)
                            annotation0 = annotation0,
                            annotation = annotation,
                            annotation1 = annotation1,
                            counts = counts0,
                            pseudo = pseudo,
                            logged = sgRNA_log_reps,
                            sgRNA_reps = sgRNA_lfc_reps,
                            sgRNA = sgRNA
                            # gene = gene
                          ))
      } else if(isNormed) {
      } else {
        rmarkdown::render(system.file("rmd/qctable.Rmd",package="MoPAC"),output_file="qc.pdf",output_dir=id,#intermediates_dir=id,knit_root_dir=id,
                          params = list( #paste0(getwd(),"/",id)
                            annotation = annotation,
                            annotation1 = annotation1,
                            pseudo = pseudo,
                            counts = counts0,
                            logged = sgRNA_log_reps,
                            sgRNA_reps = sgRNA_lfc_reps,
                            sgRNA = sgRNA
                            # gene = gene
                          ))
      }
    }
  }
  print("Done.")
  return(list(Counts=counts0, Pseudo=pseudo, Annotation=annotation0, Annotation1=annotation1,
              sgRNA_log_reps=sgRNA_log_reps, sgRNA_lfc_reps=sgRNA_lfc_reps, sgRNA_log=sgRNA_log, sgRNA_lfc=sgRNA_lfc,
              gene_log_reps=gene_log_reps, gene_log=gene_log, gene_lfc_reps=gene_lfc_reps, gene_lfc=gene_lfc))
}
