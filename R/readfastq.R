#' @title Map FASTQ files to an sgRNA library.
#'
#' @description Maps the annotation in all FASTQ files contained in a specified directory to a custom sgRNA library and generates a sample annotation file which needs to be filled out for further analysis.
#'
#' @param id working directory. All future output will be stored in a directory with this name. Defaults to "out".
#' @param library sgRNA library. It may be a data frame or a file path. It must contain a column called \code{"sgRNA"} and another called \code{"Gene"}. It may also contain a \code{"Category"} column. Only necessary if the table of counts does not include a column called "Gene".
#' @param fastq path of the directory where the FASTQ files are located. All files should have .fastq extension.
#' @param spacer_start integer value (starting from 0) specifying the starting position of the sgRNA sequences. Defaults to \code{0}.
#' @param spacer_length integer value specifying the length of the sgRNA sequences. Defaults to \code{19}.
#' @param reverse logical. Should sgRNA sequences be reversed before mapping? Defaults to \code{TRUE}.
#' @param complement logical. Should sgRNA sequences be complemented before mapping? Defaults to \code{TRUE}.
#' @param essential an optional vector of gene names to add genes categorized as essential.
#' @param nonessential an optional vector of gene names to add genes categorized as non-essential.
#' @param nontargeting an optional vector of "gene" names to add the default set of "genes" categorized as non-targeting.
#' @param shinyF logical. An internal parameter used by the graphical user interface. If using command-line mode leave as \code{NULL}.
#'
#' @details
#' The following files are produced in the working directory:
#' \itemize{
#' \item{\code{fastq_library.txt}, containing two columns: \code{"sgRNA"} and \code{"Gene"} for each sgRNA mapped to the library file, and possibly a \code{"Category"} column as well if supplied.}
#' \item{\code{fastq_annotation.csv} and \code{fastq_annotation.xlsx}, two annotation files, either of which may be chosen to fill out the following FASTQ file information: cell line, treatment, timepoint, replicate and batch. The annotation files also specify the proportion and percentage of annotation mapped, number of zero counts and gini index per FASTQ file.}
#' \item{\code{fastq_counts.txt}, containing the read counts per sgRNA for each FASTQ file.}
#' \item{\code{fastq_reads.txt}, containing the number of reads per FASTQ file.}
#' }
#'
#' @return A list containing the following components:
#' \item{Library}{The sgRNA library.}
#' \item{Annotation}{The table of fastq file annotation.}
#' \item{Counts}{The table of read counts at the FASTQ file level.}
#'
#' @author Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
#' @keywords read fastq
#'
#' @examples
#' ## Examples and sample FASTQ files available at:
#'
#' ## reads <- read.fastq(id="out",fastq="dir",library="lib.txt",
#' ##                     spacer_start=0,spacer_length=19,reverse=TRUE,complement=TRUE)
#'
#' ##--------------------------------------------------------------------------------------
#'
#' @export
#'

read.fastq <- function(id="out",library=NULL,fastq=NULL,spacer_start=0,spacer_length=19,reverse=TRUE,complement=TRUE,essential=NULL,nonessential=NULL,nontargeting=NULL,shinyF=NULL) {
  # Read library --------------------------------------------------------------------------------------------
  if(is.null(id) & is.null(library)) {
    stop("Please supply an sgRNA library or a working directory containing one.")
  } else if(!is.null(id) & is.null(library)) {
    if(file.exists(paste0(id,"/fastq_library.txt")))
      library0 <- utils::read.table(library,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
  } else if(is.data.frame(library)) {
    library0 <- library
    library0[,sapply(library0,is.factor)] <- sapply(library0[,sapply(library0,is.factor)],as.character)
  } else if(is.character(library)) {
    if(file.exists(library)) {
      if(tools::file_ext(library)=="xlsx") {
        library0 <- openxlsx::read.xlsx(library)
      } else if(tools::file_ext(library)=="csv") {
        library0 <- utils::read.csv(library,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"")
      } else library0 <- utils::read.table(library,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("The library file provided was not found in the path provided.")
  } else stop("Please provide a valid sgRNA library.")
  # Organize library:
  colnames(library0)[grep("sgRNA|sgrna|SGRNA",colnames(library0))] <- "sgRNA"
  colnames(library0)[grep("Gene|gene|GENE",colnames(library0))] <- "Gene"
  colnames(library0)[grep("Category|category|CATEGORY",colnames(library0))] <- "Category"
  if(!"sgRNA"%in%colnames(library0) | !"Gene"%in%colnames(library0))
    stop("Please make sure the library file contains a column called 'sgRNA' and another called 'Gene'.")
  if("Category"%in%colnames(library0)) {
    library0 <- library0[,c("sgRNA","Gene","Category")]
  } else {
    library0 <- library0[,c("sgRNA","Gene")]
    library0$Category <- "Other"
  }
  library0$Category[library0$Gene%in%MoPAC::essentials] <- "Essential"
  library0$Category[library0$Gene%in%MoPAC::nonessentials] <- "Nonessential"
  library0$Category[library0$Gene%in%essential] <- "Essential"
  library0$Category[library0$Gene%in%nonessential] <- "Nonessential"
  library0$Category[library0$Gene%in%nontargeting] <- "Nontargeting"
  # Create output directory and save library file inside:
  dir.create(id,showWarnings=FALSE)
  utils::write.table(library0,paste0(id,"/fastq_library.txt"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE)
  # Count sgRNAs in the FASTQ files  --------------------------------------------------------------------------------
  files <- gtools::mixedsort(list.files(path=fastq,pattern="\\.fastq$",full.names=TRUE))
  if(identical(files,character(0)))
    stop("No FASTQ files were found in this directory.")
  if(is.null(shinyF)) {
    readFASTQC(library=paste0(id,"/fastq_library.txt"),id=id,spacer_start=spacer_start,spacer_length=spacer_length,rev=as.numeric(reverse),comp=as.numeric(complement),fastqIN=files)
  } else readFASTQC_shiny(library=paste0(id,"/fastq_library.txt"),id=id,spacer_start=spacer_start,spacer_length=spacer_length,rev=as.numeric(reverse),comp=as.numeric(complement),fastqIN=files,shinyF=shinyF)
  counts0 <- utils::read.table(paste0(id,"/fastq_counts.txt"),header=TRUE,sep="",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
  annotation <- utils::read.table(paste0(id,"/fastq_reads.txt"),header=TRUE,sep="",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
  if(nrow(counts0)==0)
    stop("The sequences in the sgRNA library were not found in the FASTQ files. Please check the library.")
  # colnames(counts0) <- stringr::str_replace_all(colnames(counts0),"[^[:alnum:]]","") #remove non-alphanumeric characters in the sample names
  # annotation$File <- stringr::str_replace_all(annotation$File,"[^[:alnum:]]","") #remove non-alphanumeric characters in the sample names
  counts1 <- merge(library0,counts0,all.x=TRUE)
  counts1[is.na(counts1)] <- 0
  library1 <- counts1[,which(colnames(counts1)%in%c("sgRNA","Gene","Category"))]
  # Compute quality control information of FASTQ files  --------------------------------------------------------------
  mapped <- data.frame(File=colnames(counts1)[-which(colnames(counts1)%in%c("sgRNA","Gene","Category"))],Mapped=colSums(counts1[,-which(colnames(counts1)%in%c("sgRNA","Gene","Category"))]),stringsAsFactors=FALSE)
  annotation1 <- merge(annotation,mapped)
  annotation1[,"Mapped%"] <- annotation1[,"Mapped"]/annotation1$Total_reads*100
  # annotation1$sgRNAs <- nrow(counts1)
  # zeros <- data.frame(File=colnames(counts1)[-which(colnames(counts1)%in%c("sgRNA","Gene","Category"))],Zeros=apply(counts1[,-which(colnames(counts1)%in%c("sgRNA","Gene","Category"))],2,function(y)sum(y==0)),stringsAsFactors=FALSE)
  # gini <- data.frame(File=colnames(counts1)[-which(colnames(counts1)%in%c("sgRNA","Gene","Category"))],Gini=apply(counts1[,-which(colnames(counts1)%in%c("sgRNA","Gene","Category"))],2,ineq::Gini),stringsAsFactors=FALSE)
  # annotation1 <- merge(merge(annotation1,zeros),gini)
  annotation1 <- annotation1[gtools::mixedorder(annotation1$File),]
  # Output ----------------------------------------------------------------------------------------------------------
  # Save read counts and library:
  # utils::write.table(library1,paste0(id,"/fastq_library.txt"),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  utils::write.table(counts1,paste0(id,"/fastq_counts.txt"),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
  # Generate sample annotation1 file in .csv and .xlsx formats:
  annotation1$Replicate <- annotation1$Condition <- ""
  utils::write.csv(annotation1,paste0(id,"/fastq_annotation.csv"),quote=FALSE,row.names=FALSE)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb,"Sample annotation")
  openxlsx::writeData(wb,"Sample annotation",annotation1,rowNames=FALSE,colNames=TRUE)
  openxlsx::setColWidths(wb,"Sample annotation",cols=1:ncol(annotation1),widths="auto")
  openxlsx::freezePane(wb,"Sample annotation",firstRow=TRUE)
  openxlsx::addStyle(wb,"Sample annotation",rows=1,cols=c(1),style=openxlsx::createStyle(fgFill="black",fontColour="white"))
  openxlsx::addStyle(wb,"Sample annotation",rows=1,cols=c(2:4),style=openxlsx::createStyle(fgFill="lightgreen"))
  openxlsx::addStyle(wb,"Sample annotation",rows=1,cols=c(5:6),style=openxlsx::createStyle(fgFill="lightblue"))
  openxlsx::saveWorkbook(wb,paste0(id,"/fastq_annotation.xlsx"), overwrite = TRUE)
  print("Done.")
  if(is.null(shinyF))
    print(paste0("Please fill out the sample information in either one of the following files: ",id,"/fastq_annotation.csv or ",id,"/fastq_annotation.xlsx"))
  return(list(Library=library1,Annotation=annotation1,Counts=counts1))
}
