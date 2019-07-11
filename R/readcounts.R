#' @title Map sgRNA read counts to an sgRNA library.
#'
#' @description Maps a table of read counts to an sgRNA library and generates a sample annotation file if not provided.
#'
#' @param id an optional identifier for the analysis. All input and output will be stored in a directory with this name. No output files will be produced if not supplied.
#' @param library an optional sgRNA library. It may be a data frame or a file path. It must contain a column called \code{"sgRNA"} and another called \code{"Gene"}. It may also contain a \code{"Category"} column. Only necessary if the table of counts does not include a column called "Gene".
#' @param counts a data frame or file path containing the raw or normalized read counts. It should include a column named "sgRNA".  It may also contain a \code{"Gene"} and/or \code{"Category"} column.
#' @param annotation an optional data frame or file path containing the FASTQ file or condition annotation information. Defaults to \code{NULL}.
#' @param essential an optional vector of gene names to add genes categorized as essential.
#' @param nonessential an optional vector of gene names to add genes categorized as non-essential.
#' @param nontargeting an optional vector of "gene" names to add the default set of "genes" categorized as non-targeting.
#' @param shinyF logical. An internal parameter used by the graphical user interface. If using command-line mode leave as \code{NULL}.
#'
#' @details
#' The following files are produced in a directory with name equal to the analysis' identifier:
#' \itemize{
#'  \item{\code{library.txt}, containing two columns \code{"sgRNA"}, \code{"Gene"} (and possibly also \code{"Category"}) for each sgRNA mapped to the library file.}
#'  \item{\code{annotation.csv} and \code{annotation.xlsx}, two annotation files, either of which may be chosen to fill out the following sample information: condition, replicate and control.}
#'  \item{\code{counts.txt}, containing the read counts per sgRNA for each sample.}
#' }
#'
#' @return A list containing the following components:
#' \item{Library}{The sgRNA library.}
#' \item{Annotation}{The table of sample annotation.}
#' \item{Counts}{The table of read counts at the sample level.}
#'
#' @author Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
#' @keywords read counts
#'
#' @examples
#'
#' ## Example I. Supplying a table of read counts as a data frame:
#' ## It must include a column named "sgRNA":
#' dataset <- MoPAC::dang_cck81
#' library <- MoPAC::sgRNA_library
#' reads <- read.counts(id="test", counts=dataset, library=library)
#'
#' ## Example II. Supplying the path of the read counts file:
#' ## It may be in .xlsx, .csv or tab-separted format:
#' utils::write.table(dataset, "dataset.txt", sep='\t', quote=FALSE, row.names=FALSE)
#' utils::write.table(library, "library.txt", sep='\t', quote=FALSE, row.names=FALSE)
#' reads <- read.counts(id="test", counts="dataset.txt", library="library.txt")
#'
#' ## Example III. Without using an identifier to generate output:
#' reads <- read.counts(counts = dataset, library=library)
#'
#' ##--------------------------------------------------------------------------------------
#'
#' @export
#'

read.counts <- function(id=NULL,library=NULL,counts=NULL,annotation=NULL,essential=NULL,nonessential=NULL,nontargeting=NULL,shinyF=NULL) {
  # Read counts -------------------------------------------------------------------------
  print("Importing data..."); if(!is.null(shinyF)) shinyF(0.1,"Input","Importing data...")
  if(is.null(id) & is.null(counts)) { # insufficient information provided
    stop("Please specify a working directory or provide a table of read counts")
  } else if(!is.null(id) & is.null(counts)) { # is a table of read counts available in the working directory?
    if(file.exists(paste0(id,"/counts.txt"))) {
      counts0 <- utils::read.table(paste0(id,"/counts.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else if(file.exists(paste0(id,"/fastq_counts.txt"))) {
      counts0 <- utils::read.table(paste0(id,"/fastq_counts.txt"),header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Please provide a table of read counts or run the FASTQ alignment module first.")
  } else if(is.data.frame(counts)) { # was a data frame provided as a table of read counts?
    counts0 <- counts
    counts0[,sapply(counts0,is.factor)] <- sapply(counts0[,sapply(counts0,is.factor)],as.character)
    # colnames(counts0) <- stringr::str_replace_all(colnames(counts0),"[^[:alnum:]]","") #remove non-alphanumeric characters in the sample names
  } else if(is.character(counts)) { # was a file path provided as a table of read counts?
    if(file.exists(counts)) {
      if(tools::file_ext(counts)=="xlsx") {
        counts0 <- openxlsx::read.xlsx(counts)
      } else if(tools::file_ext(counts)=="csv") {
        counts0 <- utils::read.csv(counts,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"")
      } else counts0 <- utils::read.table(counts,header=TRUE,sep="\t",as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    } else stop("Counts file not found in the specified address.")
    # colnames(counts0) <- stringr::str_replace_all(colnames(counts0),"[^[:alnum:]]","") #remove non-alphanumeric characters in the sample names
  } else stop("Please provide a valid table of read counts.")
  colnames(counts0)[grep("sgRNA|sgrna|SGRNA",colnames(counts0))] <- "sgRNA"
  colnames(counts0)[grep("Gene|gene|GENE",colnames(counts0))] <- "Gene"
  colnames(counts0)[grep("Category|category|CATEGORY",colnames(counts0))] <- "Category"
  if(!"sgRNA"%in%colnames(counts0)) stop("Please make sure the counts contain a column called 'sgRNA'.")
  # Read library file if available --------------------------------------------------------
  library0 <- NULL
  if(!is.null(id) & is.null(library)) {
    if(file.exists(paste0(id,"/library.txt")))
      library0 <- utils::read.table(paste0(id,"/library.txt"),header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    else if(file.exists(paste0(id,"/fastq_library.txt")))
      library0 <- utils::read.table(paste0(id,"/fastq_library.txt"),header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
  } else if(!is.null(library)) {
    if(is.data.frame(library)) {
      library0 <- library
      library0[,sapply(library0,is.factor)] <- sapply(library0[,sapply(library0,is.factor)],as.character)
    } else if(is.character(library)) {
      if(file.exists(library)) {
        if(tools::file_ext(library)=="xlsx") {
          library0 <- openxlsx::read.xlsx(library)
        } else if(tools::file_ext(library)=="csv") {
          library0 <- utils::read.csv(library,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"")
        } else library0 <- utils::read.table(library,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
      } else stop("sgRNA library not found in the specified address.")
    }
  }
  # Merge counts and library if appropriate ---------------------------------------------
  if(!is.null(library0)) {
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
    if(!"Gene"%in%colnames(counts0)) counts0 <- merge(library0[,c("sgRNA","Gene")],counts0)
    if(!"Category"%in%colnames(counts0) & "Category"%in%colnames(library0)) counts0 <- merge(library0[,c("sgRNA","Gene","Category")],counts0)
    counts0 <- cbind(counts0[,which(colnames(counts0)%in%c("sgRNA","Gene","Category")),drop=F],counts0[,-which(colnames(counts0)%in%c("sgRNA","Gene","Category")),drop=F])
  }
  if(nrow(counts0)==0) stop("Please make sure that the columns in the sgRNA library are correctly labeled.")
  if(!"Gene"%in%colnames(counts0)) stop("Please make sure that either the read counts or the sgRNA library contain a column called 'Gene'.")
  counts0 <- counts0[order(counts0$Gene),]
  library0 <- counts0[,which(colnames(counts0)%in%c("sgRNA","Gene","Category"))]
  # Read annotation if available ------------------------------------------------------
  annotation0 <- NULL
  if(!is.null(id) & is.null(annotation)) {
    if(file.exists(paste0(id,"/counts.txt"))) {
      annotation0 <- utils::read.csv(paste0(id,"/annotation.csv"),header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
      if(sum(is.na(annotation0))>0)
        annotation0 <- openxlsx::read.xlsx(paste0(id,"/annotation.xlsx"),check.names=FALSE,na.strings="")
      # if(sum(is.na(annotation0))>0)
      #   stop("Please fill out the sample file annotation before proceeding.")
    } else if(file.exists(paste0(id,"/fastq_counts.txt"))) {
      annotation0 <- utils::read.csv(paste0(id,"/fastq_annotation.csv"),header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
      if(sum(is.na(annotation0))>0)
        annotation0 <- openxlsx::read.xlsx(paste0(id,"/fastq_annotation.xlsx"),check.names=FALSE,na.strings="")
      if(sum(is.na(annotation0))>0)
        stop("Please fill out the FASTQ file annotation before proceeding.")
    }
  } else if(!is.null(annotation)) {
    if(is.data.frame(annotation)) {
      annotation0 <- annotation
      annotation0[,sapply(annotation0,is.factor)] <- sapply(annotation0[,sapply(annotation0,is.factor)],as.character)
    } else if(is.character(annotation)) {
      if(file.exists(annotation)) {
        if(tools::file_ext(annotation)=="xlsx") {
          annotation0 <- openxlsx::read.xlsx(annotation,check.names=FALSE,na.strings="")
        } else if(tools::file_ext(annotation)=="csv") {
          annotation0 <- utils::read.csv(annotation,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
        } else annotation0 <- utils::read.table(annotation,header=TRUE,as.is=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"",na.strings="")
      } else stop("Annotation file not found in the specified address.")
      if(sum(is.na(annotation0))>0)
        stop("Please fill out the FASTQ file annotation before proceeding.")
    } else stop("The annotation object provided could not be read.")
  }
  # Merge fastq files if appropriate -----------------------------------------------------
  if(!is.null(annotation0)) {
    if("File"%in%colnames(annotation0)) {
      if(sum(annotation0$File%in%colnames(counts0))==nrow(annotation0)) {
        groups <- split(annotation0$File,annotation0[,c("Condition","Replicate")],drop=T)
        counts0 <- cbind(counts0[,which(colnames(counts0)%in%c("sgRNA","Gene","Category")),drop=F],do.call(cbind,lapply(groups,function(x)rowSums(counts0[,x,drop=F]))))
        annotation0$Sample <- paste0(annotation0$Condition,".",annotation0$Replicate)
        annotation0 <- unique(annotation0[,c("Sample","Condition","Replicate")])
        annotation0$Control <- ""
      } else stop("Mismatch found in the FASTQ file annotation, please re-start the analysis in an empty working directory.")
    } else if(!"Sample"%in%colnames(annotation0)) stop("Mismatch found in the sample annotation, please re-start the analysis in an empty working directory.")
  }
  # Generate new sample annotation ----------------------------------------------------------
  if(is.null(annotation0)) {
    annotation0 <- data.frame(Sample=colnames(counts0)[-which(colnames(counts0)%in%c("sgRNA","Gene","Category"))],stringsAsFactors=FALSE)
    annotation0$Control <- annotation0$Replicate <- annotation0$Condition <- ""
  }
  annotation0$Depth <- apply(counts0[,annotation0$Sample,drop=F],2,sum)
  annotation0$Zeros <- apply(counts0[,annotation0$Sample,drop=F],2,function(y)sum(y==0))
  annotation0$Gini <- apply(counts0[,annotation0$Sample,drop=F],2,ineq::Gini)
  annotation0 <- annotation0[,c("Sample","Depth","Zeros","Gini","Condition","Replicate","Control")]
  annotation0[is.na(annotation0)] <- ""
  # Output -----------------------------------------------------------------------------------
  if(!is.null(id)) {
    print("Saving output..."); if(!is.null(shinyF)) shinyF(0.1,"Input","Saving output...")
    # Create directory:
    dir.create(id,showWarnings=FALSE)
    # Save read counts and library:
    utils::write.table(counts0,paste0(id,"/counts.txt"),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
    utils::write.table(library0,paste0(id,"/library.txt"),sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE)
    # Generate sample annotation file in .csv and .xlsx formats:
    utils::write.csv(annotation0,paste0(id,"/annotation.csv"),quote=FALSE,row.names=FALSE)
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb,"Sample annotation")
    openxlsx::writeData(wb,"Sample annotation",annotation0,rowNames=FALSE,colNames=TRUE)
    openxlsx::setColWidths(wb,"Sample annotation",cols=1:ncol(annotation0),widths="auto")
    openxlsx::freezePane(wb,"Sample annotation",firstRow=TRUE)
    openxlsx::addStyle(wb,"Sample annotation",rows=1,cols=c(1),style=openxlsx::createStyle(fgFill="black",fontColour="white"))
    openxlsx::addStyle(wb,"Sample annotation",rows=1,cols=c(2:4),style=openxlsx::createStyle(fgFill="lightgreen"))
    openxlsx::addStyle(wb,"Sample annotation",rows=1,cols=c(5:7),style=openxlsx::createStyle(fgFill="lightblue"))
    openxlsx::saveWorkbook(wb,paste0(id,"/annotation.xlsx"), overwrite = TRUE)
    print("Done.")
    print(paste0("Please fill out the sample information in either one of the following files: ",id,"/annotation.csv or ",id,"/annotation.xlsx"))
  } else print("Done.")
  return(list(Library=library0,Annotation=annotation0,Counts=counts0))
}
