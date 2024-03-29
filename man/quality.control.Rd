% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality.R
\name{quality.control}
\alias{quality.control}
\title{Pre-processing and quality assessment of CRISPR screen.}
\usage{
quality.control(id = NULL, counts = NULL, annotation = NULL,
  pseudocount = 0.15, report = FALSE, shinyF = NULL)
}
\arguments{
\item{id}{an optional identifier for the analysis. All input and output will be stored in a directory with this name. No output files will be produced if not supplied.}

\item{counts}{an optional data frame or file path containing the read count information, which may be in either integer or normalized form. If \code{NULL}, the counts file in the id directory will be read. Defaults to \code{NULL}.}

\item{annotation}{an optional data frame or file path containing the sample annotation information. If \code{NULL}, the annotation file in the id directory will be read. Defaults to \code{NULL}.}

\item{pseudocount}{numeric. An optional pseudo-count for log2 transformation. If \code{NULL}, default optimization is performed over the values 0.01, 0.05, 0.1, 0.15, 0.2, 0.25 and 0.3. Defaults to 0.15.}

\item{report}{logical. Should a pdf report of quality assessment be generated? Requires the id to not be \code{NULL}, and a pandoc installation. Defaults to \code{FALSE}.}

\item{shinyF}{logical. An internal parameter used by the graphical user interface. If using command-line mode leave as \code{NULL}.}
}
\value{
A list containing the following components:
\item{Counts}{raw integer counts.}
\item{Pseudo}{optimal pseudocount and pseudocount versus variance function.}
\item{Annotation}{a data frame with sample annotation information, which is merged for fastq files belonging to the same sample.}
\item{Annotation1}{a data frame with sample annotation information, not including the control (i.e. day-0) conditions.}
\item{sgRNA_log_reps}{the log-counts obtained through the optimal or user-defined pseudocount.}
\item{sgRNA_lfc_reps}{log-fold-changes for each sgRNA of each replicate.}
\item{sgRNA_lfc}{log-fold-changes for the repliate-averaged sgRNAs.}
}
\description{
Pre-processes and assesses the quality of a CRISPR screen by computing log-fold-changes with optimized pseudo-count, and generates a pdf report.
}
\details{
Log-fold-changes at the sgRNA level are computed before and after averaging the replicates. If an identifier is supplied, the following files are produced in the directory with name equal to the identifier:
\itemize{
 \item{\code{sgRNA_log2.txt}, the log2-transformed read counts normalized by sequencing depth.}
 \item{\code{sgRNA_lfc_reps.txt} and \code{sgRNA_lfc.txt}, the log-fold-changes before and after averaging the replicates.}
 \item{\code{qc_variance.txt}, \code{qc_weights.txt}, \code{qc_scaling.txt} and \code{qc_shifting.txt}, the optimized preliminary regression parameters employing the unfiltered list of control genes whose size is equal to the most frequent gene size in the dataset.}
 \item{\code{qc.pdf}, a report for the visualization of the quality assessment.}
}
}
\examples{

## Example I. Pre-processing read counts without producing output:

## 1. Prepare counts as a data frame:
counts <- read.counts(counts = MoPAC::dang_cck81, library=MoPAC::sgRNA_library)

## 2. Prepare sample annotation as a data frame:
annotation <- counts$Annotation
annotation$Condition <- c("Plasmid0","Plasmid1",rep("DANG",4),rep("CCK81",3))
annotation$Replicate <- c("A","A","A","B","C","D","A","B","C")
annotation$Control <- c("","",rep("Plasmid0",4),rep("Plasmid1",3))

## 3. Get sgRNA-level and gene-level fold changes:
qc <- quality.control(counts=counts$Counts, annotation=annotation)


##--------------------------------------------------------------------------------------

}
\author{
Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
}
\keyword{control}
\keyword{quality}
