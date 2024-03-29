% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readcounts.R
\name{read.counts}
\alias{read.counts}
\title{Map sgRNA read counts to an sgRNA library.}
\usage{
read.counts(id = NULL, library = NULL, counts = NULL,
  annotation = NULL, essential = NULL, nonessential = NULL,
  nontargeting = NULL, shinyF = NULL)
}
\arguments{
\item{id}{an optional identifier for the analysis. All input and output will be stored in a directory with this name. No output files will be produced if not supplied.}

\item{library}{an optional sgRNA library. It may be a data frame or a file path. It must contain a column called \code{"sgRNA"} and another called \code{"Gene"}. It may also contain a \code{"Category"} column. Only necessary if the table of counts does not include a column called "Gene".}

\item{counts}{a data frame or file path containing the raw or normalized read counts. It should include a column named "sgRNA".  It may also contain a \code{"Gene"} and/or \code{"Category"} column.}

\item{annotation}{an optional data frame or file path containing the FASTQ file or condition annotation information. Defaults to \code{NULL}.}

\item{essential}{an optional vector of gene names to add genes categorized as essential.}

\item{nonessential}{an optional vector of gene names to add genes categorized as non-essential.}

\item{nontargeting}{an optional vector of "gene" names to add the default set of "genes" categorized as non-targeting.}

\item{shinyF}{logical. An internal parameter used by the graphical user interface. If using command-line mode leave as \code{NULL}.}
}
\value{
A list containing the following components:
\item{Library}{The sgRNA library.}
\item{Annotation}{The table of sample annotation.}
\item{Counts}{The table of read counts at the sample level.}
}
\description{
Maps a table of read counts to an sgRNA library and generates a sample annotation file if not provided.
}
\details{
The following files are produced in a directory with name equal to the analysis' identifier:
\itemize{
 \item{\code{library.txt}, containing two columns \code{"sgRNA"}, \code{"Gene"} (and possibly also \code{"Category"}) for each sgRNA mapped to the library file.}
 \item{\code{annotation.csv} and \code{annotation.xlsx}, two annotation files, either of which may be chosen to fill out the following sample information: condition, replicate and control.}
 \item{\code{counts.txt}, containing the read counts per sgRNA for each sample.}
}
}
\examples{

## Example I. Supplying a table of read counts as a data frame:
## It must include a column named "sgRNA":
dataset <- MoPAC::dang_cck81
library <- MoPAC::sgRNA_library
reads <- read.counts(id="test", counts=dataset, library=library)

## Example II. Supplying the path of the read counts file:
## It may be in .xlsx, .csv or tab-separted format:
utils::write.table(dataset, "dataset.txt", sep='\\t', quote=FALSE, row.names=FALSE)
utils::write.table(library, "library.txt", sep='\\t', quote=FALSE, row.names=FALSE)
reads <- read.counts(id="test", counts="dataset.txt", library="library.txt")

## Example III. Without using an identifier to generate output:
reads <- read.counts(counts = dataset, library=library)

##--------------------------------------------------------------------------------------

}
\author{
Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
}
\keyword{counts}
\keyword{read}
