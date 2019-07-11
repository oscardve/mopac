#' Sample sgRNA library.
#'
#' A sample sgRNA library.
#'
#' @format A data frame containing 74768 rows (one per sgRNA) and the following 3 columns:
#' \itemize{
#'   \item sgRNA: name of the sgRNAs.
#'   \item Gene: gene names for each sgRNA.
#'   \item ID: Hugo IDs for each gene.
#' }
#'
#' @references Aguirre, A. J. et al. Genomic Copy Number Dictates a Gene-Independent Cell Response to CRISPR/Cas9 Targeting. Cancer Discov. 6, 914–929 (2016).
#' @source https://portals.broadinstitute.org/achilles
"sgRNA_library"

#' Sample table of read counts from the Avana dataset.
#'
#' The raw integer read counts of the sgRNAs for two conditions from the Avana dataset.
#'
#' @format A data frame with 74687 rows (one per sgRNA) and the following 10 columns:
#' \itemize{
#'   \item sgRNA: name of the sgRNAs
#'   \item pDNA_0: plasmid of batch 0 (initial condition of DAN-G).
#'   \item pDNA_2: plasmid of batch 2 (initial condition of CCK-81).
#'   \item DAN-G-311Cas9 Rep A p6: replicate A of 1st condition.
#'   \item DAN-G-311Cas9 Rep B p6: replicate B of 1st condition.
#'   \item DAN-G-311Cas9 Rep C p6: replicate C of 1st condition.
#'   \item DAN-G-311Cas9 Rep D p6: replicate D of 1st condition.
#'   \item CCK-81-311cas9 Rep A p2: replicate A of 2nd condition.
#'   \item CCK-81-311cas9 Rep B p2: replicate B of 2nd condition.
#'   \item CCK-81-311cas9 Rep C p2: replicate C of 2nd condition.
#' }
#'
#' @references Meyers, R. M. et al. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nat Genet (2017). doi:10.1038/ng.3984
#' @source https://portals.broadinstitute.org/achilles
"dang_cck81"

#' Annotation for the sample table of read counts from the Avana dataset.
#'
#' An example of an annotation file for the sample from the Avana dataset.
#'
#' @format A data frame with 9 rows (one per sample) and the following 8 columns:
#' \itemize{
#'   \item Sample: name of the sample, which is equal to the column name of the dataset.
#'   \item Zeros: number of sgRNAs with zero counts in the corresponding sample.
#'   \item Gini: gini index of the corresponding sample.
#'   \item Cell_line: name of the cell line for the corresponding sample.
#'   \item Treatment: name of the treatment for the corresponding sample.
#'   \item Timepoint: name of the timepoint for the corresponding sample.
#'   \item Replicate: name of the replicate for the corresponding sample.
#'   \item Batch: name of the batch for the corresponding sample.
#' }
#'
"dang_cck81_annotation"

#' STRINGdb datasets
#'
#' @format A list of data frames containing protein interactions in homo sapiens downloaded from the STRING consortium
#'
#' @references Szklarczyk, D. et al. The STRING database in 2017: Quality-controlled protein-protein association networks, made broadly accessible. Nucleic Acids Res. 45, D362–D368 (2017).
#' @source http://www.string-db.org
"stringdb_9606"

#' Hart panessential genes used as positive controls.
#'
#' The names of the panessential genes used for categorization.
#'
#' @format A character vector
#'
#' @references Hart, T., Brown, K. R., Sircoulomb, F., Rottapel, R. & Moffat, J. Measuring error rates in genomic perturbation screens: gold standards for human functional genomics. Mol. Syst. Biol. 10, (2014).
#' @source https://www.ncbi.nlm.nih.gov/pubmed/24987113
"essentials"

#' Hart non-essential genes used as negative controls.
#'
#' The names of the non-essential genes used for categorization.
#'
#' @format A character vector
#'
#' @references Hart, T., Brown, K. R., Sircoulomb, F., Rottapel, R. & Moffat, J. Measuring error rates in genomic perturbation screens: gold standards for human functional genomics. Mol. Syst. Biol. 10, (2014).
#' @source https://www.ncbi.nlm.nih.gov/pubmed/24987113
"nonessentials"
