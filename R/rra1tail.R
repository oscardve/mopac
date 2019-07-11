#' @title One-tail robust rank aggregation
#'
#' @description Performs one-tail a-RRA analysis using the original algorithm of MAGeCK-RRA.
#'
#' @param level1 Vector containing genes.
#' @param level2 Vector containing sgRNAs.
#' @param score Vector containing sgRNA scores.
#'
#' @return A data frame with modified robust-rank-aggregation p-value.
#'
#' @author Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
#' @keywords RRA
#'
#' @references Li, W. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biol. 15, 554 (2014).
#'
#' ##--------------------------------------------------------------------------------------
#'
#' @export
#'

RRA.1tail <- function(level1=NULL,level2=NULL,score=NULL) {
  set.seed(123)
  out <- RRA_1tail(level2,level1,score)
  print("Done.")
  return(out)
}





