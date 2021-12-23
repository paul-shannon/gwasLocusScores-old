#' @title GwasLocusScores
#' @description assess RNA-binding protein binding sites, their structure and motif matches
#' @name GwasLocusScores

#' @import rtracklayer
#' @import igvR
#' @import TrenaMultiScore
#' @import TrenaProjectAD
#' @import EndophenotypeExplorer
#' @import ADvariantExplorer
#' @import trena
#' @import motifbreakR
#' @import SNPlocs.Hsapiens.dbSNP151.GRCh38
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import BiocParallel
#' @import ghdb
#'
#' @field rbp the gene symbol name of the protein
#'
#' @examples
#'   gls <- GwasLocusScores("rs4575096")
#' @export
#----------------------------------------------------------------------------------------------------
 library(R6)
 GwasLocusScores = R6Class("GwasLocusScores",
 #--------------------------------------------------------------------------------
   private = list(
      tag.snp=NULL,
      chrom=NULL,
      start.loc=NULL,
      end.loc=NULL,
      tissue.name=NULL,
      advx=NULL,      # ADvariantExplorer
      etx=NULL,        # EndophenotypeExplorer
      tbl.eqtl=NULL
      ),
 #--------------------------------------------------------------------------------
   public = list(
         #' @description
         #' Creates a new instance of this class.
         #'
         #' @param tag.snp character, an rsid identifying the locus
         #' @param chrom character, e.g., "chr1"
         #' @param start.loc numeric the start of the genomic region of interest
         #' @param end.loc numeric the end of the genomic region of interest
         #' @param tissue.name character, e.g. "GTEx.brain_hippocampus"
         #' @param

      initialize = function(tag.snp, chrom, start.loc, end.loc, tissue.name){
         private$tag.snp <- tag.snp
         private$chrom <- chrom
         private$start.loc <- start.loc
         private$end.loc <- end.loc
         private$advx <- ADvariantExplorer$new(NULL, chrom, start.loc, end.loc)
         private$tissue.name <- tissue.name
         },
      load.eqtls = function(pvalue.cutoff){
         tbl.eqtl <- private$advx$geteQTLsByLocationAndStudyID(private$chrom,
                                                               private$start.loc,
                                                               private$end.loc,
                                                               private$tissue.name,
                                                               method="tabix", simplify=TRUE)
         if(nrow(tbl.eqtl) > 0)
           tbl.eqtl <- subset(tbl.eqtl, pvalue <= pvalue.cutoff)
         private$tbl.eqtl <- tbl.eqtl
         },
      get.tbl.eqtl = function(){
         return(invisible(private$tbl.eqtl))
         }

   ) # public

 ) # class GwasLocusScores
#--------------------------------------------------------------------------------



