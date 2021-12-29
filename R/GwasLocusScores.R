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
#' @import MotifDb
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
      tbl.eqtl=NULL,
      eqtl.catalog=NULL,
      motifbreakR.results=NULL
      ),
 #--------------------------------------------------------------------------------
   public = list(

         #--------------------------------------------------------------------------------
         #' @description
         #' Creates a new instance of this class.
         #'
         #' @param tag.snp character, an rsid identifying the locus
         #' @param chrom character, e.g., "chr1"
         #' @param start.loc numeric the start of the genomic region of interest
         #' @param end.loc numeric the end of the genomic region of interest
         #' @param tissue.name character, e.g. "GTEx_V8.Brain_Hippocampus"
         #' @return an object of the GwasLocusScores class

      initialize = function(tag.snp, chrom, start.loc, end.loc, tissue.name){
         stopifnot(length(tissue.name) == 1)
         private$tag.snp <- tag.snp
         private$chrom <- chrom
         private$start.loc <- start.loc
         private$end.loc <- end.loc
         private$advx <- ADvariantExplorer$new(NULL, chrom, start.loc, end.loc)
         private$tissue.name <- tissue.name
         private$eqtl.catalog <- private$advx$geteQTLSummary()
         search.term <- sprintf("^%s$", tissue.name)
         if(length(grep(search.term, private$eqtl.catalog$unique_id)) !=1)
            stop(sprintf("'%s' does not uniquely identify a catalog study", tissue.name))
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' retrieves all eQTLs from the ebi eqtl catalogue
         #'
         #' @returns nothing
      load.eqtls = function(){
         tbl.eqtl <- private$advx$geteQTLsByLocationAndStudyID(private$chrom,
                                                               private$start.loc,
                                                               private$end.loc,
                                                               private$tissue.name,
                                                               method="tabix", simplify=TRUE)
         private$tbl.eqtl <- tbl.eqtl
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' extract all previously obtained eQTLs at or aboce the specified pvalue threshold
         #'
         #' @param pvalue.cutoff numeric, e.g., 1e-4
         #' @returns a data.frame
      get.tbl.eqtl = function(pvalue.cutoff){
         return(subset(private$tbl.eqtl, pvalue <= pvalue.cutoff))
         },

         #--------------------------------------------------------------------------------
         #' @description
         #' discover which motifs are broken by the eQTL variants, for specified gene and pval
         #'
         #' @param targetGene character, a gene symbol in the eQTL table
         #' @param pvalue.cutoff numeric, e.g., 1e-4
         #' @returns a data.frame
      breakMotifsAtEQTLs = function(targetGene, pvalue.cutoff=NA){
         tbl.sub <- subset(private$tbl.eqtl, gene==targetGene)
         if(!is.na(pvalue.cutoff))
             tbl.sub <- subset(tbl.sub, pvalue <= pvalue.cutoff)
         rsids <- tbl.sub$rsid
         browser()
         motifs.selected <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco-core-A"))
         snps.gr <- snps.from.rsid(rsid = rsids,
                              dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                              search.genome=BSgenome.Hsapiens.UCSC.hg38)
         bpparam <- MulticoreParam(workers=4)
         private$motifbreakR.results <- motifbreakR(snpList = snps.gr,
                                                    filterp = TRUE,
                                                    pwmList = motifs.selected,
                                                    show.neutral=FALSE,
                                                    method = c("ic", "log", "notrans")[1],
                                                    bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                                    BPPARAM = bpparam,
                                                    verbose=TRUE)
         tbl.breaks <- as.data.frame(private$motifbreakR.results, row.names=NULL)
         tbl.breaks <- subset(tbl.breaks, effect=="strong")
         tbl.breaks$pctDelta <- with(tbl.breaks, pctAlt - pctRef)
         return(invisible(tbl.breaks))
         } # breakMotifsAtEQTLs

   ) # public

 ) # class GwasLocusScores
#--------------------------------------------------------------------------------



