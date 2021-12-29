library(RUnit)
library(GwasLocusScores)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_eqtls()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    tag.snp.chrom <- "chr1"
    tag.snp.loc <- 161185602
    shoulder <- 60000
    gls <- GwasLocusScores$new("rs4575096", "chr1", tag.snp.loc-shoulder, tag.snp.loc+shoulder,
                               "GTEx_V8.Brain_Hippocampus")
    checkTrue(all(c("GwasLocusScores", "R6") %in% class(gls)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_eqtls <- function()
{
    tag.snp.chrom <- "chr1"
    tag.snp.loc <- 161185602
    shoulder <- 60000

    checkException(
        gls <- GwasLocusScores$new("rs4575096", "chr1", tag.snp.loc-shoulder, tag.snp.loc+shoulder,
                                   tissue.name="intentional.error"),
        silent=TRUE
        )


    gls <- GwasLocusScores$new("rs4575096", "chr1", tag.snp.loc-shoulder, tag.snp.loc+shoulder,
                               tissue.name="GTEx_V8.Brain_Hippocampus")

    gls$load.eqtls()
    tbl.eqtl <- gls$get.tbl.eqtl(5e-4)
    checkTrue(nrow(tbl.eqtl) > 50)
    checkTrue(nrow(tbl.eqtl) < 100)
    checkEquals(colnames(tbl.eqtl), c("rsid", "pvalue", "gene", "total.alleles", "beta", "id"))
    gene.count <- length(unique(tbl.eqtl$gene))
    checkTrue(gene.count > 5 & gene.count < 10)  # 8 on 28 dec 2021)

    tbl.eqtl <- gls$get.tbl.eqtl(1e-5)
    checkTrue(nrow(tbl.eqtl) > 5 & nrow(tbl.eqtl) <= 10) # 7 on 28 dec 2021

} # test_eqtls
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
