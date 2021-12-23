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
    gls <- GwasLocusScores$new("rs4575096", "chr1", tag.snp.loc-shoulder, tag.snp.loc+shoulder)
    checkTrue(all(c("GwasLocusScores", "R6") %in% class(gls)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_eqtls <- function()
{
    tag.snp.chrom <- "chr1"
    tag.snp.loc <- 161185602
    shoulder <- 60000
    gls <- GwasLocusScores$new("rs4575096", "chr1", tag.snp.loc-shoulder, tag.snp.loc+shoulder,
                               tissue.name="GTEx.brain_hippocampus")

    gls$load.eqtls(5e-6)
    tbl.eqtl <- gls$get.tbl.eqtl()
    checkEquals(nrow(tbl.eqtl), 0)

    gls$load.eqtls(1e-4)
    tbl.eqtl <- gls$get.tbl.eqtl()
    dim(tbl.eqtl)
    checkEquals(nrow(tbl.eqtl), 0)

} # test_eqtls
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
