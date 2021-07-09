## ----dependencies, warning=FALSE, message=FALSE-------------------------------
library(bnbc)

## ----dataLoad-----------------------------------------------------------------
data(cgEx)
cgEx

## ----create-------------------------------------------------------------------
cgEx <- ContactGroup(rowData=rowData(cgEx),
                     contacts=contacts(cgEx),
                     colData=colData(cgEx))

## ----print--------------------------------------------------------------------
cgEx

## ----contactTriang------------------------------------------------------------
mat <- matrix(1:9, nrow = 3, ncol = 3)
mat[lower.tri(mat)] <- 0
mat
## Now we fill in the lower triangular matrix with the upper triangular
mat[lower.tri(mat)] <- mat[upper.tri(mat)]
mat

## ----data_to_bnbc, eval=FALSE, echo=TRUE--------------------------------------
#  ## Example not run
#  ## Convert upper triangles to symmetry matrix
#  MatsList <- lapply(upper.mats.list, function(M) {
#      M[lower.tri(M)] <- M[upper.tri(M)]
#  })
#  ## Use ContactGroup constructor method
#  cg <- ContactGroup(rowData = LociData, contacts = MatsList, colData = SampleData)

## ----cooler_get_genome_index--------------------------------------------------
coolerDir <- system.file("cooler", package = "bnbc")
cools <- list.files(coolerDir, pattern="cool$", full.names=TRUE)

step <- 4e4

bin.ixns.list <- bnbc:::getGenomeIdx(cools[1], step)
bin.ixns <- bin.ixns.list$bin.ixns

## ----cooler_get_cg------------------------------------------------------------
data(cgEx)
cool.cg <- bnbc:::getChrCGFromCools(bin.ixns,
                                    files = cools,
                                    chr = "chr22",
                                    colData = colData(cgEx)[1:2,])
all.equal(contacts(cgEx)[[1]], contacts(cool.cg)[[1]])

## ----band_example-------------------------------------------------------------
mat.1 <- contacts(cgEx)[[1]]
mat.1[1000:1005, 1000:1005]
b1 <- band(mat=mat.1, band.no=2)
band(mat=mat.1, band.no=2) <- b1 + 1
mat.1[1000:1005, 1000:1005]

## ----logcpm-------------------------------------------------------------------
cgEx.cpm <- logCPM(cgEx)

## ----smoothing----------------------------------------------------------------
cgEx.smooth <- boxSmoother(cgEx.cpm, h=5)
## or
## cgEx.smooth <- gaussSmoother(cgEx.cpm, radius=3, sigma=4)

## ----bnbc---------------------------------------------------------------------
cgEx.bnbc <- bnbc(cgEx.smooth, batch=colData(cgEx.smooth)$Batch,
                  threshold=1e7, step=4e4, nbands=11, verbose=FALSE)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

