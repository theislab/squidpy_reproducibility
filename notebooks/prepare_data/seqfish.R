library(SingleCellExperiment)
library(zellkonverter)

counts <- readRDS("/Users/giovanni.palla/Datasets/SpatialMouseGastro/counts.rds")
meta <- readRDS("/Users/giovanni.palla/Datasets/SpatialMouseGastro/metadata.rds")

sce <- SingleCellExperiment(list(counts=counts),
                            colData=meta[, 1:12]
)

writeH5AD(sce, file = "/Users/giovanni.palla/Datasets/SpatialMouseGastro/mouse_gastro.h5ad")
