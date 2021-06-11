library(zellkonverter)
library(SpatialExperiment)

dir <- system.file(
  file.path("extdata", "10xVisium", "section1"),
  package = "SpatialExperiment")

# read in counts
fnm <- file.path(dir, "raw_feature_bc_matrix")
sce <- DropletUtils::read10xCounts(fnm)

# read in image data
img <- readImgData(
  path = file.path(dir, "spatial"),
  sample_id="foo")

# read in spatial coordinates
fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
xyz <- read.csv(fnm, header = FALSE,
                col.names = c(
                  "barcode", "in_tissue", "array_row", "array_col",
                  "pxl_row_in_fullres", "pxl_col_in_fullres"))

# construct observation & feature metadata
rd <- S4Vectors::DataFrame(
  symbol = rowData(sce)$Symbol)

# construct 'SpatialExperiment'
(spe <- SpatialExperiment(
  assays = list(counts = assay(sce)),
  colData = colData(sce), rowData = rd, imgData = img,
  spatialData=DataFrame(xyz),
  spatialCoordsNames=c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  sample_id="foo"))

dir <- system.file(
  file.path("extdata", "10xVisium"),
  package = "SpatialExperiment")

sample_ids <- c("section1", "section2")
samples <- file.path(dir, sample_ids)

(spe <- read10xVisium(samples, sample_ids,
                      type = "sparse", data = "raw",
                      images = "lowres", load = FALSE))

writeH5AD(spe, "/Users/giovanni.palla/Projects/squidpy_reproducibility/notebooks/reproducibility/spa_exp.h5ad")

# try to load anndata with spatial data

ad <- readH5AD("/Users/giovanni.palla/.cache/squidpy/visium_hne_adata_crop.h5ad")
# doesn't work, throw error
