library(zellkonverter)
library(Giotto)

sce  <- readH5AD("/Users/giovanni.palla/.cache/squidpy/imc.h5ad")

gobject <- createGiottoObject(raw_exprs = sce@assays@data$X, 
                         spatial_locs = sce@int_colData$reducedDims$spatial,
                         cell_metadata = sce@colData$`cell type`)
  
gobject <- createSpatialNetwork(gobject,
                     name = "spatial_network",
                     method = 'Delaunay',)

cell_proximities <- cellProximityEnrichment(
  gobject,
  cluster_column = 'V1',
  spatial_network_name = 'spatial_network',
  adjust_method = 'fdr',
  number_of_simulations = 1000
)

cellProximityHeatmap(gobject, cell_proximities)
