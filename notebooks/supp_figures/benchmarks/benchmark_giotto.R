library(zellkonverter)
library(Giotto)
library(microbenchmark)

sce_seqfish  <-
  readH5AD("/Users/giovanni.palla/.cache/squidpy/seqfish_R.h5ad")
sce_imc  <-
  readH5AD("/Users/giovanni.palla/.cache/squidpy/imc_R.h5ad")
sce_visium  <-
  readH5AD("/Users/giovanni.palla/.cache/squidpy/visium_R.h5ad")
sce_mibitof  <-
  readH5AD("/Users/giovanni.palla/.cache/squidpy/mibitof_R.h5ad")
sce_4i  <-
  readH5AD("/Users/giovanni.palla/.cache/squidpy/fouri_R.h5ad")
sce_slideseqv2  <-
  readH5AD("/Users/giovanni.palla/.cache/squidpy/slideseq_R.h5ad")
sce_merfish  <-
  readH5AD("/Users/giovanni.palla/.cache/squidpy/merfish_R.h5ad")

data_list <- list(sce_imc , sce_visium, sce_mibitof, sce_seqfish, sce_4i, sce_slideseqv2, sce_merfish)
id_list <- list("imc", "visium", "mibitof", "seqfish", "4i", "slideseq", "merfish")
time_list <- list(5,5,5,1,1,1,1)

bench_graph <- function(sce, id, t) {
  print(paste0(id, " ", t))
  gobject <- createGiottoObject(
    raw_exprs = sce@assays@data$X,
    spatial_locs = sce@int_colData$reducedDims$spatial,
    cell_metadata = sce@colData$`cluster`
  )
  results_network <- microbenchmark(
    net_delaunay = createSpatialNetwork(gobject, name = "spatial_network", method = 'Delaunay',),
    net_knn  = createSpatialNetwork(gobject, name = "spatial_network", method = 'kNN',),
    unit = "s",
    times = t
  )
  gobject <- createSpatialNetwork(gobject,
                                  name = "spatial_network",
                                  method = 'Delaunay',)
  results_nhood <- microbenchmark(
    cellproxy = cellProximityEnrichment(
      gobject,
      cluster_column = 'V1',
      spatial_network_name = 'spatial_network',
      adjust_method = 'fdr',
      number_of_simulations = 1000
    ),
    unit = "s",
    times = t
  )
  
  df <- rbind(summary(results_network), summary(results_nhood))
  df$dataset <- id
  return(df)
}


result <- mapply(bench_graph, data_list, id_list,time_list,  SIMPLIFY = F)
final_df <- do.call("rbind", result)
write.csv(final_df, "/Users/giovanni.palla/Projects/squidpy_reproducibility/notebooks/supp_figures/benchmark_giotto_results.csv")
