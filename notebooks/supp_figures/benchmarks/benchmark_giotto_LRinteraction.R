library(zellkonverter)
library(Giotto)
library(microbenchmark)
library(readr)

sce_visium  <-
  readH5AD("/Users/giovanni.palla/.cache/squidpy/visium_R.h5ad")

sce <- sce_visium
gobject <- createGiottoObject(
  raw_exprs = sce@assays@data$X,
  norm_expr = sce@assays@data$X,
  spatial_locs = sce@int_colData$reducedDims$spatial,
  cell_metadata = sce@colData$`cluster`
)

LR_data = data.table::fread(system.file("extdata", "mouse_ligand_receptors.txt", package = 'Giotto'))
LR_data[, ligand_det := ifelse(mouseLigand %in% gobject@gene_ID, T, F)]
LR_data[, receptor_det := ifelse(mouseReceptor %in% gobject@gene_ID, T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$mouseLigand
select_receptors = LR_data_det$mouseReceptor

data.frame(select_ligands) %>% write_csv(
  "/Users/giovanni.palla/Projects/squidpy_reproducibility/notebooks/supp_figures/benchmarks/ligand_mouse_giotto.csv.gz"
)
data.frame(select_receptors) %>% write_csv(
  "/Users/giovanni.palla/Projects/squidpy_reproducibility/notebooks/supp_figures/benchmarks/receptor_mouse_giotto.csv.gz"
)

iter_list <- list(100, 200, 500, 1000)
time_list <- list(1, 1, 1, 1)

bench_ligrec <- function(iterations, t) {
  print(paste0(iterations, " ", t))
  result_ligrec <- microbenchmark(
    interact = exprCellCellcom(
      gobject = gobject,
      verbose = F,
      cluster_column = 'V1',
      random_iter = iterations,
      gene_set_1 = select_ligands,
      gene_set_2 = select_receptors,
    ),
    unit = "s",
    times = t
  )
  
  df <- summary(result_ligrec)
  df$iterations <- iterations
  return(df)
}


result <-
  mapply(bench_ligrec, iter_list, time_list,  SIMPLIFY = F)

final_df <- do.call("rbind", result)

write.csv(
  final_df,
  "/Users/giovanni.palla/Projects/squidpy_reproducibility/notebooks/supp_figures/benchmarks/benchmark_giotto_ligrec.csv"
)
