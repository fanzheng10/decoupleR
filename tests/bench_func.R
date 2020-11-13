# 1. Load Prerequisites ----
library(here)
library(ggplot2)
library(pheatmap)
library(jsonlite)

# options
opts <- list(
  scira = list(),
  pscira = list(),
  mean = list(.likelihood = NULL),
  viper = list(options = list(verbose = FALSE, minsize=0)),
  gsva = list(options = list(verbose = FALSE))
)

# network locations
dorothea_loc <- here("inst/testdata/inputs", "input_dorothea_full.rds")
chea3_loc <- here("inst/testdata/inputs/ChEA3", "ChEA3_realTFs_only.RDS")
regnet_loc  <- here("inst/testdata/inputs/regnetwork", "regnework_realTFs_only.rds")


# statitics
statistics <- c(
  "scira",
  "pscira",
  "mean",
  "viper",
  "gsva"
)

# benchmark datasets
knock_expr <- here("inst/testdata/inputs", "knockTF_gene_expr.rds")
knock_meta <- here("inst/testdata/inputs", "knockTF_meta.rds")

dbm_expr <- here("inst/testdata/inputs", "input-dorothea_bench_expr.rds")
dbm_meta <- here("inst/testdata/inputs", "input-dorothea_bench_meta.rds")



# 2. Prepare design ----
design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta,
  "dorothea_knock", dorothea_loc, "A", "tf", "target", statistics, knock_expr, knock_meta,
  "dorothea_dbm", dorothea_loc, "A", "tf", "target", statistics, dbm_expr, dbm_meta,
  "chEA3_knock", chea3_loc, "remap_chip_seq", "tf", "target", statistics, knock_expr, knock_meta,
  "chEA3_dbm", chea3_loc, "remap_chip_seq", "tf", "target", statistics, dbm_expr, dbm_meta,
  "regnet_knock", regnet_loc, "High", "tf", "target", statistics, knock_expr, knock_meta,
  "regnet_dbm", regnet_loc, "High", "tf", "target", statistics, dbm_expr, dbm_meta,
)
design_tibble

write(toJSON(design_tibble, pretty=T),
      here("inst/testdata/inputs/", "design_tibble.json"))




# 3. Run benchmark ----
bench_test <- load_design(here("inst/testdata/inputs/",
                               "design_tibble.json")) %>%
  bench_couple(., opts) %>%
  bench_format()  # %>%
  # mutate(roc = activity %>% map(calc_roc_curve))

# saveRDS(bench_test, here("inst/testdata/outputs", "full_run_test.rds"))

# Estimate ROC
bench_roc <- bench_test[-c(14:17),] %>%
  mutate(roc = activity %>% map(calc_roc_curve))

bench_test[14,]$activity

bench_test[14:17,] %>%
  mutate(roc = activity %>% map(calc_roc_curve))

# bench_test[13,]$activity
# bench_test[14,]$activity


# 4. Get Summary and Plots ----
bench_test

bench_summary <- bench_sumplot(bench_roc,
                                 title = "",
                                 coverage = TRUE)
bench_summary$auroc

?pheatmap

# 5. Save Benchmark ----
# saveRDS(bench_test, here("inst/testdata/outputs", "chEA3_cbd_results.rds"))


