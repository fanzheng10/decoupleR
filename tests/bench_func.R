# 1. Load Prerequisites ----
library(here)
library(ggplot2)
library(pheatmap)

opts <- list(
  scira = list(),
  pscira = list(),
  mean = list(.likelihood = NULL),
  viper = list(options = list(verbose = FALSE, minsize=0)),
  gsva = list(options = list(verbose = FALSE))
)

# 2. Run benchmark ----
bench_test <- load_design(here("inst/testdata/inputs/",
                               "design_tibble.json")) %>%
  bench_couple(., opts) %>%
  bench_format() %>%
  mutate(roc = activity %>% map(calc_roc_curve))



# 3. Get Summary and Plots ----
bench_summary <- bench_sumplot(bench_test,
                                 title = "",
                                 coverage = TRUE)
bench_summary$auroc_heat

# 4. Save Benchmark ----
# saveRDS(bench_test, here("inst/testdata/outputs", "chEA3_cbd_results.rds"))


