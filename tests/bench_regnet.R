library(here)
# Load Prerequisites ----
# Locations for the network resource, benchmark gene expression, and meta data
net_loc <- here("inst/testdata/inputs/regnetwork", "regnework_realTFs_only.rds")
# bench_expr <- here("inst/testdata/inputs", "knockTF_ex.rds")
# bench_meta <- here("inst/testdata/inputs", "knockTF_meta.rds")
bench_expr <- here("inst/testdata/inputs", "input-dorothea_bench_example.rds")
bench_meta <- here("inst/testdata/inputs", "input-dorothea_bench_meta.rds")

# Libraries/Confidence levels in each gene set (or anything that )
regs <- c("High", "Medium", "Low")

# Statistics
statistics <- c(
 "scira",
 "pscira",
  "mean",
  "viper",
  "gsva"
)

# Define options for each statistic
opts <- list(
  scira = list(),
  pscira = list(),
  mean = list(.likelihood = NULL),
  viper = list(options = list(verbose = FALSE, minsize=4)),
  gsva = list(options = list(verbose = FALSE))
)

# Design
design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "regnet", net_loc, regs[1], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "regnet", net_loc, regs[2], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "regnet", net_loc, regs[3], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "regnet", net_loc, regs[1:2], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "regnet", net_loc, regs[1:3], "tf", "target", statistics, bench_expr, bench_meta, opts,
)


# 3.2 Check prereq
design_tibble <- design_tibble %>%
  mutate(net_bln = net_loc %>% check_prereq(),
         expr_bln = bnch_expr %>% check_prereq(),
         meta_bln = bench_meta %>% check_prereq())
design_tibble


# 3.3 Test with redesigned design tibble ----
bench_test <- design_tibble %>%
  mutate(activity = pmap(., function(name, net_loc, regs,
                                     gene_source, target, statistics,
                                     bnch_expr, bench_meta, opts,
                                     net_bln, expr_bln, meta_bln){

    # Check conditions and load prerequisites
    if(!net_bln){
      .GlobalEnv$network <- readRDS(net_loc)
    }
    if(!expr_bln){
      .GlobalEnv$gene_expression <- readRDS(bnch_expr) %>% as.matrix()
    }
    if(!meta_bln){
      .GlobalEnv$meta_data <- readRDS(bench_meta)
    }

    # filter network
    network_filtered <- network %>%
      dplyr::filter(confidence %in% regs) %>%
      distinct()

    # Obtain Activity with decouple and format
    decouple(mat = gene_expression, network = network_filtered,
             .source = all_of(gene_source), .target = target,
             statistics = statistics,
             .options = opts) %>%
      dplyr::rename(id=condition) %>%
      inner_join(meta_data, by="id") %>%
      dplyr::select(-c(.data$run_id, .data$p_value)) %>%
      group_split(statistic, .keep=T) %>%
      as.list()
  }))
bench_test

# format bench_test
bench_test2 <- bench_test %>%
  rowwise() %>%
  dplyr::mutate(statistics =
                  list(flatten_chr(.$activity[[1]] %>%
                                     map(function(tib)
                                       unique((tib)$statistic))))) %>%
  unnest(c(activity, statistics)) %>%
  rowwise() %>%
  mutate(name = paste(name, statistics,
                      sep = "_", (paste(unlist(regs), collapse = "")))) %>%
  ungroup() %>%
  mutate(roc = activity %>% map(calc_roc_curve)) %>%
  select(name, activity, roc)
bench_test2



# Get roc results
library(ggplot2)
bench_dt_plots <- plot_roc_auroc(bench_test2,
                                 title = "",
                                 coverage = TRUE)



# 3.4 Then save object with results and summary ----
out <- list(bench_test2, bench_dt_plots)
names(out) <- c("benchmark", "plots")
out$benchmark
out$plots$auroc_plot

saveRDS(out, here("inst/testdata/outputs", "regnet_cbd_results.rds"))
