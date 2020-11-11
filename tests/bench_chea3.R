library(here)
# Load Prerequisites ----
# Locations for the network resource, benchmark gene expression, and meta data
net_loc <- here("inst/testdata/inputs/ChEA3", "ChEA3_realTFs_only.RDS")
# bench_expr <- here("inst/testdata/inputs", "knockTF_ex.rds")
# bench_meta <- here("inst/testdata/inputs", "knockTF_meta.rds")
bench_expr <- here("inst/testdata/inputs", "input-dorothea_bench_expr.rds")
bench_meta <- here("inst/testdata/inputs", "input-dorothea_bench_meta.rds")

x <- readRDS(here("inst/testdata/inputs", "input-dorothea_bench_expr.rds"))
m <- readRDS(here("inst/testdata/inputs", "input-dorothea_bench_meta.rds"))
#
# x2 <- readRDS(here("inst/testdata/inputs", "knockTF_ex.rds"))
# m2 <- readRDS(here("inst/testdata/inputs", "knockTF_meta.rds"))

# head(x)
# head(x2)
#
# head(m)
# head(m2)

# Libraries/Confidence levels in each gene set (or anything that )
regs <- c("archs4_coexpression", "encode_chip_seq", "enrichr_queries",
          "gtex_coexpression", "literature_chip_seq", "remap_chip_seq" )



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


design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "chEA3", net_loc, regs[1], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "chEA3", net_loc, regs[2], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "chEA3", net_loc, regs[3], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "chEA3", net_loc, regs[4], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "chEA3", net_loc, regs[5], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "chEA3", net_loc, regs[6], "tf", "target", statistics, bench_expr, bench_meta, opts,
)
design_tibble[1,3]

# #  Load bench and network data
# gene_expression <- readRDS(here("inst/testdata/inputs/", "input-dorothea_bench_example.rds")) %>%
#   as.matrix()
# meta_data <- readRDS(here("inst/testdata/inputs/", "input-dorothea_bench_meta.rds"))
# dorothea_genesets <- readRDS(here("inst/testdata/inputs/", "input_dorothea_full.rds"))



# 3.2 Function that checks preceding vector element ----
check_prereq <- function(loc_vector){
  loc_tib <- tibble(current=loc_vector, behind=lag(loc_vector))

  pmap_lgl(loc_tib, function(behind, current){
    if(is.na(behind) || behind!=current){
      FALSE
    } else{
      TRUE
    }
  })
}

# Check prereq
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
