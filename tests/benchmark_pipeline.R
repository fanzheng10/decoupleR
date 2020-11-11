# 0. Load data ----
library(here)
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

dorothea_genesets <- dorothea::dorothea_hs %>%
  # head(20000) %>%
  mutate(likelihood = 1)
saveRDS(dorothea_genesets, here("inst/testdata/inputs/", "input_dorothea_full.rds"))


mat <- file.path(input_dir, "input-expr_matrix.rds") %>%
  readRDS()

knock_ex <- here("inst/testdata/inputs/", "knockTF_ex.rds") %>%
  readRDS() %>%
  as.matrix()

knock_full <-here("inst/testdata/inputs/knockTF_gene_expr.rds") %>%
  readRDS() %>%
  as.matrix()

meta_expr <- here("inst/testdata/inputs/", "knockTF_meta.rds") %>%
  readRDS()


# *. Run decouple with knock_tf example ----

# Available statistics
statistics <- c(
  "scira",
  "pscira",
  "mean",
  "viper",
  "gsva"
)

# decouple_knock_ex <- decouple(
#   mat = knock_ex,
#   network = dorothea_genesets,
#   .source = tf,
#   .target = target,
#   statistics = statistics,
#   .options = list(
#     scira = list(),
#     pscira = list(),
#     mean = list(.likelihood = NULL),
#     viper = list(options = list(verbose = FALSE)),
#     gsva = list(options = list(verbose = FALSE))
#   )
# ) %>%
#   dplyr::select(-.data$run_id) %>%
#   dplyr::arrange(.data$statistic, .data$tf, .data$condition)



# 1. Prepare Benchmark pipeline steps ----
# 1.1 Preps for benchmark pipeline ----
design_tibble = tribble(
  ~net_name, ~confidence_levels, ~gene_source, ~target, ~statistics,
  "dorothea", c("A"), "tf", "target", statistics,
  "dorothea", c("A", "B"), "tf", "target", statistics,
  # "dorothea", c("A", "B", "C"), "tf", "target", statistics,
)
design_tibble

bench_test <- design_tibble %>%
  mutate(activity = pmap(., function(net_name, confidence_levels,
                                        gene_source, target, statistics){
    # filter network
    network <- dorothea_genesets %>%
      filter(confidence %in% confidence_levels)

    # calculate activity
    decouple(mat = knock_full, network = network,
             .source = all_of(gene_source), .target = target,
             statistics = statistics,
             .options = list(
               scira = list(),
               pscira = list(),
               mean = list(.likelihood = NULL),
               viper = list(options = list(verbose = FALSE, minsize=0)),
               gsva = list(options = list(verbose = FALSE)))) %>%
      dplyr::select(-c(.data$run_id, .data$p_value)) %>%
      dplyr::rename(id=condition)
  }))
bench_test

# 1.2 format activity ----
bench_test$activity[[1]]

activity_list.format <- bench_test$activity[[1]] %>%
  group_split(statistic) %>%
  set_names(levels(as.factor(bench_test$activity[[1]]$statistic))) %>%
  map(function(l) l[-1]) # remove statistic column
activity_list.format

# 1.3 join activity with meta data by experiment ----
activity_list.join <- activity_list.format %>%
  map(function(l) (l %>%
                     inner_join(., meta_expr, by="id") # %>%
                     # rename(experiment_id = id)
                   ))
activity_list.join


# 1.4 Combine 1.1 to 1.3----
bench_combined <- bench_test %>%
  mutate(statistics = list(levels(as.factor(activity[[1]]$statistic))))  %>%
  mutate(activity = activity %>%
           map(function(tib) tib %>%
                 inner_join(meta_expr, by="id") %>%
                 group_split(statistic, .keep=FALSE) %>%
                 as.list()
           ))  %>%
  unnest(c(activity, statistics)) %>%
  rowwise() %>%
  mutate(name = paste(net_name, statistics,
                          paste0(confidence_levels, collapse = ""),
                          sep = "_")) %>%
  ungroup()
bench_combined



# 1.5 Get ROC ----
bench_roc <- bench_combined %>%
  mutate(roc = activity %>% map(calc_roc_curve))
bench_roc


# 1.6 Auroc Plots ----
bench_plots <- plot_roc_auroc(bench_roc, title = "Knock_TF", coverage = TRUE)
bench_plots$coverage

# 2. Pipeline miscs  ----
# 2.1 Import Design tibble as JSON ----
library(jsonlite)
library(here)

design_tibble = tribble(
  ~net_name, ~confidence_levels, ~gene_source, ~target, ~statistics,
  "dorothea", c("A"), "tf", "target", statistics,
  "dorothea", c("A", "B"), "tf", "target", statistics,
  # "dorothea", c("A", "B", "C"), "tf", "target", statistics,
)

write(toJSON(design_tibble, pretty=T),
      here("inst/testdata/inputs/", "design_tibble.json"))

##

design_tibble_test <- fromJSON(here("inst/testdata/inputs/",
                                    "design_tibble.json"))
design_tibble_test

gene_expression <- readRDS(here("inst/testdata/inputs/", "knockTF_ex.rds"))
meta_data <- readRDS(here("inst/testdata/inputs/", "knockTF_meta.rds"))
dorothea_genesets <- dorothea::dorothea_hs %>%
  # head(20000) %>%
  mutate(likelihood = 1)


bench_dt_test <- design_tibble %>%
  mutate(activity = pmap(., function(net_name, confidence_levels,
                                     gene_source, target, statistics){
    # filter network
    network <- dorothea_genesets %>%
      filter(confidence %in% confidence_levels)

    # calculate activity
    decouple(mat = mat, network = network,
             .source = all_of(gene_source), .target = target,
             statistics = statistics,
             .options = list(
               scira = list(),
               pscira = list(),
               mean = list(.likelihood = NULL),
               viper = list(options = list(verbose = FALSE, minsize=0)),
               gsva = list(options = list(verbose = FALSE)))) %>%
      dplyr::select(-c(.data$run_id, .data$p_value)) %>%
      dplyr::rename(id=condition)
  }))

bench_dt_test <- bench_dt_test  %>%
  mutate(statistics = list(levels(as.factor(activity[[1]]$statistic))))  %>%
  mutate(activity = activity %>%
           map(function(tib) tib %>%
                 inner_join(meta_expr, by="id") %>%
                 group_split(statistic, .keep=FALSE) %>%
                 as.list()
           ))  %>%
  unnest(c(activity, statistics)) %>%
  rowwise() %>%
  mutate(name = paste(net_name, statistics,
                      paste0(confidence_levels, collapse = ""),
                      sep = "_")) %>%
  ungroup()

bench_dt_test <- bench_dt_test %>%
  mutate(roc = activity %>% map(calc_roc_curve))


bench_dt_plots <- plot_roc_auroc(bench_roc, title = "Knock_TF", coverage = TRUE)
bench_dt_plots$coverage
bench_dt_plots$roc_plot
bench_dt_plots$auroc_plot


# 2.2 Redesign JSON ----
library(here)

# Locations for the network resource, benchmark gene expression, and meta data
net_loc <- here("inst/testdata/inputs", "input_dorothea_full.rds")
bench_expr <- here("inst/testdata/inputs", "knockTF_ex.rds")
bench_meta <- here("inst/testdata/inputs", "knockTF_meta.rds")


# Libraries/Confidence levels in each gene set (or anything that )
regs <- c("A", "B", "C", "D", "E")


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
  viper = list(options = list(verbose = FALSE, minsize=0)),
  gsva = list(options = list(verbose = FALSE))
  )


design_tibble = tribble(
  ~net_name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "dorothea", net_loc, regs[1], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "dorothea", net_loc, regs[1:2], "tf", "target", statistics, bench_expr, bench_meta, opts,
  # "dorothea", net_loc, regs[1:3], "tf", "target", statistics, bench_expr, bench_meta, opts,
  # "dorothea", net_loc, regs, "tf", "target", statistics, bench_expr, bench_meta, opts,
)
design_tibble


# Test
gene_expression <- readRDS(here("inst/testdata/inputs/", "input-dorothea_bench_example.rds")) %>%
  as.matrix()
meta_data <- readRDS(here("inst/testdata/inputs/", "input-dorothea_bench_meta.rds"))
dorothea_genesets <- readRDS(here("inst/testdata/inputs/", "input_dorothea_full.rds"))

# 2.3 Test with redesigned deisgn tibble ----

bench_test <- design_tibble %>%
  mutate(activity = pmap(., function(net_name, net_loc, regs,
                                     gene_source, target, statistics,
                                     bnch_expr, bench_meta, opts){

    # Check conditions and load prerequisites
    network <- dorothea_genesets %>%
      dplyr::filter(confidence %in% regs)

    # Obtain Activity with decouple
    decouple(mat = gene_expression, network = network,
             .source = all_of(gene_source), .target = target,
             statistics = statistics,
             .options = list(
               scira = list(),
               pscira = list(),
               mean = list(.likelihood = NULL),
               viper = list(options = list(verbose = FALSE, minsize=0)),
               gsva = list(options = list(verbose = FALSE))
             )) %>%
      dplyr::select(-c(.data$run_id, .data$p_value)) %>%
      dplyr::rename(id=condition)
  }))


# format tibble and calculate ROC
bench_test2 <- bench_test %>%
  mutate(statistics = list(levels(as.factor(activity[[1]]$statistic))))  %>%
  mutate(activity = activity %>%
           map(function(tib) tib %>%
                 inner_join(meta_data, by="id") %>%
                 group_split(statistic, .keep=F) %>%
                 as.list()
           )) %>%
  unnest(c(activity, statistics)) %>%
  rowwise() %>%
  mutate(name = paste(net_name, statistics,
                      paste0(regs, collapse = ""),
                      sep = "_")) %>%
  ungroup() %>%
  mutate(roc = activity %>% map(calc_roc_curve))



# Get roc results
library(ggplot2)
bench_dt_plots <- plot_roc_auroc(bench_test2, title = "", coverage = TRUE)
bench_dt_plots$auroc_plot


bench_dt_plots$coverage
bench_dt_plots$roc_plot
bench_dt_plots$auroc_plot

# Then save object with results and summary


# 2.4 Function that checks preceding vector element ----
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

check_prereq(design_tibble$net_loc)



# 3. Assemble into a well-formatted pipeline ----
library(here)
# 3.1 Load Prerequisites ----
# Locations for the network resource, benchmark gene expression, and meta data
net_loc <- here("inst/testdata/inputs", "input_dorothea_full.rds")
bench_expr <- here("inst/testdata/inputs", "knockTF_ex.rds")
bench_meta <- here("inst/testdata/inputs", "knockTF_meta.rds")

bench_expr2 <- here("inst/testdata/inputs", "input-dorothea_bench_example.rds")
bench_meta2 <- here("inst/testdata/inputs", "input-dorothea_bench_meta.rds")


# Libraries/Confidence levels in each gene set (or anything that )
regs <- c("A", "B", "C", "D", "E")

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
  viper = list(options = list(verbose = FALSE, minsize=0)),
  gsva = list(options = list(verbose = FALSE))
)


design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "dorothea", net_loc, regs[1], "tf", "target", statistics, bench_expr, bench_meta, opts,
  "dorothea", net_loc, regs[1:3], "tf", "target", statistics, bench_expr, bench_meta, opts,
  # "dorothea", net_loc, regs[1:3], "tf", "target", statistics, bench_expr, bench_meta, opts,
  # "dorothea", net_loc, regs, "tf", "target", statistics, bench_expr, bench_meta, opts,
  "dorothea_cbd", net_loc, regs[1:3], "tf", "target", statistics, bench_expr2, bench_meta2, opts,
)
design_tibble

#  Load bench and network data
gene_expression <- readRDS(here("inst/testdata/inputs/", "input-dorothea_bench_example.rds")) %>%
  as.matrix()
meta_data <- readRDS(here("inst/testdata/inputs/", "input-dorothea_bench_meta.rds"))
dorothea_genesets <- readRDS(here("inst/testdata/inputs/", "input_dorothea_full.rds"))



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
      dplyr::filter(confidence %in% regs)

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
bench_test$activity[[1]]

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
out$plots
