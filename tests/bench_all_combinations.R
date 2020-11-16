# 1. Format libraries in a standardized manner ----
library(here)
library(tidyverse)
library(magrittr)

# Load Real TFs
tfs_real <- read_delim(here("inst/benchdata/inputs/",
                           "TF_annotations_Lambert_2018.csv"),
                      delim = ";") %>%
  select(tf = Name, real = "Is TF?") %>%
  mutate(real = ifelse(real == "Yes", TRUE, FALSE)) %>%
  filter(real == TRUE) %>%
  pull(tf)

# 1.1 Filter Dorothea ====
dorothea_raw <- dorothea::dorothea_hs %>%
  mutate(likelihood = 1)

# Filter
dorothea_filtered <- dorothea_raw %>%
  filter((tf %in% tfs_real))

# Save
saveRDS(dorothea_filtered, here("inst/benchdata/inputs/",
                                "networks/dorothea_filtered.rds"))

# 1.2. Filter chEA3 libs ====
# (chea3_raw was generated from gmts using a decoupler-benchmark function)
chea3_raw <- readRDS(here("inst/benchdata/inputs/",
                          "networks/chea3_raw.rds"))

# Filter
chea3_filtered <- chea3_raw %>%
  filter((tf %in% tfs_real))

# Save
saveRDS(chea3_filtered, here("inst/benchdata/inputs/",
                                "networks/chea3_filtered.rds"))


# 1.3. Filter regnetwork ====
regnetwork_raw <- readRDS(here("inst/benchdata/inputs/",
                               "networks/regnetwork_raw.rds"))

regnetwork_filtered <- regnetwork_raw %>%
  filter((tf %in% tfs_real))


saveRDS(regnetwork_filtered, here("inst/benchdata/inputs/",
                             "networks/regnetwork_filtered.rds"))

# 2. Format Benchmark data sets -----
# 2.1. Filter DBD data ====
expr <- readRDS(here("inst/testdata/inputs", "dorothea_benchmark_data.rds")) %>%
  filter(organism == "human")

# Spread data (one column per experiment id)
gene_expression <- expr %>%
  select(gene, id, expression = t) %>%
  spread(id, expression, fill=0) %>%
  na.omit() %>%
  data.frame(row.names=1, stringsAsFactors = F, check.names = F)
head(gene_expression)
saveRDS(gene_expression, here("inst/benchdata/inputs/", "input-dorothea_bench_expr.rds"))

# save example data
saveRDS(gene_expression[,1:10], here("inst/benchdata/inputs/", "input-dorothea_bench_example.rds"))

# Save meta-data in separate data frame
meta_expr <- expr %>%
  dplyr::select(one_of("id", "accession", "tf", "platform", "info", "effect",
                       "source","sign")) %>%
  distinct() %>%
  rename(target = tf)
meta_expr
saveRDS(meta_expr, here("inst/benchdata/inputs/", "input-dorothea_bench_meta.rds"))

# 2.2. Filter KTF data  ====
expr = read.table(here("inst/testdata/inputs", "knockTF.txt"), sep = "\t", header = TRUE) %>%
  select("Sample_ID", "TF", "Gene", "Log2FC") %>%
  select(expression = Log2FC, id = Sample_ID,
         gene = Gene, tf = TF, everything()) %>%
  group_by(id, gene, tf) %>%
  summarize(expression = mean(expression)) %>%
  ungroup()

# spread data (one column per experiment id)
gene_expression = expr %>%
  select(gene, id, expression) %>%
  spread(id, expression, fill = 0) %>%
  na.omit() %>%
  data.frame(row.names=1, stringsAsFactors = F, check.names = F)
saveRDS(gene_expression, here("inst/benchdata/inputs/", "KnockTF_gene_expr.rds"))

# Save Example
saveRDS(gene_expression[,1:10], here("inst/benchdata/inputs/", "KnockTF_gene_example.rds"))


# save meta-data in separate data frame
meta_expr = expr %>%
  select(one_of("id", "tf")) %>%
  distinct() %>%
  rename(target = tf) %>%
  mutate(sign = -1) #knock-down or knock-out of all TFs
saveRDS(meta_expr, here("inst/benchdata/inputs/", "KnockTF_meta.rds"))

# clear env.
rm(list = ls())


# 3. Create Design JSONs ----

# locations
# benchmark datasets
knock_expr <- here("inst/benchdata/inputs", "KnockTF_gene_expr.rds")
knock_meta <- here("inst/benchdata/inputs", "KnockTF_meta.rds")

dbd_expr <- here("inst/benchdata/inputs", "input-dorothea_bench_example.rds")
dbd_meta <- here("inst/benchdata/inputs", "input-dorothea_bench_meta.rds")

# network locations
dorothea_loc <- here("inst/benchdata/inputs/networks", "dorothea_filtered.rds")
chea3_loc <- here("inst/benchdata/inputs/networks", "chea3_filtered.rds")
regnet_loc  <- here("inst/benchdata/inputs/networks", "regnetwork_filtered.rds")

# available stats
statistics <- c(
  "mean",
  "pscira",
  "scira",
  "viper",
  "gsva",
  "fgsea"
)

# options
opts <- list(
  scira = list(),
  pscira = list(),
  mean = list(),
  viper = list(options = list(verbose = FALSE)),
  gsva = list(options = list(verbose = FALSE)),
  fgsea = list(options = list())
)


# 3.1. Dorothea + Dorothea Benchmark Data (DBD) ====
regs <- list(c("A"),
          c("A", "B", "C"),
          c("A", "B", "C", "D", "E"))

design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta,
  "dorothea_dbd", dorothea_loc, regs[[1]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "dorothea_dbd", dorothea_loc, regs[[2]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "dorothea_dbd", dorothea_loc, regs[[3]], "tf", "target", statistics, dbd_expr, dbd_meta,
  )

saveRDS(design_tibble, here("inst/benchdata/inputs/designs",
                            "dorothea_dbd_design.rds"))





# 3.2. Dorothea + Knock_TF Data (KTF) ====
design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta,
  "dorothea_ktf", dorothea_loc, regs[[1]], "tf", "target", statistics, knock_expr, knock_meta,
  "dorothea_ktf", dorothea_loc, regs[[2]], "tf", "target", statistics, knock_expr, knock_meta,
  "dorothea_ktf", dorothea_loc, regs[[3]], "tf", "target", statistics, knock_expr, knock_meta,
)

saveRDS(design_tibble, here("inst/benchdata/inputs/designs",
                            "dorothea_ktf_design.rds"))


# 3.3. ChEA3 + DBD ====
ch_regs <- list("archs4_coexpression", "encode_chip_seq", "enrichr_queries",
          "gtex_coexpression", "literature_chip_seq", "remap_chip_seq")

design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta,
  "ChEA3_dbd", chea3_loc, ch_regs[[1]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "ChEA3_dbd", chea3_loc, ch_regs[[2]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "ChEA3_dbd", chea3_loc, ch_regs[[3]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "ChEA3_dbd", chea3_loc, ch_regs[[4]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "ChEA3_dbd", chea3_loc, ch_regs[[5]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "ChEA3_dbd", chea3_loc, ch_regs[[6]], "tf", "target", statistics, dbd_expr, dbd_meta,
)
design_tibble$regs %<>% as.list

saveRDS(design_tibble, here("inst/benchdata/inputs/designs",
                            "ChEA3_dbd_design.rds"))

# 3.4. ChEA3 + KTF ====
design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta,
  "ChEA3_ktf", chea3_loc, ch_regs[[1]], "tf", "target", statistics, knock_expr, knock_meta,
  "ChEA3_ktf", chea3_loc, ch_regs[[2]], "tf", "target", statistics, knock_expr, knock_meta,
  "ChEA3_ktf", chea3_loc, ch_regs[[3]], "tf", "target", statistics, knock_expr, knock_meta,
  "ChEA3_ktf", chea3_loc, ch_regs[[4]], "tf", "target", statistics, knock_expr, knock_meta,
  "ChEA3_ktf", chea3_loc, ch_regs[[5]], "tf", "target", statistics, knock_expr, knock_meta,
  "ChEA3_ktf", chea3_loc, ch_regs[[6]], "tf", "target", statistics, knock_expr, knock_meta,
)
design_tibble$regs %<>% as.list

saveRDS(design_tibble, here("inst/benchdata/inputs/designs",
                            "ChEA3_ktf_design.rds"))


# 3.5. RegNetwork + DBD ====
rn_regs <- list(c("High"),
                  c("High", "Medium"),
                  c("High", "Medium", "Low"))


design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta,
  "regnet_dbd", regnet_loc, rn_regs[[1]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "regnet_dbd", regnet_loc, rn_regs[[2]], "tf", "target", statistics, dbd_expr, dbd_meta,
  "regnet_dbd", regnet_loc, rn_regs[[3]], "tf", "target", statistics, dbd_expr, dbd_meta,
)

saveRDS(design_tibble, here("inst/benchdata/inputs/designs",
                            "regnet_dbd_design.rds"))

# 3.6. RegNetwork + KTF ====
design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta,
  "regnet_ktf", regnet_loc,  rn_regs[[1]], "tf", "target", statistics, knock_expr, knock_meta,
  "regnet_ktf", regnet_loc,   rn_regs[[2]], "tf", "target", statistics, knock_expr, knock_meta,
  "regnet_ktf", regnet_loc,   rn_regs[[3]], "tf", "target", statistics, knock_expr, knock_meta,
)

saveRDS(design_tibble, here("inst/benchdata/inputs/designs",
                            "regnet_ktf_design.rds"))


# 3.7. Combine all designs ====
designs <-
  list.files(here("inst/benchdata/inputs/designs")) %>%
  map(function(file_loc)
    readRDS(here("inst/benchdata/inputs/designs", file_loc))) %>%
  bind_rows()

saveRDS(designs,
        here("inst/benchdata/inputs/",
             "all_combinations.rds"))

designs_ktf <- designs %>%
  filter(str_detect(name, "ktf"))

saveRDS(designs_ktf,
        here("inst/benchdata/inputs/",
             "all_ktf.rds"))


designs_dbd <- designs %>%
  filter(str_detect(name, "dbd"))


saveRDS(designs_dbd,
        here("inst/benchdata/inputs/",
             "all_dbd.rds"))


# 4. Benchmark ----
opts <- list(
  scira = list(),
  pscira = list(),
  mean = list(),
  viper = list(options = list(verbose = FALSE)),
  gsva = list(options = list(verbose = FALSE)),
  fgsea = list(options = list())
)



dbd_test_dor <- run_benchmark(here("inst/benchdata/inputs/designs",
                                   "dorothea_dbd_design.rds"), opts)



dbd_test_dor_format <- dbd_test_dor %>%
  rowwise() %>%
  dplyr::mutate(statistics =
                  list(flatten_chr(.$activity[[1]] %>%
                                     map(function(tib)
                                       unique((tib)$statistic))))) %>%
  unnest(c(activity, statistics)) %>%
  mutate(stime = activity %>% map(function(tib) unique(tib$stime))) %>%
  mutate(rtime = activity %>% map(function(tib) unique(tib$rtime))) %>%
  unnest(c(stime, rtime)) %>%
  arrange(stime) %>%
  select(name, statistics, regs, ctime, stime, rtime)


# ctime <- dbd_test_dor_format[1,]$ctime
# rtime1 <- dbd_test_dor_format[1,]$rtime
# rtime2 <- dbd_test_dor_format[14,]$rtime
#
# difftime(rtime1, ctime)
# difftime(rtime2, rtime1)


ctime_value <- unique(dbd_test_dor_format$ctime)
rtime_vector <- unique(dbd_test_dor_format$rtime)

tibble(current=rtime_vector, behind=lag(rtime_vector))

tib_loc <- tibble(current=vector_loc, behind=lag(vector_loc))

pmap_lgl(tib_loc, function(behind, current){
  ifelse(is.na(behind) || behind!=current, FALSE, TRUE)
})

# # 4.2 Run DBD ====
# library(here)
# dbd_all_combs <- run_benchmark(here("inst/benchdata/inputs/",
#                    "all_dbd.rds"), opts)
#
# saveRDS(dbd_all_combs, here("inst/benchdata/outputs/dbd_all_perf.rds"))
# # 4.3 Run Knock TF ====
# readRDS(here("inst/benchdata/inputs/",
#              "all_dbd.rds"))


# 5. Summarize and Visualize output ----

