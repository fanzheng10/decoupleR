# 0. Load Prerequisites ----
library(magrittr)
library(tidyverse)

bench_input <- system.file("benchdata", "inputs", package = "decoupleR")
bench_output <- system.file("benchdata", "outputs", package = "decoupleR")

# Load Real TFs from Lambert
tfs_real <- read_delim(file.path(bench_input,
                                 "TF_annotations_Lambert_2018.csv"),
                       delim = ";") %>%
  select(tf = Name, real = "Is TF?") %>%
  mutate(real = ifelse(real == "Yes", TRUE, FALSE)) %>%
  filter(real == TRUE) %>%
  pull(tf)


# 1. Format libraries in a standardized manner ----
# 1.1 Filter Dorothea ====
dorothea_filtered <- dorothea::dorothea_hs %>%
  mutate(likelihood = 1) %>%
  filter((tf %in% tfs_real))

# Save
saveRDS(dorothea_filtered, file.path(bench_input,
                                     "networks/dorothea_filtered.rds"))

# 1.2. Filter RegNetwork ====
high_confidence = read.csv(file.path(bench_input, "raw_networks/regnetwork",
                                "high_confidence.csv"),
                           header = T)
medium_confidence = read.csv(file.path(bench_input, "raw_networks/regnetwork/",
                                       "medium_confidence.csv"),
                             header = T)
low_confidence = read.csv(file.path(bench_input, "raw_networks/regnetwork",
                                    "low_confidence.csv"),
                          header = T)

# bind and format
regnetwork = rbind(high_confidence, medium_confidence, low_confidence) %>%
  filter(!str_detect(regulator_symbol, "miR|mir|hsa")) %>%
  filter(!str_detect(target_symbol, "miR|mir|hsa")) %>%
  rename(tf = regulator_symbol, target = target_symbol) %>%
  mutate(likelihood = 1, mor = 1) %>%
  select(tf, target, confidence,  likelihood, mor)
head(regnetwork)

saveRDS(regnetwork, file.path(bench_input, "raw_networks",
                              "regnetwork_raw.rds"))


regnetwork_raw <- readRDS(file.path(bench_input, "raw_networks",
                               "regnetwork_raw.rds"))

regnetwork_filtered <- regnetwork_raw %>%
  filter((tf %in% tfs_real))

saveRDS(regnetwork_filtered, file.path(bench_input, "networks",
                                  "regnetwork_filtered.rds"))



# chEA3
# 1.3. Filter chEA3 libs ====
load_chea3 <- function(confidence_list, source = "local") {

  if (source == "local"){
    library(here)
    # Load local copy of libraries
    paths <- c(file.path(bench_input, "raw_networks/chea3", "ARCHS4_Coexpression.gmt"),
               file.path(bench_input, "raw_networks/chea3", "ENCODE_ChIP-seq.gmt"),
               file.path(bench_input, "raw_networks/chea3", "Enrichr_Queries.gmt"),
               file.path(bench_input, "raw_networks/chea3", "GTEx_Coexpression.gmt"),
               file.path(bench_input, "raw_networks/chea3", "Literature_ChIP-seq.gmt"),
               file.path(bench_input, "raw_networks/chea3", "ReMap_ChIP-seq.gmt"))
  } else if (source == "web"){
    # Load libraries from ChEA3 website
    paths <- c("https://amp.pharm.mssm.edu/chea3/assets/tflibs/ARCHS4_Coexpression.gmt",
               "https://amp.pharm.mssm.edu/chea3/assets/tflibs/ENCODE_ChIP-seq.gmt",
               "https://amp.pharm.mssm.edu/chea3/assets/tflibs/Enrichr_Queries.gmt",
               "https://amp.pharm.mssm.edu/chea3/assets/tflibs/GTEx_Coexpression.gmt",
               "https://amp.pharm.mssm.edu/chea3/assets/tflibs/Literature_ChIP-seq.gmt",
               "https://amp.pharm.mssm.edu/chea3/assets/tflibs/ReMap_ChIP-seq.gmt")
  }

  all_chea3 <-  c("archs4_coexpression", "encode_chip_seq", "enrichr_queries",
                  "gtex_coexpression", "literature_chip_seq", "remap_chip_seq")

  names(paths) <- all_chea3

  if (confidence_list == "all") confidence_list <- all_chea3

  regulons <- map_dfr(paths[confidence_list], function(path) {

    # load file in gmt format (fill unequal row lengths with NA)
    num_col <-max(count.fields(path, sep = "\t"))
    col_names <- paste("X",1:num_col,sep="")
    df <- read.table(path, fill = TRUE, quote = "", col.names = col_names) %>%
      as_tibble()

    # extract library name from file path
    lib <- path %>%
      basename() %>%
      str_remove(".gmt") %>%
      str_to_lower() %>%
      str_replace("-", "_")

    # tidy gmt file
    tidy_df <- df %>%
      rename(tf = X1) %>%
      separate(tf, into = c("tf"), sep="_", extra = "drop") %>%
      group_by(tf) %>%
      mutate(experiment = paste0(tf,"_", row_number())) %>%
      ungroup() %>%
      select(-tf) %>%
      pivot_longer(-experiment, names_to = "key", values_to = "target") %>%
      mutate(tf = str_remove(experiment, "_.*"),
             experiment = str_remove(experiment, ".*_")) %>%
      select(tf, experiment, target,-key) %>%
      mutate(confidence = lib, likelihood = 1, mor = 1) %>%
      na_if("") %>%
      drop_na() %>%
      type_convert()
  })
}


# get chEA3 network
chea3_network <- load_chea3("all")  %>%
  select(-experiment)
saveRDS(chea3_network, file.path(bench_input, "raw_networks", "chea3_raw.rds"))


# (chea3_raw was generated from gmts using a decoupler-benchmark function)
chea3_raw <- readRDS(file.path(bench_input, "raw_networks" ,"chea3_raw.rds"))

# Filter
chea3_filtered <- chea3_raw %>%
  filter((tf %in% tfs_real))

# Save
saveRDS(chea3_filtered, file.path(bench_input, "networks",
                             "chea3_filtered.rds"))
file.path(bench_input, "dorothea_benchmark_data.rds")


# 2. Format Benchmark data sets -----
# 2.1. Filter DBD data ====
expr <- readRDS(file.path(bench_input, "raw_benchdata",
                          "dorothea_benchmark_data.rds")) %>%
  filter(organism == "human")

# Spread data (one column per experiment id)
gene_expression <- expr %>%
  select(gene, id, expression = t) %>%
  spread(id, expression, fill=0) %>%
  na.omit() %>%
  data.frame(row.names=1, stringsAsFactors = F, check.names = F)
head(gene_expression)
saveRDS(gene_expression, file.path(bench_input, "dorothea_bench_expr.rds"))

# save example data
saveRDS(gene_expression[,1:10], file.path(bench_input, "dorothea_bench_example.rds"))

# Save meta-data in separate data frame
meta_expr <- expr %>%
  dplyr::select(one_of("id", "accession", "tf", "platform", "info", "effect",
                       "source","sign")) %>%
  distinct() %>%
  rename(target = tf)
meta_expr
saveRDS(meta_expr, file.path(bench_input,"dorothea_bench_meta.rds"))


# 2.2. a) Filter KTF data  ====
expr = read.table(file.path(bench_input, "raw_benchdata",
                            "knockTF.txt"), sep = "\t", header = TRUE) %>%
  select("Sample_ID", "TF", "Gene", "Log2FC") %>%
  select(expression = Log2FC, id = Sample_ID,
         gene = Gene, tf = TF, everything()) %>%
  ungroup()

# spread data (one column per experiment id)
gene_expression = expr %>%
  select(gene, id, expression) %>%
  spread(id, expression, fill = 0) %>%
  na.omit() %>%
  data.frame(row.names=1, stringsAsFactors = F, check.names = F)
saveRDS(gene_expression, file.path(bench_input, "KnockTF_gene_expr.rds"))

# Save Example
saveRDS(gene_expression[,1:10], file.path(bench_input, "KnockTF_gene_example.rds"))

# save meta-data in separate data frame
meta_expr = expr %>%
  select(one_of("id", "tf")) %>%
  distinct() %>%
  rename(target = tf) %>%
  mutate(sign = -1) #knock-down or knock-out of all TFs
saveRDS(meta_expr, file.path(bench_input, "KnockTF_meta.rds"))

# clear env.
rm(list = ls())
bench_input <- system.file("benchdata", "inputs", package = "decoupleR")

# 2.2. b) KTF p-val as activity ====
expr = read.table(file.path(bench_input, "raw_benchdata",
                            "knockTF.txt"),
                  sep = "\t", header = TRUE)

expr_pvalue <- expr  %>%
  select("Sample_ID", "TF", "Gene", "Log2FC", "P_value") %>%
  select(id = Sample_ID, gene = Gene, tf = TF, everything()) %>%
  mutate(expression = abs(log2(as.double(P_value))) * sign(Log2FC)) %>%
  ungroup() %>%
  drop_na() # gene obsvs dropped in 70 experiments

# spread data (one column per experiment id)
gene_expression = expr_pvalue %>%
  select(gene, id, expression) %>%
  spread(id, expression, fill = 0) %>%
  na.omit() %>%
  data.frame(row.names=1, stringsAsFactors = F, check.names = F)
saveRDS(gene_expression, file.path(bench_input,
                              "KnockTF_pvalue_expr.rds"))

# Save Example
saveRDS(gene_expression[,1:10], file.path(bench_input,
                                     "KnockTF_pval_example.rds"))


# 2.2. c) KTF rank as activity ====
expr_rank <- expr %>%
  select("Sample_ID", "TF", "Gene", "Log2FC", "Rank") %>%
  select(id = Sample_ID, gene = Gene, tf = TF, everything()) %>%
  group_by(id) %>%
  add_count() %>%
  mutate(expression = abs(log2(Rank/n)) * sign(Log2FC)) %>%
  ungroup()

# spread data (one column per experiment id)
gene_expression = expr_rank %>%
  select(gene, id, expression) %>%
  spread(id, expression, fill = 0) %>%
  na.omit() %>%
  data.frame(row.names=1, stringsAsFactors = F, check.names = F)
saveRDS(gene_expression, file.path(bench_input, "KnockTF_rank_expr.rds"))

# Save Example
saveRDS(gene_expression[,1:10], file.path(bench_input,
                                     "KnockTF_rank_example.rds"))

# 3. Create Design Tibbles ----
# locations
bench_input <- system.file("benchdata", "inputs", package = "decoupleR")
bench_output <- system.file("benchdata", "outputs", package = "decoupleR")
# benchmark datasets
knock_expr <- file.path(bench_input,"KnockTF_gene_expr.rds")
knock_meta <- file.path(bench_input, "KnockTF_meta.rds")

dbd_expr <- file.path(bench_input, "dorothea_bench_expr.rds")
dbd_meta <- file.path(bench_input, "dorothea_bench_meta.rds")

# network locations
dorothea_loc <- file.path(bench_input, "networks/dorothea_filtered.rds")
chea3_loc <- file.path(bench_input, "networks/chea3_filtered.rds")
regnet_loc  <- file.path(bench_input, "networks/regnetwork_filtered.rds")

# available stats
statistics <- c(
  "mean",
  "pscira",
  "scira",
  "viper",
  "gsva" #,
  # "fgsea"
)

# options
opts <- list(
  scira = list(),
  pscira = list(),
  mean = list(),
  viper = list(verbose = FALSE, minsize=0),
  gsva = list(verbose = FALSE)
  # fgsea = list(options = list())
)


# 3.1. Dorothea + Dorothea Benchmark Data (DBD) ====
regs <- list(c("A"),
          c("A", "B", "C"),
          c("A", "B", "C", "D", "E"))

design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "dorothea_dbd", dorothea_loc, regs[[1]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "dorothea_dbd", dorothea_loc, regs[[2]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "dorothea_dbd", dorothea_loc, regs[[3]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  )

saveRDS(design_tibble, file.path(bench_input,
                            "designs/dorothea_dbd_design.rds"))


# 3.2. Dorothea + Knock_TF Data (KTF) ====
design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "dorothea_ktf", dorothea_loc, regs[[1]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "dorothea_ktf", dorothea_loc, regs[[2]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "dorothea_ktf", dorothea_loc, regs[[3]], "tf", "target", statistics, knock_expr, knock_meta, opts,
)

saveRDS(design_tibble, file.path(bench_input,
                            "designs/dorothea_ktf_design.rds"))


# 3.3. ChEA3 + DBD ====
ch_regs <- list("archs4_coexpression", "encode_chip_seq", "enrichr_queries",
          "gtex_coexpression", "literature_chip_seq", "remap_chip_seq")

design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "ChEA3_dbd", chea3_loc, ch_regs[[1]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "ChEA3_dbd", chea3_loc, ch_regs[[2]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "ChEA3_dbd", chea3_loc, ch_regs[[3]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "ChEA3_dbd", chea3_loc, ch_regs[[4]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "ChEA3_dbd", chea3_loc, ch_regs[[5]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "ChEA3_dbd", chea3_loc, ch_regs[[6]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
)
design_tibble$regs %<>% as.list

saveRDS(design_tibble, file.path(bench_input, "designs",
                            "ChEA3_dbd_design.rds"))

# 3.4. ChEA3 + KTF ====
design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "ChEA3_ktf", chea3_loc, ch_regs[[1]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "ChEA3_ktf", chea3_loc, ch_regs[[2]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "ChEA3_ktf", chea3_loc, ch_regs[[3]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "ChEA3_ktf", chea3_loc, ch_regs[[4]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "ChEA3_ktf", chea3_loc, ch_regs[[5]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "ChEA3_ktf", chea3_loc, ch_regs[[6]], "tf", "target", statistics, knock_expr, knock_meta, opts,
)
design_tibble$regs %<>% as.list

saveRDS(design_tibble, file.path(bench_input, "designs",
                            "ChEA3_ktf_design.rds"))


# 3.5. RegNetwork + DBD ====
rn_regs <- list(c("High"),
                  c("High", "Medium"),
                  c("High", "Medium", "Low"))


design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "regnet_dbd", regnet_loc, rn_regs[[1]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "regnet_dbd", regnet_loc, rn_regs[[2]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
  "regnet_dbd", regnet_loc, rn_regs[[3]], "tf", "target", statistics, dbd_expr, dbd_meta, opts,
)

saveRDS(design_tibble, file.path(bench_input,
                            "designs/regnet_dbd_design.rds"))

# 3.6. RegNetwork + KTF ====
design_tibble = tribble(
  ~name, ~net_loc, ~regs, ~gene_source, ~target, ~statistics, ~bnch_expr, ~bench_meta, ~opts,
  "regnet_ktf", regnet_loc,  rn_regs[[1]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "regnet_ktf", regnet_loc,   rn_regs[[2]], "tf", "target", statistics, knock_expr, knock_meta, opts,
  "regnet_ktf", regnet_loc,   rn_regs[[3]], "tf", "target", statistics, knock_expr, knock_meta, opts,
)

saveRDS(design_tibble, file.path(bench_input,
                            "designs/regnet_ktf_design.rds"))

# 3.7. Combine all designs ====
designs <-
  list.files(file.path(bench_input, "designs")) %>%
  map(function(file_loc)
    readRDS(file.path(bench_input, "designs", file_loc))) %>%
  bind_rows() %>%
  rename(row_name = name, # change names for clarity
         lvls = regs)

saveRDS(designs, file.path(bench_input, "design_all.rds"))


# 4. Generate Random Networks ----
# 4.1. Random regulons for Dorothea (Resampled) ====
dorothea_filtered <- readRDS(file.path(bench_input, "networks",
                                       "dorothea_filtered.rds"))


# get unique target gene pool
gene_pool <- dorothea_filtered %>%
  pull(target) %>%
  unique()

# Get Regulon Example
ADNP_regulon <- dorothea_filtered %>%
  group_by(tf) %>%
  add_count() %>%
  filter(tf == "ADNP")

# Get genes not part of the regulon (not target genes)
nt_gene_pool <- subset(gene_pool, !(gene_pool %in% ADNP_regulon$target))
summary(!(gene_pool %in% ADNP_regulon$target))
regulon_length <- ADNP_regulon$n


# Sample the same number of genes from random non-target genes
nt_gene_pool[(sample.int(length(nt_gene_pool), regulon_length))]


# Get random network
random_network <- dorothea_filtered %>%
  group_by(tf) %>%
  add_count() %>%
  select(tf, n, target) %>%
  nest(real_targets = target) %>%
  distinct() %>%
  rowwise() %>%
  mutate(available_genes = real_targets %>%
           map(function(tg) gene_pool[!(gene_pool %in% tg)])) %>%
  ungroup() %>%
  mutate(target = {
    map2(n, available_genes,
         .f = function(.x, .y){
           set.seed(69420)
           .y[(sample.int(length(.y), .x))]
         })
  }) %>%
  select(tf, target) %>%
  unnest(target)


# check if reproducible (xd is second run)
# all(xd==random_ADNP_regulon) # TRUE

random_ADNP_regulon <- random_network %>%
group_by(tf) %>%
  add_count() %>%
  filter(tf == "ADNP")
# check if appoints regulon with non-intersect gene of the same size
setdiff(ADNP_regulon$target, random_ADNP_regulon$target)
length(ADNP_regulon$target)
length(random_ADNP_regulon$target)

# assign remaining columns
random_network <- arrange(random_network)
random_network$confidence <- dorothea_filtered$confidence
random_network$likelihood <- dorothea_filtered$mor
random_network$mor <- dorothea_filtered$mor

# save network
saveRDS(random_network, file.path(bench_input, "networks",
                                  "dorothea_random.rds"))



# 4.2. Random regulons for Dorothea with BiRewire ----
library(BiRewire)
dorothea_net <- readRDS(file.path(bench_input,
                                  "networks","dorothea_filtered.rds"))

dorothea_bi <- dorothea_net %>%
  select(tf, mor, target)


# Induced bipartite and SIF builder
dorothea_dsg = birewire.induced.bipartite(dorothea_bi,
                                          delimitators = list(negative = '-1', positive = '1'))
dorothea_sif = birewire.build.dsg(dorothea_dsg,
                                  delimitators = list(negative = '-1', positive = '1'))

# Rewire dsg
random_dsg = birewire.rewire.dsg(dsg = dorothea_dsg)
random_sif = birewire.build.dsg(random_dsg,
                                delimitators = list(negative = '-1', positive = '1'))

# Jacard dsg
birewire.similarity.dsg(dorothea_dsg, random_dsg)


# Check Results
dorothea_sif <- dorothea_sif %>%
  rename(tf = source,
         mor = sign) %>%
  as_tibble()

ADNP_reg <- dorothea_sif %>%
  filter(tf == "ADNP")

ADNP_reg2 <- dorothea_net %>%
  filter(tf == "ADNP")

ADNP_reg3 <- random_sif %>%
  rename(tf = source,
         mor = sign) %>%
  filter(tf == "ADNP")

setdiff(ADNP_reg$target, ADNP_reg2$target)
setdiff(ADNP_reg$target, ADNP_reg3$target)

# remodel random network as Dorothea
random_net <- random_sif %>%
  rename(tf = source,
         mor = sign) %>%
  as_tibble() %>%
  arrange(tf)

# join confidence and likelihood
all(dorothea_net$tf == random_net$tf)

random_net <- random_net %>%
  mutate(likelihood=1,
         confidence=dorothea_net$confidence,
         mor=as.numeric(mor)) %>%
  select(tf, confidence, target, mor, likelihood)


# there seems to be duplicated TF-target edges in the BiRewire network
random_net <- random_net %>%
  slice(which(!duplicated(random_net %>% select(tf, target))))

saveRDS(random_net, file.path(bench_input,
                              "networks", "dorothea_birewire.rds"))



# 5. Kinases Test ----
# Dorothea BD
dbd_ge <- readRDS(file.path(bench_input, "dorothea_bench_example.rds"))
dbd_me <- readRDS(file.path(bench_input,"dorothea_bench_meta.rds"))


# Load kbd
load(file.path(bench_input, "KinaseData", "KSets.rdata"))
load(file.path(bench_input, "KinaseData", "phosphosite_KinaseTargets.rdata"))
kinase_network <- KSets
rm(KSets)

# Reform Kinase Bench data
kinase_bd <- read_csv(file.path(bench_input, "KinaseData", "Kinase_BM_data.csv")) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("X1")
saveRDS(kinase_bd, file.path(bench_input, "kprep", "kinase_bd.rds"))



# Reform Kinase Meta data
kinase_meta <- read.csv(file.path(bench_input, "KinaseData", "KinaseConditionPairs.csv")) %>%
  rename(id="Condition",
         target="Kinase") %>%
  mutate(sign = ifelse(Regulation=="up", 1, -1))

saveRDS(kinase_meta, file.path(bench_input, "kprep", "kinase_meta.rds"))

# Reform Network
readRDS(file.path(bench_input,
                  "networks", "dorothea_birewire.rds"))


kn <- kinase_network %>%
  select(Regulator, hgnc_symbol, Evidence) %>%
  mutate(mor = 1, likelihood = 1)
saveRDS(kn, file.path(bench_input, "kprep", "kinase_network.rds"))
