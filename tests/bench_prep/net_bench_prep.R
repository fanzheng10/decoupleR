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

