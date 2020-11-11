
# 1. Convert Dorothea and KnockTF Bnch ----
# ** Dorothea BNCH ====
expr <- readRDS(here("inst/testdata/inputs", "dorothea_benchmark_data.rds")) %>%
  filter(organism == "human")
head(expr)


# Spread data (one column per experiment id)
gene_expression <- expr %>%
  select(gene, id, expression = t) %>%
  spread(id, expression, fill=0) %>%
  na.omit() %>%
  data.frame(row.names=1, stringsAsFactors = F, check.names = F)
head(gene_expression)
saveRDS(gene_expression, here("inst/testdata/inputs/", "input-dorothea_bench_expr.rds"))

# save example data
saveRDS(gene_expression[,1:10], here("inst/testdata/inputs/", "input-dorothea_bench_example.rds"))

# Save meta-data in separate data frame
meta_expr <- expr %>%
  dplyr::select(one_of("id", "accession", "tf", "platform", "info", "effect",
                       "source","sign")) %>%
  distinct() %>%
  rename(target = tf)
meta_expr
saveRDS(meta_expr, here("inst/testdata/inputs/", "input-dorothea_bench_meta.rds"))

# ** KnockTF ====
expr = read.table(here("inst/testdata/inputs", "knockTF.txt"), sep = "\t", header = TRUE)  %>%
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
saveRDS(gene_expression, here("inst/testdata/inputs/", "KnockTF_gene_expr.rds"))

# Save Example
saveRDS(gene_expression[,1:10], here("inst/testdata/inputs/", "KnockTF_gene_example.rds"))


# save meta-data in separate data frame
meta_expr = expr %>%
  select(one_of("id", "tf")) %>%
  distinct() %>%
  rename(target = tf) %>%
  mutate(sign = -1) #knock-down or knock-out of all TFs
saveRDS(meta_expr, here("inst/testdata/inputs/", "KnockTF_meta_data.rds"))



# 2. Format and Save libs  ----
# 2.1 ChEA3 ----
load_chea3 <- function(confidence_list, source = "local") {

  if (source == "local"){
    # Load local copy of libraries
    paths <- c(here("inst/testdata/inputs", "ChEA3", "ARCHS4_Coexpression.gmt"),
               here("inst/testdata/inputs", "ChEA3", "ENCODE_ChIP-seq.gmt"),
               here("inst/testdata/inputs", "ChEA3", "Enrichr_Queries.gmt"),
               here("inst/testdata/inputs", "ChEA3", "GTEx_Coexpression.gmt"),
               here("inst/testdata/inputs", "ChEA3", "Literature_ChIP-seq.gmt"),
               here("inst/testdata/inputs", "ChEA3", "ReMap_ChIP-seq.gmt"))
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


# List of real TFs
tf_list <- read_delim(here("inst/testdata/inputs/",
                           "TF_annotations_Lambert_2018.csv"),
                      delim = ";") %>%
  select(tf = Name, real = "Is TF?") %>%
  mutate(real = ifelse(real == "Yes", TRUE, FALSE))
head(tf_list)

# get chEA3 network
chea3_network <- load_chea3("all")
head(chea3_network)

# Filter putative TFs
tfs_true <- tf_list %>% filter(real == TRUE) %>% pull(tf)
length(tfs_true)

chea3_network_filtered <- chea3_network %>%
  filter(tf %in% tfs_true) %>%
  select(-experiment)
head(chea3_network_filtered)

# Save ChEA3 network
saveRDS(chea3_network_filtered, here("inst/testdata/inputs/ChEA3",
                                     "ChEA3_realTFs_only.RDS"))
