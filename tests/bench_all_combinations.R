# 1. Format libraries in a standardized manner ----

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
  filter((tf %in% tfs_real)) %>%
  group_by(tf, confidence) %>%
  add_count() %>%
  filter(n >= 10)  %>%
  select(-n)

# Save
saveRDS(dorothea_filtered, here("inst/benchdata/inputs/",
                                "networks/dorothea_filtered.rds"))

# 1.2. Filter chEA3 libs ====
# (chea3_raw was generated from gmts using a decoupler-benchmark function)
chea3_raw <- readRDS(here("inst/benchdata/inputs/",
                          "networks/chea3_raw.rds"))

# Filter
chea3_filtered <- chea3_raw %>%
  filter((tf %in% tfs_real)) %>%
  group_by(tf, confidence) %>%
  add_count() %>%
  filter(n >= 10) %>%
  select(-n)

# Save
saveRDS(chea3_filtered, here("inst/benchdata/inputs/",
                                "networks/chea3_filtered.rds"))


# 1.3. Filter regnetwork ====
regnetwork_raw <- readRDS(here("inst/benchdata/inputs/",
                               "networks/regnetwork_raw.rds"))

regnetwork_filtered <- regnetwork_raw %>%
  filter((tf %in% tfs_real)) %>%
  group_by(tf, confidence) %>%
  add_count() %>%
  filter(n >= 10) %>%
  select(-n)

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
saveRDS(meta_expr, here("inst/benchdata/inputs/", "KnockTF_meta_data.rds"))

# clear env.
rm(list = ls())

# 3. Create Design JSONs ----




# 3.1. Dorothea + Dorothea Benchmark Data (DBD) ====
# 3.2. Dorothea + Knock_TF Data (KTF) ====

# 3.3. ChEA3 + DBD ====
# 3.4. ChEA3 + KTF ====

# 3.5. RegNetwork + DBD ====
# 3.6. RegNetwork + KTF ====

# 4. Benchmark ----

# 5. Summarize and Visualize output ----
