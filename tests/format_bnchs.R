
# 5.* Convert Dorothea and KnockTF Bnch ----
# ** Dorothea ====
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

