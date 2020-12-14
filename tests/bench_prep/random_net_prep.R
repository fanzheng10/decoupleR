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
