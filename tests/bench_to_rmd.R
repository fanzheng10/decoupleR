# 0. Load prerequisites ----
# Directories
input_dir <- system.file("extdata", package = "decoupleR")
bench_input <- system.file("benchdata", "inputs", package = "decoupleR")
bench_output <- system.file("benchdata", "outputs", package = "decoupleR")

# Load RDS with all designs
design_all <- readRDS(file.path(input_dir, "design_all.rds"))


# 1. Run Dorothea DBD with all stats (- fgsea) and network combinations ----
# Filter runs only with dorothea benchmark
design_dbd <- design_all %>%
  filter(endsWith(row_name,"dbd"))

# design_dbd_res <- run_benchmark(design_dbd)
# saveRDS(design_dbd_res, file.path(bench_output, "design_dbd_res.rds"))

# exchange for example data
# design_d$bnch_expr <- file.path(bench_input, "input-dorothea_bench_example.rds")
# xd <- run_benchmark(design_d[-3,])

design_dbd_res <- readRDS(file.path(bench_output, "design_dbd_res.rds"))

print(design_dbd_res@bench_res, n=100)
print(design_dbd_res@summary$auroc_summary, n=100)
design_dbd_res@summary$auroc_heat



# 2. Run Knock TF benchmark with all combs ----
design_ktf <- design_all %>%
  filter(endsWith(row_name,"ktf"))

# check design
design_ktf

# Run benchmark
# In this case we run it without formatting, as mean, normalized mean, and
# scira contain infinite values that need to be filtered
# design_ktf_res <- run_benchmark(design_ktf, .form = F, .perform = F)
# saveRDS(design_ktf_res, file.path(bench_output, "design_ktf_res.rds"))

# load results
design_ktf_res <- readRDS(file.path(bench_output, "design_ktf_res_inf.rds"))

# Check infinite values
row_infinites <- design_ktf_res@bench_res %>%
  bench_format() %>%
  mutate(ininite_no = activity %>%
           map(function(tib) sum(is.infinite(tib$score)))
  ) %>%
  select(row_name, lvls, statistic, ininite_no) %>%
  unnest(ininite_no)

print(row_infinites, n=100)
# Looks like the numbers are quite low, given that there are between
# 60,000 and 950,000 values in the different rows/runs
# As such, these infinites filter will be introduced in the benchmark
# Maybe it's worth thinking what causes them

# format and replace
design_ktf_res@bench_res <- design_ktf_res@bench_res %>%
  bench_format() %>%
  mutate(activity = activity %>%
           map(function(tib) tib %>%
                 mutate_at(vars(score), ~replace(., is.infinite(.), 0))
           )) %>%
  mutate(roc = activity %>% map(calc_roc_curve))



# assign summary
design_ktf_res@summary <- design_ktf_res@bench_res %>%
  bench_sumplot()
# Save Object
# saveRDS(design_ktf_res, file.path(bench_output, "design_ktf_res.rds"))

# read results
design_ktf_res <- readRDS(file.path(bench_output, "design_ktf_res.rds"))

# Show Result summary
print(design_ktf_res@bench_res, n=100)
print(design_ktf_res@summary$auroc_summary, n=100)
design_ktf_res@summary$auroc_heat




# 3. Run Knock TF (Rank-based Activity) benchmark with all combs ----
design_ktf_rank <- design_all %>%
  filter(endsWith(row_name,"ktf"))

# Check design
design_ktf_rank

# re-assign to p-value-rank-based TF activity
design_ktf_rank$bnch_expr <- file.path(bench_input, "KnockTF_rank_expr.rds")

# run
# design_ktf_rank_res <- run_benchmark(design_ktf_rank)
# saveRDS(design_ktf_rank_res, file.path(bench_output, "design_ktf_rank_res.rds"))

# Read results
design_ktf_rank_res <- readRDS(file.path(bench_output, "design_ktf_rank_res.rds"))
print(design_ktf_rank_res@bench_res, n=100)
print(design_ktf_rank_res@summary$auroc_summary, n=100)
design_ktf_rank_res@summary$auroc_heat


# 4. Create Random regulons for each TF ----
# Random regulons should have the same gene coverage and exclude the TP genes

# 4.1. Random regulons for Dorothea (My approach) ====
dorothea_filtered <- readRDS(file.path(bench_input,
                                       "networks/dorothea_filtered.rds"))


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

# check if appoints regulon with non-intersect gene of the same size
intersect(ADNP_regulon$target, random_ADNP_regulon$target)
length(ADNP_regulon$target)
length(random_ADNP_regulon$target)

random_network <- readRDS(file.path(bench_input, "networks/random_dorothea.rds"))



random_network <- arrange(random_network)
random_network$confidence <- dorothea_net$confidence
random_network$likelihood <- dorothea_net$mor
random_network$mor <- dorothea_net$mor


saveRDS(random_network, file.path(bench_input, "networks/random_dorothea.rds"))


# 4.2.1. Dorothea vs Random Dorothea DBD RUN ----
# Filter run dorothea only with both benchmarks
design_doro <- design_all %>%
  filter(str_detect(row_name,"dorothea"))

# check if correct
design_doro$row_name

# replace network
design_rand <- design_doro %>%
  mutate(net_loc = file.path(bench_input, "networks/random_dorothea.rds"),
         row_name = list("dor_random_dbd", "dor_random_dbd", "dor_random_dbd",
                         "dor_random_ktf", "dor_random_ktf", "dor_random_ktf")) %>%
  unnest(row_name)


# modify design_dbd
design_rand_res <- run_benchmark(design_rand)
saveRDS(design_rand_res, file.path(bench_output, "design_ktf_res.rds"))


design_rand@bench_res
design_rand@summary$auroc_summary


# 4.3. Random regulons for Dorothea with BiRewire ----
library(BiRewire)

dorothea_net <- readRDS("/home/dbdimitrov/Repos/decoupleR/inst/benchdata/inputs/networks/dorothea_filtered.rds")

xd <- dorothea_net %>%
  select(tf, target, mor) %>%
  pivot_wider(names_from = target, values_from = mor)

xd[is.na(xd)] <- 0

dorothea_bi <- xd %>% as.data.frame()

# -> replace BRCA

###################################################
### code chunk number 2: GetABipartiteGraph
###################################################
data(BRCA_binary_matrix)##loads an binary genomic event matrix for the
##breast cancer dataset
g=birewire.bipartite.from.incidence(BRCA_binary_matrix)##models the dataset
## as igraph bipartite graph


###################################################
### code chunk number 3: PerformAnalisys
###################################################
step=5000
max=100*sum(BRCA_binary_matrix)
scores<-birewire.analysis.bipartite(BRCA_binary_matrix,step,
                                    verbose=FALSE,max.iter=max,n.networks=5,display=F)





###################################################
### code chunk number 4: PerformAnalisysUndirected
###################################################
g.und<-erdos.renyi.game(directed=F,loops=F,n=1000,p.or.m=0.01)
m.und<-get.adjacency(g.und,sparse=FALSE)
step=100
max=100*length(E(g.und))
scores.und<-birewire.analysis.undirected(m.und,step=step,
                                         verbose=FALSE,max.iter=max,n.networks=5)



###################################################
### code chunk number 5: Rewire
###################################################
m2<-birewire.rewire.bipartite(BRCA_binary_matrix,verbose=FALSE)
g2<-birewire.rewire.bipartite(g,verbose=FALSE)

g2

###################################################
### code chunk number 6: RewireUndirected
###################################################
m2.und<-birewire.rewire.undirected(m.und,verbose=FALSE)
g2.und<-birewire.rewire.undirected(g.und,verbose=FALSE)

###################################################
### code chunk number 7: Similarity
###################################################
sc=birewire.similarity(BRCA_binary_matrix,m2)
sc=birewire.similarity(BRCA_binary_matrix,t(m2))#also works


