# 0. Load prerequisites ----
# Directories
library(decoupleR)
library(ggplot2)
library(pheatmap)

bench_input <- system.file("benchdata", "inputs", package = "decoupleR")
bench_output <- system.file("benchdata", "outputs", package = "decoupleR")

# Load RDS with all designs
design_all <- readRDS(file.path(bench_input, "design_all.rds"))


# 1. Run Dorothea DBD with all stats (- fgsea) and network combinations ----
# Filter runs only with dorothea benchmark
design_dbd <- design_all %>%
  filter(endsWith(row_name,"dbd"))

# dbd_res <- run_benchmark(design_dbd)
# saveRDS(dbd_res, file.path(bench_output, "dbd_res.rds"))

# read Res
dbd_res <- readRDS(file.path(bench_output, "dbd_res.rds"))
print(dbd_res@bench_res, n=100)
print(dbd_res@summary$auroc_summary, n=100)
dbd_res@summary$auroc_heat



# 2. Run Knock TF benchmark with all combs ----
design_ktf <- design_all %>%
  filter(endsWith(row_name,"ktf"))

# check design
design_ktf

# Run benchmark
# ktf_res <- run_benchmark(design_ktf)
# saveRDS(ktf_res, file.path(bench_output, "ktf_res.rds"))
ktf_res <- readRDS(file.path(bench_output, "ktf_res.rds"))
print(ktf_res@bench_res, n=100)
print(ktf_res@summary$auroc_summary, n=100)
ktf_res@summary$auroc_heat
# Note that infinite filtering was introduced because  mean, normalized mean
# and scira return infinite values with KTF that need to be filtered

## 2.1. Example where I run it without excluding infinities as part of the
# pipeline, but rather handle them ater the run
# 2.* Handle issues post-run example ====
# load results
ktf_res_inf <- readRDS(file.path(bench_output, "ktf_res_inf.rds"))

# Check infinite values
row_infinites <- ktf_res_inf@bench_res %>%
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
ktf_res@bench_res <- ktf_res_inf@bench_res %>%
  bench_format() %>%
  mutate(activity = activity %>%
           map(function(tib) tib %>%
                 mutate_at(vars(score), ~replace(., is.infinite(.), 0))
           )) %>%
  mutate(roc = activity %>% map(calc_roc_curve))

# assign summary
ktf_res@summary <- ktf_res@bench_res %>%
  bench_sumplot()
# Save Object
# saveRDS(ktf_res, file.path(bench_output, "ktf_res.rds"))

# read results
ktf_res <- readRDS(file.path(bench_output, "ktf_res.rds"))
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
# ktf_rank_res <- run_benchmark(design_ktf_rank)
# saveRDS(ktf_rank_res, file.path(bench_output, "ktf_rank_res.rds"))

# Read results
ktf_rank_res <- readRDS(file.path(bench_output, "ktf_rank_res.rds"))
print(ktf_rank_res@bench_res, n=100)
print(ktf_rank_res@summary$auroc_summary, n=100)
ktf_rank_res@summary$auroc_heat


# 4. Run with Random regulons for Dorothea regulons ----
# Random regulons should have the same gene coverage and exclude the TP genes

# 4.1. Dorothea vs Random Dorothea DBD RUN ----
# Filter run dorothea only with both benchmarks
design_doro <- design_all %>%
  filter(str_detect(row_name,"dorothea"))

# check if correct
design_doro$row_name

# replace network
design_rand <- design_doro %>%
  mutate(net_loc = file.path(bench_input, "networks",
                             "dorothea_random.rds"),
         row_name = list("dor_random_dbd", "dor_random_dbd", "dor_random_dbd",
                         "dor_random_ktf", "dor_random_ktf", "dor_random_ktf")) %>%
  unnest(row_name)

# check if OK
View(design_rand)

# # modify design_dbd
# dor_rand_res <- run_benchmark(design_rand)
# saveRDS(dor_rand_res, file.path(bench_output, "dor_rand_res.rds"))

# Read Results
dor_rand_res <- readRDS(file.path(bench_output, "dor_rand_res.rds"))
print(dor_rand_res@bench_res, n=100)
print(dor_rand_res@summary$auroc_summary, n=100)
dor_rand_res@summary$auroc_heat



# 4.2. Run with BiRewire Net ----
# Filter run dorothea only with both benchmarks
design_doro <- design_all %>%
  filter(str_detect(row_name,"dorothea"))
 run_benchmark(design_birewire)
# check if correct

# replace network and names
design_birewire <- design_doro %>%
  mutate(net_loc = file.path(bench_input,
                             "networks", "dorothea_birewire.rds"),
         row_name = list("dor_birew_dbd", "dor_birew_dbd", "dor_birew_dbd",
                         "dor_birew_ktf", "dor_birew_ktf", "dor_birew_ktf")) %>%
  unnest(row_name)

# dor_birewire_res <- run_benchmark(design_birewire)
# saveRDS(dor_birewire_res, file.path(bench_output, "dor_birewire_res.rds"))
dor_birewire_res <- readRDS(file.path(bench_output, "dor_birewire_res.rds"))
print(dor_birewire_res@bench_res, n=100)
dor_birewire_res@summary$auroc_heat
