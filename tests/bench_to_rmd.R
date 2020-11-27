# Load prerequisites ----
input_dir <- system.file("extdata", package = "decoupleR")
bench_dir <- system.file("benchdata", "inputs", package = "decoupleR")
bench_output <- system.file("benchdata", "outputs", package = "decoupleR")

# Load RDS with all designs  ----
design_all <- readRDS(file.path(input_dir, "design_all.rds"))


# 1. Run Dorothea DBD with all stats (- fgsea) and network combinations ----
# Filter runs only with dorothea benchmark
design_dbd <- design_all %>%
  filter(endsWith(row_name,"dbd"))

# design_dbd_res <- run_benchmark(design_dbd)
# saveRDS(design_dbd_res, file.path(bench_output, "design_dbd_res.rds"))

# exchange for example data
# design_d$bnch_expr <- file.path(bench_dir, "input-dorothea_bench_example.rds")
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
design_ktf_rank$bnch_expr <- file.path(bench_dir, "KnockTF_rank_expr.rds")

# run
# design_ktf_rank_res <- run_benchmark(design_ktf_rank)
# saveRDS(design_ktf_rank_res, file.path(bench_output, "design_ktf_rank_res.rds"))

# Read results
design_ktf_rank_res <- readRDS(file.path(bench_output, "design_ktf_rank_res.rds"))
print(design_ktf_rank_res@bench_res, n=100)
print(design_ktf_rank_res@summary$auroc_summary, n=100)
design_ktf_rank_res@summary$auroc_heat


# 4. Investigate the effect of Regulon size and confidence levels
