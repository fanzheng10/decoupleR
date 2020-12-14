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
