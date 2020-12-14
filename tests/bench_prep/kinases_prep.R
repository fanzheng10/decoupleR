# 5. Kinases Prep ----
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
