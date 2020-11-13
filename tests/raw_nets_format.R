# RegNetwork ----
high_confidence = read.csv(here("inst/testdata/inputs/regnetwork",
                                "high_confidence.csv"), header = T)
medium_confidence = read.csv(here("inst/testdata/inputs/regnetwork",
                                  "medium_confidence.csv"), header = T)
low_confidence = read.csv(here("inst/testdata/inputs/regnetwork",
                               "low_confidence.csv"), header = T)

# bind and format
regnetwork = rbind(high_confidence, medium_confidence, low_confidence) %>%
  filter(!str_detect(regulator_symbol, "miR|mir|hsa")) %>%
  filter(!str_detect(target_symbol, "miR|mir|hsa")) %>%
  rename(tf = regulator_symbol, target = target_symbol) %>%
  mutate(likelihood = 1, mor = 1) %>%
  select(tf, target, confidence,  likelihood, mor)
head(regnetwork)

data_loc <- "~/PhD/data/decoupleR_bench_data/inputs/networks/"
saveRDS(regnetwork, paste(data_loc, "regnetwork_raw.rds", sep=""))

# chEA3
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

# get chEA3 network
chea3_network <- load_chea3("all")  %>%
  select(-experiment)
saveRDS(chea3_network, paste(data_loc, "chea3_raw.rds", sep=""))
head(chea3_network)
head(regnetwork)
