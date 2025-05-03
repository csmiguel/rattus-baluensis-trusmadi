###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com 
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# May 2024
###.............................................................................
#GOAL: rename alleles from Trusmadi based on alleles from Camacho-Sanchez et al. 2018.
# plus merge tidy results from Trusmadi and Kinabalu.
#PROJECT: rbaluensis_trusmadi
###.............................................................................
library(dplyr)
library(tidyr)
library(stringr)
library(tidyGenR)

# load tidy results from trusmadi
gen <- readRDS("data/raw/genotypes.rds")

# tidy data from Trusmadi
tidy_trusmadi <-
  gen |>
  dplyr::mutate(population = "trusmadi") |>
  dplyr::select(population, sample, locus, allele, allele_no, sequence)

# read and format allele and genotypes data from Camacho-Sanchez et al. 2018
# read allele sequences
alleles2018 <- "data/raw/camacho-sanchez-2018/Sequences_introns_rbal.csv"
kinabalu_al <-
  read.csv(alleles2018, sep = ";") %>%
  dplyr::mutate(locus =
                  stringr::str_remove(string = allele, pattern = ".$") %>%
                  tolower()) %>%
  tibble::as_tibble()

# read genotypes
genotypes2018 <- "data/raw/camacho-sanchez-2018/genotypes.csv"
kinabalu_geno <-
  read.csv(genotypes2018, sep = ";") |>
  dplyr::rename(sample = ID, population = Location) |>
  dplyr::mutate(population = plyr::mapvalues(
    x = population,
    from = c("TAMBOYUCON", "KINABALU"),
    to = c("tambuyukon", "kinabalu"))) |>
  dplyr::select(-Extraction) |>
  pivot_longer(cols = -c(population, sample),
               names_to = "locus",
               values_to = "allele",
               values_drop_na = T) |>
  dplyr::mutate(locus = str_replace(locus, "_.$", "") |>
                  tolower()) 
# tidy kinabalu data
tidy_kinabalu <-
  dplyr::left_join(kinabalu_geno, kinabalu_al,
                  by = c("locus", "allele")) |>
  plyr::ddply(~sample + locus, function(x) mutate(x, allele_no = 1:nrow(x))) |>
  as_tibble

# rename alleles from trusmadi based on kinabalu
tidy_trusmadi_renamed <-
  tidyGenR::rename_alleles(df_alleles = tidy_trusmadi, 
                           ref_alleles = tidy_kinabalu)

# bind kinabalu and trusmandi
tidy_rbal <-
  rbind(tidy_trusmadi_renamed, tidy_kinabalu)
saveRDS(tidy_rbal, "output/genotypes_rbal.rds")

# remove monomorphic and fetub
tidy_rbal_poly <-
  tidyGenR::remove_monomorphic(tidy_rbal) |>
  dplyr::filter(locus != "fetub")
# save polimorphic loci without fetub.
saveRDS(tidy_rbal_poly, "output/genotypes_rbal_poly.rds")
