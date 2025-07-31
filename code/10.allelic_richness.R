###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: allelic richness
#DESCRIPTION: allelic richness in R
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(pegas)
library(dplyr)

# genotypes 
gen <- readRDS("output/genotypes_rbal_poly.rds")

fx1 <- function(x) {
  y <- paste0(x, collapse = "/")
  return(y)
}
# transform tidy to hierstat
gen_hierf <-
  gen |>
  tidyGenR::gen_tidy2integers() |>
  select(population, sample, locus, allele) |>
  tidyr::pivot_wider(id_cols = c("population", "sample"),
                     names_from = "locus",
                     values_from = "allele",
                     values_fn = fx1) |>
  select(-sample) |>
  mutate_all(as.factor) |>
  as.loci()

# rarefy
rar_loc <-
  pegas::allelicrichness(gen_loci, method = "rarefaction") |>
  apply(1, function(x) round(sum(x), 1))
# raw allelic richness
raw_loc <-
  pegas::allelicrichness(gen_loci, method = "raw") |>
  apply(1, sum)

# save results
sink("output/allelic-richness.txt")
cat("Raw allelic richness:\n")
raw_loc

cat("\nRarefied allelic richness: \n")
rar_loc
sink()


