###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: compute genetic diversity
#DESCRIPTION: number of segregating sites, number of haplotypes, haplotype diversity,
# nucleotide diversity.
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(tidyverse)
library(pegas)
library(Biostrings)
library(DECIPHER)
library(xlsx)

# function to remove indels from msa
remove_indels_msa <- function(msa) {
  msa_m <- as.matrix(msa)
  tf_v <- apply(msa_m, 2, function(x) !any(grepl("-", x)))
  msa_mfilt <- msa_m[, tf_v]
  DNAStringSet(apply(msa_mfilt, 1, paste0, collapse = ""))
}

# genotypes 
gen <- readRDS("output/genotypes_rbal_poly.rds")
loci <- unique(gen$locus)
pop <- unique(gen$population)
msa <- readRDS("output/msa.rds")

# Per POPULATION ####
# compute hap nuc diversity per population
per_pop <-
  lapply(pop, function(y) {
    gen_pop <- filter(gen, population == y)
    lapply(loci, function(x) {
      j <-
        filter(gen_pop, locus == x) %>%
        pull(sequence) %>%
        Biostrings::DNAStringSet() %>%
        DECIPHER::AlignSeqs(gapOpening = c(-3, -1))
      z <- ape::as.DNAbin(j)
      w <- pegas::hap.div(z)
      y <- pegas::nuc.div(z)
      S <- length(ape::seg.sites(z, trailingGapsAsN = F, strict = T))
      k <-
        remove_indels_msa(j) %>%
        ape::as.DNAbin()
      S2 <-
        length(ape::seg.sites(k, trailingGapsAsN = F, strict = T))
      return(list(S = S,
                  S2 = S2,
                  hap_div = w,
                  nuc_div = y))
    }) %>%
      setNames(loci)
    }) %>%
  setNames(pop)
# number of haplotypes per population
nhap_pop <-
  gen |>
  group_by(population, locus) %>%
  summarize(nhap = length(unique(allele)))

# reformat hap, nuc div
hapnuc_div <-
  data.frame(value = unlist(per_pop)) %>%
  tibble::rownames_to_column("temp") %>%
  tidyr::separate_wider_delim(cols = temp, delim = ".",
                              names = c("population", "locus", "stat")) %>%
  tidyr::pivot_wider(id_cols = c(population,locus),
                       names_from = stat,
                       values_from = value) %>%
      as_tibble()
# join and reformat number of haplotypes, hap. div and nuc. div.
div_pop <-
  dplyr::left_join(hapnuc_div, nhap_pop) %>%
  plyr::dlply(~population, function(x) {
    x %>%
      dplyr::mutate(hap_div = round(hap_div, 2),
                    nuc_div1000 = round(nuc_div * 1000, 2)) %>%
      dplyr::select(locus, S, S2, nhap, hap_div, nuc_div1000)
  })

# GLOBAL ####
# compute hap nuc diversity
global_div <-
    lapply(loci, function(x) {
      j <-
        filter(gen, locus == x) %>%
        pull(sequence) %>%
        Biostrings::DNAStringSet() %>%
        DECIPHER::AlignSeqs(gapOpening = c(-3, -1))
      z <- ape::as.DNAbin(j)
      w <- pegas::hap.div(z)
      y <- pegas::nuc.div(z)
      S <- length(ape::seg.sites(z, trailingGapsAsN = F, strict = T))
      k <-
        remove_indels_msa(j) %>%
        ape::as.DNAbin()
      S2 <-
        length(ape::seg.sites(k, trailingGapsAsN = F, strict = T))
      return(list(S = S,
                  S2 = S2,
                  hap_div = w,
                  nuc_div = y))
      
    }) %>%
      setNames(loci)
# number of haplotypes
nhap_glob <-
  gen |>
  group_by(locus) |>
  summarize(nhap = length(unique(allele)))

# reformat hap, nuc div
hapnuc_div_global <-
  data.frame(value = unlist(global_div)) |>
  tibble::rownames_to_column("temp") |>
  tidyr::separate_wider_delim(cols = temp,
                              delim = ".",
                              names = c("locus", "stat")) |>
  tidyr::pivot_wider(id_cols = c(locus),
                     names_from = stat,
                     values_from = value) |>
  as_tibble()
# join and reformat number of haplotypes, hap. div and nuc. div.
div_pop_global <-
  dplyr::left_join(hapnuc_div_global, nhap_glob) |>
  dplyr::mutate(hap_div = round(hap_div, 2),
                    nuc_div1000 = round(nuc_div * 1000, 2)) |>
  dplyr::select(locus, S, S2, nhap, hap_div, nuc_div1000)

# write results to excel ####
excel_div <- "output/pegas_div.xlsx"
file.remove(excel_div)

# per pop results
names(div_pop) %>%
  lapply(function(x) {
    write.xlsx(div_pop[[x]],
               file = excel_div,
               sheetName = x,
               row.names = FALSE,
               append = T)
  })

# global results
write.xlsx(x = as.data.frame(div_pop_global),
               file = excel_div,
               sheetName = "all_pops",
               row.names = FALSE,
               append = T)
