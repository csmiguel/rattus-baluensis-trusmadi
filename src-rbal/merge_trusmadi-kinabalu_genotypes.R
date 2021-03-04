###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: merge genotypes from Camacho-Sanchez et al. 2018 and trusmadi.
#PROJECT: rattus-baluensis-trusmadi
###.............................................................................
library(dplyr)

#read and format genotypes from trusmadi
trusmadi <-
  readRDS("data/intermediate/genotypes_notnull.rds")
trusmadi$asvs <-
  trusmadi$asvs %>%
    dplyr::mutate(population = "trusmadi",
                  allele = as.character(allele)) %>%
    dplyr::select(population, sample, locus, allele)


#read and format data from Kinabalu
path_alleles <- "etc/camacho-sanchez-2018/Sequences_introns_rbal.csv"
path_genotypes <- "etc/camacho-sanchez-2018/genotypes.csv"
# read allele sequences
kinabalu_al_seq <-
  read.csv(path_alleles, sep = ";") %>%
  dplyr::mutate(locus =
    stringr::str_remove(string = allele, pattern = ".$") %>% tolower()) %>%
  tibble::as_tibble() %>%
  setNames(names(trusmadi$al_seq))

# read genotypes
kinabalu_genotypes <-
  read.csv(path_genotypes, sep = ";") %>%
  dplyr::rename(sample = ID, population = Location) %>%
  dplyr::mutate(population = plyr::mapvalues(
    x = population,
    from = c("TAMBOYUCON", "KINABALU"),
    to = c("tambuyukon", "kinabalu"))) %>%
  dplyr::select(-Extraction) %>%
  reshape2::melt(id = c("population", "sample"),
                 value.name = "allele",
                 variable.name = "locus") %>%
  dplyr::mutate(locus = stringr::str_replace(locus, "_.$", "") %>%
                  tolower()) %>%
  dplyr::filter(!is.na(allele)) %>%
  tibble::as_tibble()
#join data
kinabalu <- dplyr::left_join(kinabalu_genotypes, kinabalu_al_seq,
                             by = c("locus", "allele"))
#set the equivalence between alleles from Camacho-Sanchez et al. 2018
#and alleles from trusmadi
allele_equivalence <-
  dplyr::full_join(kinabalu_al_seq, trusmadi$al_seq,
                 by = c("seq", "locus"),
                 suffix = c("_kinabalu", "_trusmadi")) %>%
  dplyr::arrange(locus) %>%
  dplyr::filter(!is.na(allele_trusmadi) & !is.na(allele_kinabalu))

# replace alleles from trusmadi that already exist
# in Kinabalu with Kinabalu's alleles' names.
trusmadi_al_seqs <-
  sort(unique(trusmadi$al_seq$locus)) %>%
  lapply(function(locusi) {
    x <- dplyr::filter(trusmadi$al_seq, locus == locusi)
    y <- dplyr::filter(allele_equivalence, locus == locusi)
  x$allele_renamed <-
    plyr::mapvalues(
      x = x$allele,
      from = y %>% dplyr::pull(allele_trusmadi),
      to = y %>% dplyr::pull(allele_kinabalu)
      )
  x
  }) %>% do.call(what = rbind)


trusmadi_renamed <-
  dplyr::left_join(trusmadi$asvs, trusmadi_al_seqs,
                 by = c("locus", "allele")) %>%
  dplyr::select(-allele) %>%
  dplyr::rename(allele = allele_renamed)

#bind kinabalu and trusmandi
genotypes <- rbind(trusmadi_renamed, kinabalu)

#remove loci which were not concatenated in trusmadi
loc_2_remove <-
  genotypes %>%
  dplyr::filter(grepl(pattern = "NNN", x = seq)) %>%
  pull(locus) %>% unique()

#filter locus to remove
genotypes <-
  genotypes %>%
  dplyr::filter(!(locus %in% loc_2_remove))
#save genotypes
saveRDS(genotypes, "output/trusmadikinabalu_genotypes.rds")
