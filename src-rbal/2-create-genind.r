###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: merge genotypes from Camacho-Sanchez et al. 2018 and trusmadi into
# a genind object
#PROJECT: rattus-baluensis-trusmadi
###.............................................................................
library(dplyr)
library(adegenet)

#read tidy gentoypes
gen <-
  readRDS("output/trusmadikinabalu_genotypes.rds") %>%
  dplyr::select(-seq)

ploidy <- 2
#reshape genotype data to:
h <- reshape2::acast(data = gen,
                     formula = population + sample ~ locus + allele,
                     fun.aggregate = length,
                     value.var = "allele")
colnames(h) <- stringr::str_replace(colnames(h), "_", ".")


#create population vector
# assign a ficticious unique population to the data if no pop data is provided
pops <-
    rownames(h) %>%
    stringr::str_remove("_.*$")
#store genotypes in genind object
hh <- new("genind",
          tab = h,
          type = "codom",
          ploidy = ploidy,
          pop = pops)

saveRDS(hh, "output/genind-trusmadikinabalu.rds")
