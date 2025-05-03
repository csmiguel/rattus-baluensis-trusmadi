###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: add genotype data to sf1.xls
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(tidyGenR)
library(tidyverse)
library(xlsx)

data("primers")
pr <-
  primers %>%
  arrange(locus) %>%
  {.["locus"]}

# number of raw called alleles  per locus
called_alleles_rawgen <-
  readRDS("output/genotypes_rbal.rds") %>%
  filter(population == "trusmadi") %>%
  group_by(locus) %>%
  dplyr::summarize(raw_alleles = n())

# number of called alleles per locus used for popgen
# ie: after removing monomorphic and fetub.
called_alleles_20filt <-
  readRDS("output/genotypes_rbal_poly.rds") %>%
  filter(population == "trusmadi") %>%
  group_by(locus) %>%
 dplyr::summarize(filt_alleles = n())


# join data
df <- 
  dplyr::left_join(pr, called_alleles_rawgen, by = "locus") %>%
  dplyr::left_join(called_alleles_20filt, by = "locus") 
  
# write to excel
excel_path <- "output/suppltables.xlsx"
write.xlsx(x = df,
           file =  excel_path,
           sheetName = "ST5.tracked-genotypes",
           append = T)

# append genotypes
gen <-
  readRDS("output/genotypes_rbal_poly.rds") %>%
  tidyGenR::gen_tidy2wide()

write.xlsx(x = gen,
           file =  excel_path,
           sheetName = "ST6.genotypes",
           append = T)
# append alleles
vdigest <- Vectorize(digest::digest)
al <-
  readRDS("output/genotypes_rbal_poly.rds") |>
  dplyr::mutate(md5 = vdigest(sequence, "md5"),
                nt = nchar(sequence)) |>
  dplyr::select(locus, allele, nt, md5, sequence) |>
  dplyr::distinct()
write.xlsx(x = al,
           file =  excel_path,
           sheetName = "ST7.alleles",
           append = T)
