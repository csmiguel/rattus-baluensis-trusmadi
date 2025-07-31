###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# June 2025
###.............................................................................
#GOAL: Venn diagram with shared / unique alleles per population
#PROJECT: rbaluensis-trusmadi
###.............................................................................

library(ggVennDiagram)
library(ggplot2)
library(tibble)
library(dplyr)

# path to clean polymorphic genotypes

p_gen <- "output/genotypes_rbal_poly.rds"

# wrangle genotypes
gen <-
  readRDS(p_gen) |>
  select(population, locus, allele) |>
  distinct() |>
  rowwise() |>
  mutate(al_named = paste(locus, allele, collapse = "_")) |>
  select(population, al_named)

# vector with alleles for each population
tru <-
  filter(gen, population == "trusmadi") |>
  pull(al_named)
kin <-
  filter(gen, population == "kinabalu") |>
  pull(al_named)
tam <-
  filter(gen, population == "tambuyukon") |>
  pull(al_named)

list_al <-
  list(Trusmadi = tru,
       Kinabalu = kin,
       Tambuyukon = tam)
# plot
p1 <-
  ggVennDiagram::ggVennDiagram(list_al, label = "count") +
  scale_fill_gradient(low = "#F4FAFE", high = "grey20")
ggsave(filename = "output/venn_alleles.pdf",
       p1,
       width = 5,
       height = 5)
