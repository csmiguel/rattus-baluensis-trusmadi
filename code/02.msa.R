###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: create multiple sequence alignments for all the alleles in each locus
#DESCRIPTION: input is tidy genotypes, output are a fasta msa for each locus written
# to 'output' and an a list of R DNAStringSet with alignments for each locus.
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(plyr)
library(dplyr)
library(DECIPHER)
library(Biostrings)

# tidy results for all individuals (kinabalu + trusmadi)
gen <-
  readRDS("output/genotypes_rbal_poly.rds")
# create dir to write msa
fp <- "data/intermediate/msa"
file.remove(list.files(fp, full.names = T))
dir.create(fp)

loci <- gen[["locus"]] %>% unique()

# save msa of all alleles for each locus
msa_al <-
  loci |>
  lapply(function(x) {
    cat(x, "\n")
    msal <-
      filter(gen, locus == x) %>%
      mutate(nameseq = paste(sample, population, locus, allele, allele_no, sep = "_")) %>%
      pull(sequence, name = nameseq) %>%
      Biostrings::DNAStringSet(use.names = T) %>%
      DECIPHER::AlignSeqs(gapOpening = c(-3, -1))
    Biostrings::writeXStringSet(msal,
                                filepath = file.path(fp, paste0(x, ".fasta")))
    return(msal)
  }) |>
  setNames(loci)
# save html of msa of unique allele per locus
loci |>
  lapply(function(x) {
    cat("Save hmtl per MSA, only unique seqs.", x, "\n")
      filter(gen, locus == x) |>
      dplyr::distinct(locus, allele, sequence) |>
      pull(sequence, name = allele) %>%
        {
          if (length(.) == 1) {
            NULL
            cat(x, "has one seq.")
            } else if (length(.) > 1) {
            cat(x, "has more than one seq.")
              Biostrings::DNAStringSet(., use.names = T) |>
              DECIPHER::AlignSeqs(gapOpening = c(-3, -1)) |>
              DECIPHER::BrowseSeqs(colWidth = 100,
                htmlFile = file.path(fp, paste0(x, ".html")),
                openURL = F,
                title = x)
          }
          }
      })

# save 
saveRDS(msa_al, "output/msa.rds")
