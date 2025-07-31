##

library(tidyGenR)
library(Biostrings)
library(dplyr)
# read genotypes from polymorphic locy
gen <-
  readRDS("output/genotypes_rbal_poly.rds")

# list of loci
loci <- unique(gen$locus)
msa <- readRDS("output/msa.rds")

# start writing text file
p_dir_popabc <- "data/intermediate/popabc"
if(!dir.exists(p_dir_popabc))
  dir.create(p_dir_popabc)
# read msa
msa <- readRDS("output/msa.rds")
# list of loci
loci <- unique(gen$locus)

lapply(loci, function(x) {
  msa_i <- msa[[x]]
  # remove gaps
  vector_selection <- !apply(as.matrix(msa_i) == "-", 2, any)
  seqs <- apply(as.matrix(msa_i)[, vector_selection], 1, paste0, collapse = "")
  msa2 <- DNAStringSet(seqs)
  # get populations
  pops <- stringr::str_split_i(names(msa2), "_", 2)
  # create output file
  p_nex_i <- file.path(p_dir_popabc, paste0(x, ".nex"))
  if(file.exists(p_nex_i))
    file.remove(p_nex_i)
  t1 <- tempfile()
  out_popart(msa2, blocks = "DATA", outnex = t1)
  tx <- readLines(t1)
  tx <- stringr::str_remove(tx, stringr::fixed(" MISSING=? GAP=- INTERLEAVE=NO"))
  writeLines(tx, p_nex_i)
  con <- file(p_nex_i, open = "a")  # 'a' = append mode
  # loop to write SETS block with pop info
  writeLines(
    text = c("BEGIN SETS;",
             paste0("  TAXSET ", "'", levels(as.factor(pops))[1], "' = ",
                    paste(which(pops == levels(as.factor(pops))[1]),
                          collapse = " "), ";"),
             paste0("  TAXSET ", "'", levels(as.factor(pops))[2], "' = ",
                    paste(which(pops == levels(as.factor(pops))[2]),
                          collapse = " "), ";"),
             paste0("  TAXSET ", "'", levels(as.factor(pops))[3], "' = ",
                    paste(which(pops == levels(as.factor(pops))[3]),
                          collapse = " "), ";"),
             "END;"),
    con = con)
  close(con)
  })
