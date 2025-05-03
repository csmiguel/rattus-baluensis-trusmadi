###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: format loci MSA for PopART
#DESCRIPTION: input are tidy data and MSA.
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(ape)
library(tidyGenR)
library(tidyverse)

# genotypes t+k+t
gen <-
  readRDS("output/genotypes_rbal_poly.rds")

# column names to use in fasta headers
seqname_vars <- c("sample", "locus", "allele", "allele_no")

# fasta header constructrion
fh <- paste0("{",
            paste(seqname_vars, collapse = "}_{"), 
            "}"
            )
# add col 2 df
df1 <-
  gen %>%
  mutate(seqname = glue::glue(fh) %>% as.character)

# biostring MSA plus correct names to match with traits
msa <-
  readRDS("output/msa_poli.rds") %>%
  lapply(function(x) {
    nx <-
      names(x) %>%
      stringr::str_remove("_tambuyukon|_kinabalu|_trusmadi")
  names(x) <- nx
  return(x)
})
msa <-
  readRDS("output/msa_poli.rds") %>%
  lapply(function(x) {
    nx <-
      names(x) %>%
      stringr::str_remove("_tambuyukon|_kinabalu|_trusmadi")
    names(x) <- nx
    return(x)
  })

# per locus
destdir <- "output/nexpopart"
dir.create(destdir)
names(msa) %>%
  lapply(function(x) {
    cat(x, "\n")
    out_popart(msa = msa[[x]],
               x = filter(df1, locus == x),
               sname = "seqname",
               xgroups = "population",
               outnex = file.path(destdir, paste0(x, ".nex"))
    )
    }
    )
