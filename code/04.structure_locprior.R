###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: write structure input files and prepare parameters files for 
# structure_threader
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(tidyGenR)
library(glue)
library(stringr)
library(adegenet)

# 1. write genotypes to str file
meta <- read.csv("data/raw/metadata.csv")
gen <- readRDS("output/genotypes_rbal_poly.rds")
str_outdir <- "data/intermediate/structure-locprior"
dir.create(str_outdir)
gen_tidy2wide(gen) |>
  tidyGenR::gen_wide2structure(write_out = file.path(str_outdir, "structure.str"),
                             popdata = meta)

# RUN params for structure_threader
K = 7 # from 1 to k
infile = file.path(getwd(), str_outdir, "structure.str")
outpath = file.path(getwd(), str_outdir)
runsPerK = 4
MCMC = 200000
burnin = 50000
num_of_threads = 5
path_to_structure = "/usr/local/bin/structure"
# model
admixture = 1

# generate empty params files
system(glue("structure_threader params -o {outpath}"))

# edit params
# main params
mainparams <- file.path(outpath, "mainparams")
mp <- readLines(mainparams)
mp_edited <-
  mp |>
  stringr::str_replace("#define NUMINDS.*CHANGEME",
                       paste0("#define NUMINDS", "\t",
                              length(unique(gen$sample)))) %>%
  stringr::str_replace("#define NUMLOCI.*CHANGEME",
                       paste0("#define NUMLOCI", "\t",
                              length(unique(gen$locus)))) %>%
  stringr::str_replace("#define MAXPOPS.*CHANGEME",
                       paste0("#define MAXPOPS", "\t", K)) |>
  stringr::str_replace("#define BURNIN\\s*[0-9]*",
                       paste0("#define BURNIN", "\t", format(burnin, scientific = FALSE))) |>
  stringr::str_replace("#define NUMREPS\\s*[0-9]*",
                       paste0("#define NUMREPS", "\t", format(MCMC, scientific = FALSE)))

writeLines(mp_edited, con = mainparams)

# extraparams
extraparams <- file.path(outpath, "extraparams")
ep <- readLines(extraparams)
ep_edited <-
  ep |>
  stringr::str_replace("#define LOCPRIORINIT\\s*0",
                       paste0("#define LOCPRIORINIT", "\t", "1.0")) |>
  stringr::str_replace("#define LOCPRIOR\\s*0",
                       paste0("#define LOCPRIOR", "\t", "1")) |>
  stringr::str_replace("#define FREQSCORR\\s*0",
                       paste0("#define FREQSCORR", "\t", "1")) |>
  stringr::str_replace("#define INFERALPHA\\s*0",
                       paste0("#define INFERALPHA.", "\t", "1")) |>
  stringr::str_replace("#define RANDOMIZE\\s*0",
                       paste0("#define RANDOMIZE", "\t", "1")) |>
  stringr::str_replace("#define LOCISPOP\\s*0",
                       paste0("#define LOCISPOP", "\t", "1"))

writeLines(ep_edited, con = extraparams)

# run
runstr <- glue("structure_threader run -K {K} -R {runsPerK} -i {infile} -o {outpath} -t {num_of_threads} -st {path_to_structure}")
writeLines(runstr, con = file.path(outpath, "run-exec"))

