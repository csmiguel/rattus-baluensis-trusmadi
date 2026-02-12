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
library(dplyr)

# 1. write genotypes to str file
meta <- read.csv("data/raw/metadata.csv")
gen <- readRDS("output/genotypes_rbal_poly.rds")

# subset samples into subsets 1 and 2
samples_trusmadi <-
  filter(gen, population == "trusmadi") |> pull(sample) |> unique()

subset1 <-
  sample(samples_trusmadi, length(samples_trusmadi) / 2, replace = F)
subset2 <-
  samples_trusmadi[!samples_trusmadi %in% subset1]
assertthat::assert_that(length(intersect(subset1, subset2)) == 0)

subsets <-
  list(subset1 = subset1,
     subset2 = subset2)

# create structure files for each subset
seq_along(subsets) |>
  lapply(function(i) {
    # name vars
    names_i <- names(subsets)[i]
    subset_i <- subsets[[names_i]]
    # create subset dir
    str_outdir <- file.path("data/intermediate/structure-subsets", names_i)
    dir.create(str_outdir)
    # format str file
    gen_i<-
      gen |>
      filter((population != "trusmadi") | (sample %in% subset_i))
    gen_tidy2wide(gen_i) |>
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
    path_to_structure = "/home/miguel/miniconda3/envs/structure_env/bin/structure"
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
                                  length(unique(gen_i$sample)))) %>%
      stringr::str_replace("#define NUMLOCI.*CHANGEME",
                           paste0("#define NUMLOCI", "\t",
                                  length(unique(gen_i$locus)))) %>%
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
    runstr <- glue("structure_threader run -K {K} -R {runsPerK} -i $PWD/{infile} -o $PWD/{outpath} -t {num_of_threads} -st {path_to_structure} &")
    writeLines(runstr, con = file.path(outpath, "run-exec"))
    
  })



