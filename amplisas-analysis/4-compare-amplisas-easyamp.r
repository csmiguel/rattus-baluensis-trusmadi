setwd("/projects/rattus-baluensis-trusmadi/")
#read ampliSAS genotypes and compare them to those from easyamp
library(dplyr)
library(magrittr)
#load functions
source("src/functions/genotyping.r")
source("src/parameters/parameters.r")
source("src/functions/synonymize_alleles.r")
source("src/functions/generate_output1.r")

# files with amplisas output
amplisas_output <- list.files("amplisas-analysis/results_amplisas/allseqs",
                              ".txt$",
                              full.names = T)
# names of loci
loci <-
  amplisas_output %>% 
  stringr::str_remove("amplisas-analysis/results_amplisas/allseqs/") %>%
  stringr::str_remove(".txt")
# names of alleles
al_names_all <- as.character(101:999)
  
#1. read per-marker genotypes and tidy data for amplisas
amplisas <-
  loci %>% #for each locus
    lapply(function(locus) {
      h <-
        read.table(paste0("amplisas-analysis/results_amplisas/allseqs/", locus, ".txt"),
                 sep = "\t",
                 header = T) %>% # read amplisas output
        #reformat amplisas output
        dplyr::select(-MD5, -LENGTH, -DEPTH, -SAMPLES) %>%
        dplyr::rename(alleles = NAME,
                      seq = SEQ) %>%
        dplyr::mutate(locus = stringr::str_split(alleles, "-")[[1]][1],
                      allele = stringr::str_split(alleles, "-")[[1]][2]) %>%
        dplyr::select(-alleles) %>%
        dplyr::as_tibble()
      h %>%
        dplyr::select(-seq, -locus, -allele) %>%
        t() %>%
        tidyr::replace_na(0) %>%
        {attr(., "dimnames")[[2]] <- h$seq; .} %>%
        # apply filtering as in easyamp
        #clean low freq variants
        clean_lowfrequencyvariants() %>%
        # rename alleles
        rename_alleles() %>%
        # remove paralogues
        remove_nontarget(locus_i = locus) %>%
        # tidy genotypes
        tidy_asvs() %>%
        # remove_extra_alleles
        remove_extra_alleles(locus_i = locus) %>%
        # call homozygous
        call_homozygous() %>%
        # add locus name to data
        add_locus_name(locus_i = locus)
    }) %>% setNames(loci)

# reformar genotypes amplisas: from lists of loci to dataframes of alleles and genotypes
gen <-
  amplisas %>%
  plyr::ldply(function(x) x$asvs, .id = NULL) %>% tibble::as_tibble()
#bind allele data
al <-
  amplisas %>%
  plyr::ldply(function(x) x$al_seq, .id = NULL) %>%
  dplyr::rename(allele = name) %>% tibble::as_tibble()
# create list
genotypes_amplisas <- list(asvs = gen, al_seq = al)

#2. compare amplisas to easyamp
# read results from easyamp
genotypes_easyamp <-
  readRDS("data/intermediate/genotypes_notnull.rds")
# reformat
genotypes_easyamp$asvs %<>%
  dplyr::select(-read_count) %>%
  dplyr::mutate(allele = as.character(allele))

#synonymize alleles
#input genotypes to merge and returns genotypes with synonymized alleles
datasets <- c("amplisas", "easyamp")
synonymized_genotypes <-
  synonymize_alleles(genotypes_amplisas,
                     genotypes_easyamp,
                     ab_names = datasets)
#3. write output genotype tables
genotypes_hierfstat <-
  datasets %>%
  lapply(function(gen_data) {
    temp1 <-
      list(asvs = synonymized_genotypes$asvs %>%
                    dplyr::filter(set == gen_data) %>%
                    dplyr::select(-set),
           al_seq = synonymized_genotypes$al_seq)
    generate_output1(temp1, paste0("amplisas-analysis/output1_", gen_data, ".csv"))
  }) %>% setNames(c("amplisas", "easyamp"))

#4. compute similarity between genotypes
# assert the genotype matrices to compare have the same loci and samples
assertthat::assert_that(all(
  names(genotypes_hierfstat$amplisas) == names(genotypes_hierfstat$easyamp)),
  msg = "Amplisas and easyamp produce genotype matrices with different loci")
assertthat::assert_that(all(
    rownames(genotypes_hierfstat$amplisas) == rownames(genotypes_hierfstat$easyamp)),
  msg = "Amplisas and easyamp produce genotype matrices with different samples")
# T/F matrix with genotypes
out_file_dif <- "amplisas-analysis/output1_differences.csv"
diff_gen <- as.data.frame(genotypes_hierfstat$amplisas == genotypes_hierfstat$easyamp)
write.csv(diff_gen, file = out_file_dif)

# compute proportion of similar genotype calls over the total of evaluated non-NA genotypes
temp2 <-
  diff_gen %>%
  dplyr::select(-pop)
same_genotypes <- sum(temp2 == T, na.rm = T)
total_nonNA <- sum(!is.na(temp2))
same_genotypes / total_nonNA
cat("\nA",
    round(same_genotypes / total_nonNA, 3),
    "of the genotype calls were silimar between EasyAmp and ampliSAS over the total of evaluated non-NA genotypes. A detailed output is written at", out_file_dif, "\n")
  

