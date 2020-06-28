###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: dada2
#PROJECT: amplicon-genotyping
###.............................................................................

library(dada2)
library(dplyr)

# filters
# allele balance
ab <- 0.3
cov <- 11
ploidy <- 2
all_samples <- readLines("data/intermediate/samples-list")
loci <- readLines("data/raw/loci")

clean_lowfrequencyvariants <- function(table_ab) {
  #table_ab is a matrix, array or dataframe with cols being alleles, rows being samples and cells being variant count.
  h <-
    apply(table_ab, 1, function(x){
    h <- max(x)
    x[x/h < ab | x < cov] <- 0
    x
  }) %>% as.data.frame() %>% setNames(colnames(table_ab))
  #when there is only one allele
  if(ncol(table_ab) > 1) h <- t(h)
  #after initial cleaning, some alleles are spurious, and must be removed
  #for that reaseon. I will remove cols which are all 0
  h[, apply(h, 2, sum) > 0] %>%
    as.data.frame(row.names = rownames(table_ab)) %>%
      {setNames(object = ., colnames(table_ab)[1:ncol(.)])}
}

#for all loci
genotypes <-
  loci[24:30] %>%
  lapply(function(locus) {
    fnFs <- sort(list.files("data/intermediate/",
                            pattern = paste0("*", locus, ".*1.fastq.gz"),
                            full.names = TRUE))
    fnRs <- sort(list.files("data/intermediate",
                            pattern = paste0("*", locus, ".*2.fastq.gz"),
                            full.names = TRUE))
    sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
    #if there are no files with the locus (ie, not amplified), then return 00 genotypes
    if(length(fnRs) == 0 | length(fnFs) == 0){
      hh <- rep("00", length(all_samples))
      names(hh) <- all_samples
      hh
    }else{

    #create path for sequences that will be filtered downstream
    filt_path <- "data/intermediate/filtered"
    if(!file_test("-d", filt_path)) dir.create(filt_path)

    filtFs <- file.path(filt_path, paste0(sample.names, locus, "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_path, paste0(sample.names, locus, "_R_filt.fastq.gz"))

    #filter. There cannot be ambiguous positions
    out <- filterAndTrim(fnFs, filtFs,
                         fnRs, filtRs,
                         maxN = 0, maxEE = c(3, 6),
                         compress=TRUE, multithread=TRUE)
    #if filtering returns no reads, then assign 00 to each genotype
    if(sum(out[, 2]) == 0){
      hh <- rep("00", length(all_samples))
      names(hh) <- all_samples
      hh
      }else{
    #error matrices

    errF <- dada2::learnErrors(filtFs, multithread = TRUE)
    errR <- dada2::learnErrors(filtRs, multithread = TRUE)

    #dereplication
    derepFs <-
      dada2::derepFastq(filtFs, verbose = TRUE) %>% setNames(sample.names)
    derepRs <-
      dada2::derepFastq(filtRs, verbose = TRUE) %>% setNames(sample.names)

    #infer ASV for pooled samples
    dadaFs <- dada2::dada(derepFs, err = errF, multithread = TRUE, pool = T)
    dadaRs <- dada2::dada(derepRs, err = errR, multithread = TRUE, pool = T)

    #Merge paired reads. Spurious sequence variants are further reduced by merging overlapping reads.
    mergers <-
      dada2::mergePairs(dadaFs, derepFs,
                        dadaRs, derepRs,
                        verbose = TRUE, minOverlap = 10, maxMismatch = 0)
    #if reads do not merge, just concanetate them
    if(mergers %>% sapply(function(x) nrow(x) == 0) %>% all()){
      mergers <-
        dada2::mergePairs(dadaFs, derepFs,
                          dadaRs, derepRs,
                          verbose = TRUE, justConcatenate = T)
    }

    #Construct sequence table from the provided list of samples:
    seqtab <- makeSequenceTable(mergers)
    #name alleles
    colnames(seqtab) <- letters[1:ncol(seqtab)]
    #check table
    seqtab
    #we can see some very low frequency alleles in some samples. I will remove those
    #Remove bimeras:
    seqtabNoC <- dada2::removeBimeraDenovo(seqtab)

    #remove low frequency variants
    hh <- clean_lowfrequencyvariants(seqtabNoC)

    #convert to 2 column dataframe
    hh %>%
      apply(1, function(x){
        y <- names(x)[x > 0]
        if(length(y) == 0) "00"
        else if(length(y) == 1) paste0(y, y, collapse = "")
        else paste0(y, collapse = "")
      }) %>% #I will convert it to a data frame and do a left join with all samples in case there is anyone missing
      as.data.frame() %>%
      tibble::rownames_to_column(var = "id") %>%
      {dplyr::left_join(x = data.frame(id = all_samples), y = ., by = "id")} %>%
      tibble::column_to_rownames(var = "id")
    }
  }
    }) %>%
  {do.call(what = cbind, args = .)}
