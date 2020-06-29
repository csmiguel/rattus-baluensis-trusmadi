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

#load parameters for genotype calling
source("code/parameters/parameters.r")
#load function to cleas spurious ASVs
source("code/functions/remove-spurious-ASVs.r")

#samples
all_samples <- readLines("data/intermediate/samples-list")
#loci
loci <- readLines("data/raw/loci")
loci_set <- loci[1:2]
#for all loci
genotypes <-
  loci_set %>%
  lapply(function(locus) {
    fnFs <- sort(list.files("data/intermediate",
                            pattern = paste0("*", locus, ".*1.fastq.gz"),
                            full.names = TRUE))
    fnRs <- sort(list.files("data/intermediate",
                            pattern = paste0("*", locus, ".*2.fastq.gz"),
                            full.names = TRUE))
  #samples from the fastq files
    sample.names <- sapply(strsplit(basename(fnFs), "\\."), `[`, 1)
    #if there are no files with the locus (ie, not amplified),
    if (length(fnRs) == 0 | length(fnFs) == 0) {
      cat("no files for", locus)
    }else{
    #create path for sequences that will be filtered downstream
    filt_path <- "data/intermediate/filtered"
    if (!file_test("-d", filt_path)) dir.create(filt_path)

    filtFs <- file.path(filt_path, paste0(sample.names, locus,
      "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_path, paste0(sample.names, locus,
      "_R_filt.fastq.gz"))

    #filter. There cannot be ambiguous positions
    out <- dada2::filterAndTrim(fnFs, filtFs,
                         fnRs, filtRs,
                         maxN = 0, maxEE = c(3, 6),
                         compress = TRUE, multithread = TRUE)
    #if filtering returns no reads, then assign 00 to each genotype
    if (sum(out[, 2]) == 0) {
      cat("filtering removed all reads for", locus)
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

    #Merge paired reads. Spurious sequence variants are further reduced
    #by merging overlapping reads.
    mergers <-
      dada2::mergePairs(dadaFs, derepFs,
                        dadaRs, derepRs,
                        verbose = TRUE, minOverlap = 10, maxMismatch = 0)
    #if reads do not merge, just concanetate them
    if (mergers %>% sapply(function(x) nrow(x) == 0) %>% all()){
      mergers <-
        dada2::mergePairs(dadaFs, derepFs,
                          dadaRs, derepRs,
                          verbose = TRUE, justConcatenate = T)
    }

    #Construct sequence table from the provided list of samples:
    seqtab <- dada2::makeSequenceTable(mergers) %>%
            dada2::removeBimeraDenovo() #Remove bimeras
    #remove spurious ASVs
    hh <- clean_lowfrequencyvariants(seqtab)
    #allele names
    al_names <- letters[1:ncol(hh)]
    #create data frame with allele sequences
    al_seq <- data.frame(name = al_names,
                     seq = names(hh))
    #name alleles
    colnames(hh) <- al_names
    #convert to 2 column dataframe
    geno <-
      hh %>%
        apply(1, function(x) {
          y <- names(x)[x > 0]
          if (length(y) == 0) "00"
          else if (length(y) == 1) paste0(y, y, collapse = "")
          else paste0(y, collapse = "")
        }) %>%
        #I will convert it to a data frame and do a left join with all samples
        #in case there is anyone missing
        as.data.frame() %>%
        tibble::rownames_to_column(var = "id") %>%
        {dplyr::left_join(x = data.frame(id = all_samples), y = ., by = "id")} %>%
        tibble::column_to_rownames(var = "id")
    geno[is.na(geno), 1] <- "00"
    list(genotypes = geno, allele_seq = al_seq)
    }
  }
    }) %>%
    setNames(loci_set)
#create data frame with genotypes
genotypes <-
  genotyping %>%
    sapply(function(x){
      x$genotypes
      }) %>%
      do.call(what = cbind) %>%
      {row.names(.) <- sample.names; .} %>%
      as.data.frame() %>%
      setNames(names(genotyping))
#create data frame with allele sequences
al_seqs <-
  genotyping %>%
    lapply(function(x) x$allele_seq) %>%
    reshape2::melt() %>%
    dplyr::rename(allele = name,
                  sequence = seq,
                  locus = L1) %>%
    dplyr::select(locus, allele, sequence) %>%
    dplyr::arrange(locus, allele)

#write alleles to fasta
write.fasta(sequences = as.list(al_seqs$sequence),
            names = paste(al_seqs$locus, al_seqs$allele, sep = "_"),
            file.out = "output/alleles.fasta")

saveRDS(genotypes, file = "data/intermediate/genotypes.rds")
saveRDS(al_seqs, file = "data/intermediate/allele_sequences.rds")
