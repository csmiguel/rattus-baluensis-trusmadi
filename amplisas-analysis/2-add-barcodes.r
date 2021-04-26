# add barcodes to raw sequences
library(dplyr)
library(ShortRead)
library(Biostrings)
library(DNABarcodes)

# read samples
samples <-
  read.table("data/intermediate/samples-list", col.names = "sample") %>%
  dplyr::pull(sample)
#generate barcodes
nnucl <- 6
barcod <- create.dnabarcodes(nnucl)
assertthat::assert_that(length(barcod) >length(samples),
                        msg = "The number of barcodes is less than the number of samples. Increase the n paramter")
barcodes <-
  barcod[seq_along(samples)] %>%
  {data.frame(sample = samples, barcode = .)}
# quality associated to barcode
highestqual <- rep("I", nnucl) %>% paste0(collapse = "")

# read fastq demultiplexed by individual
c("1.fastq.gz", "2.fastq.gz") %>%
  lapply(function(fr) {
    samples %>%
      lapply(function(sample) {
        path2reads <- list.files("data/raw", pattern = paste0(sample, ".", fr), full.names = T)
        # read fastq file
        h <- ShortRead::readFastq(path2reads)
        barcode1 <- barcodes$barcode[barcodes$sample == sample]
        # create new reads with barcodes
        ShortReadQ(
          quality = Biostrings::quality(h)@quality %>% #add vector with highest quality to 5' end.
            as.character() %>%
            {paste0(highestqual, .)} %>%
            BStringSet(),
          sread = ShortRead::sread(h) %>% #add barcode to 5' end
            as.character() %>%
            {paste0(barcode1, .)} %>%
            Biostrings::DNAStringSet(),
          id = h@id #same fastq headers
        ) %>%
          ShortRead::writeFastq(file = file.path("amplisas-analysis", fr),
                                mode = "a") #append to existing file
      })
    })

# export barcode data
barcodes %>%
  dplyr::mutate(barcode_r = barcode) %>%
  dplyr::rename(">sample" = sample,
                barcode_f = barcode) %>%
write.table(file = "amplisas-analysis/amplicon-data.csv",
          quote = F,
          row.names = F,
          sep = ",",
          append = T)
