# genetic diversity of pops
###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: write population genetic statistics.
# private alleles, histogram with alleles per locus, no of samples.
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(tidyverse)
library(scales)
library(tidyGenR)

# 1. histogram plot of number of alleles per locus ####
n_hap <-
  readRDS("output/genotypes_rbal.rds") %>%
  dplyr::filter(locus != "fetub") %>%
  dplyr::group_by(locus) %>%
  dplyr::summarise(n_loc = length(unique(allele)))

p_hap <-
  ggplot(n_hap, aes(x = n_loc)) +
  geom_histogram(binwidth=1, colour="black", fill="grey", position="identity") +
  scale_x_continuous(breaks= breaks_pretty()) +
  scale_y_continuous(breaks= breaks_pretty()) +
  xlab("Number of alleles") +
  ylab("Count") +
  annotate("text", x = 7, y = 5, size = 3,
           label = paste(sum(n_hap$n_loc == 1), "monomorphic loci")) +
  annotate("text", x = 7, y = 4, size = 3, 
           label = paste(sum(n_hap$n_loc > 1), "polymorphic loci")) +
  theme_classic()
ggsave("output/haplotype_number.pdf",
       plot = p_hap, width = 3, height = 3)

### 2. private alleles per population ####
fun_priv_al <- function(data) {
  # 1..3 if private for trusmadi, tambuyukon, kinabalu, respectively.
  pops <- c("trusmadi", "tambuyukon", "kinabalu")
  select(data, all_of(pops)) %>% 
    apply(1, function(x) {
      x <- as.numeric(x)
      if(identical(x, c(1L, 0, 0)))
        pops[1]
      else if (identical(x, c(0, 1L, 0)))
        pops[2]
      else if (identical(x, c(0, 0, 1L)))
        pops[3]
      else
        "non_private"
    })
}
priv_al <-
  readRDS("output/genotypes_rbal.rds") %>%
  dplyr::filter(!locus %in% c("fetub")) %>%
  remove_monomorphic() %>%
  pivot_wider(id_cols = c(locus, allele),
              names_from = population,
              values_from = sequence,
              values_fn = function(x) length(unique(x)),
              values_fill = 0) %>%
  mutate(priv = fun_priv_al(data = .)) %>%
  pull(priv) %>%
  table() %>%
  as.data.frame() %>%
  setNames(c("population", "n"))

# 3. n samples per pop
nsamples <-
  readRDS("output/genotypes_rbal.rds") %>%
  dplyr::select(population, sample) %>%
  plyr::daply(~population, function(x) {
    length(unique(x$sample))
  }) %>%
  as.data.frame() %>%
  setNames("n_samples")
  
sink(file = "output/genetic_diversity.txt")
cat("Private alleles: \n")
priv_al
cat("\n\nNumber of samples:\n")
nsamples
cat("\nNumber of haplotypes:\n")
as.data.frame(n_hap)
sink()
      
