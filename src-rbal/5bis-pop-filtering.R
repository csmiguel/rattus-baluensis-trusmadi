###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: reformat and filter genotypes from population genetics
#PROJECT: amplicon-genotyping
###.............................................................................
library(dplyr)
library(pegas)

#load genotypes
gen_genind <- readRDS("output/genind-trusmadikinabalu.rds")

#load filtering parameters
source("src/parameters/pop-filtering.r")
#load filtering functions
source("src/functions/pop-filt.r")
#read ploidy
# I do a trick to prevent populating the working env with extra variables
get_ploidy <- function() {
  source("src/parameters/parameters.r", local = T)
  ploidy
}
ploidy <- get_ploidy()

sink("output/log.txt", append = T)
cat(
  as.character(Sys.time()),
  "\nLOG on src/5-pop-filtering.R\n")
#start filtering
library(dplyr)
library(pegas)
library(ggplot2)
#load genotypes
gen_genind <- readRDS("output/genind-trusmadikinabalu.rds")

#load filtering parameters
source("src/parameters/pop-filtering.r")
#load filtering functions
source("src/functions/pop-filt.r")
#read ploidy
# I do a trick to prevent populating the working env with extra variables
get_ploidy <- function() {
  source("src/parameters/parameters.r", local = T)
  ploidy
}
ploidy <- get_ploidy()

#start filtering

temp <- gen_genind %>%
  # filter out individuals according to their missing data
  filter_ind_missingness() %>%
  # filter out loci according to their missing data
  missingness_locus() %>%
  # remove monorphic loci
  remove_monomorphic()

#imputate data and estimate heterozyosity and hwe tests
hwe <-
  imputation_hw(temp) %>% adegenet::seppop() %>%
 lapply(function(x) {
   div <- adegenet::summary(x)
   z <- data.frame(Hobs = div$Hobs, Hexp = div$Hexp)
   h <-
     hw.test(x) %>%
     as.data.frame() %>%
     dplyr::select(Pr.exact)
   list(het = z, hwtest = h)
   }) %>%
 setNames(levels(pop(temp)))
#dataframe with heterozyosity
div <- names(hwe) %>%
  lapply(function(popul) {
    x <- hwe[[popul]]
    x$het %>%
      tibble::rownames_to_column("locusi") %>%
      dplyr::mutate(pop = popul)
  }) %>%
  do.call(what = rbind) %>% tibble::as_tibble()
#probabilies from hwe
probs <- lapply(hwe, function(x) x$hwtest) %>%
  do.call(what = cbind) %>%
  setNames(levels(pop(temp))) %>%
  tibble::rownames_to_column("locusi") %>%
  reshape2::melt(id = "locusi",
    value.name = "prob.not.hwe",
    variable.name = "pop") %>% tibble::as_tibble()
div_plot <-
  dplyr::left_join(div, probs,
    by = c("locusi", "pop"))

ggplot(div_plot, aes(x = Hobs, y = Hexp, label = locusi)) +
    geom_point(color = ifelse(div_plot$prob.not.hwe < 0.01, "red", "black")) +
    ggrepel::geom_text_repel() +
    facet_grid(~pop) +
    theme_classic()
ggsave(width = 10, height = 4, filename = "output/hwe_plots.pdf")


#export to structure
hierfstat::genind2hierfstat(pop_filtered_genotypes,
                            pop = adegenet::indNames(pop_filtered_genotypes)) %>%
  hierfstat::write.struct(fname = "output/output5-structure.str",
                          ilab = indNames(pop_filtered_genotypes),
                          MARKERNAMES = T,
                          pop = pop(pop_filtered_genotypes))
#save filtered genotypes
saveRDS(pop_filtered_genotypes, "data/output/output6-pop-genotypes-genind.rds")
