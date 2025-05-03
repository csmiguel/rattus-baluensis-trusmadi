###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: haplotype networks
#PROJECT: rbaluensis-trusmadi
###.............................................................................

library(geneHapR)
library(dplyr)

gen <-
  readRDS("output/genotypes_rbal_poly.rds")
msa <-
  readRDS("output/msa.rds")

# metadata for pie colors
meta <-
  dplyr::select(gen, population, sample, locus)

# to all, create hap network and plot
hap_plots <-
  names(msa) |>
  lapply(function(x) {
    cat(paste("\nbegin", x, "\n"))
    # add a constant base as 1st column of MSA because seqs2hap discard indels at 1st column.
    cmsa <-
      cbind(matrix(rep("A", length(msa[[x]])), ncol = 1),
            as.matrix(msa[[x]]))  %>%
      apply(1, function(x) paste(x, collapse = "")) %>%
      Biostrings::DNAStringSet()
    
    h <- seqs2hap(cmsa,
                  hapPrefix = "H",
                  hetero_remove = TRUE,
                  na_drop = TRUE,
                  maxGapsPerSeq = .25,
                  replaceMultiAllele = TRUE) %>%
      hap_summary()
    accinfo1 <- data.frame(filter(meta, locus == x),
                           row.names = names(msa[[x]]))
    # rownames of accINFO must match names in msa_poli
    hapN <- get_hapNet(h,
                       AccINFO = accinfo1,
                       groupName = "population")
    plotHapNet(hapN,
               size = "freq",                   # circle size
               cex = 0.8,                       # size of hap symbol
               col.link = 2,                    # link colors
               link.width = 2,                  # link widths
               show.mutation = 2,               # mutation types one of c(0,1,2,3)
               legend = T,
               main = x)
    cat(paste("\nend", x, "\n"))
    p <- recordPlot()
    pp <- cowplot::plot_grid(p)
    ggsave(filename = file.path("output", "geneHapR",
                                paste0(x,"_hap_network.pdf")), plot = pp)
    return(p)
  })

# arrange plots in one page
 library(cowplot)
 n <- length(hap_plots)
 nCol <- floor(sqrt(n))
# 
 pp <- cowplot::plot_grid(plotlist = hap_plots, ncol = nCol)
 
 ggsave(filename = "output/geneHapR/all_loci.svg", width = 40, height = 50, units = "cm")
