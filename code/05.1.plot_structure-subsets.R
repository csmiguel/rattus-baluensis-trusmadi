###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: plot structure results
#DESCRIPTION:
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(tidyGenR)
library(starmie)
library(pophelper)
library(ggplot2)
library(gridExtra)
library(dplyr)

# pop data
lev_str <- c("Kinabalu", "Tambuyukon", "Trusmadi")
meta <-
  read.csv("data/raw/metadata.csv") |> # sort by label for structure
  mutate(population = factor(population, levels = lev_str)) |>
  arrange(population)

# subsets
subsets <- c("subset1", "subset2")

lapply(subsets, function(subset_i) {
  # list output from structure files
  dirstr <- file.path("data/intermediate/structure-subsets", subset_i)
  fp <- list.files(dirstr, pattern = "str_", full.names = T)
  qlist <- pophelper::readQStructure(fp, indlabfromfile = T)
  # reorder rows according to labels str
  old_order <- rownames(qlist[[1]])
  new_order <- order(meta$population[match(old_order, meta$sample)])
  
  # apply to all, reorder samples according to pop
  qlist_ord <- 
    lapply(qlist, function(x) {
      x[new_order,, drop = F]
    })
  
  
  # pop labels
  onelabset1 <- 
    meta[match(rownames(qlist_ord[[1]]), meta$sample)
         ,"population",drop=FALSE] |>
    dplyr::mutate(population = as.character(population))
  
  
  # tabulate
  tq <- pophelper::tabulateQ(qlist_ord)
  sq <- pophelper::summariseQ(tq)
  
  evanno <- evannoMethodStructure(data = sq,
                                  returndata = T,
                                  exportplot = F,
                                  returnplot = T,
                                  imgtype = "pdf"
  )
  pevanno <- grid.arrange(evanno$plot)
  ggsave(filename = paste0("output/evanno_structure-", subset_i, ".pdf"),
         plot = pevanno,
         width = 5, height = 5)
  
  # plot
  # align K with local function
  alk <-
    tidyGenR:::align_matrices(qlist_ord) %>%
    as.qlist()
  # set1 
  set1 <- grepl("K[2-4]", names(alk))
  p1 <-
    plotQ(alk[set1],
          imgoutput="join",
          returnplot = T,
          sharedindlab = T,
          #sortind = "all",
          exportplot = F,
          showindlab = T,
          useindlab = T,
          splab = paste0("K = ",sapply(alk[set1],ncol)),
          #indlabwithgrplab = T,
          grplab = onelabset1)
  # set 2
  set2 <- grepl("K[5-7]", names(alk))
  p2 <-
    plotQ(alk[set2],
          imgoutput="join",
          returnplot = T,
          sharedindlab = T,
          #sortind = "all",
          exportplot = F,
          showindlab = T, 
          useindlab = T,
          splab = paste0("K = ",sapply(alk[set2],ncol)),
          #indlabwithgrplab = T,
          grplab = onelabset1)
  pall <- grid.arrange(p1$plot[[1]],
                       p2$plot[[1]],
                       ncol = 2)
  ggsave(plot = pall,
         filename = paste0("output/kall-", subset_i, ".pdf"),
         width = 12,
         height = 7)
  
  # get averages for the 4 reps and save to plot
  # for each K get average cluster assignation
  Ks <- unique(stringr::str_extract(names(alk), "K[0-9]"))
  alk_m <- lapply(alk, as.matrix)
  
  # after checking the individual runs it seems one run at K3 converged to another solution
  # with less likelihood. I will remove this run for estimating the q-average.
  # select matrices
  alk_m_sel <- alk_m[names(alk_m) %in% tq$file[c(5:8, 9, 11, 12)]]
  # average
  alk_av_selection <-
    c("K2", "K3") %>%
    lapply(function(x) {
      z <- alk_m_sel[grepl(x, names(alk_m_sel))]
      starmie::averageQ(z) %>%
        as.data.frame()
    }) %>%
    as.qlist()
  
  # plot 2 to 3
  pav23 <-
    plotQ(alk_av_selection,
          imgoutput="join",
          returnplot = T,
          sharedindlab = T,
          #sortind = "all",
          exportplot = F,
          showindlab = T, 
          useindlab = T,
          splab = paste0("K = ", sapply(alk_av_selection, ncol)),
          #indlabwithgrplab = T,
          grplab = onelabset1)
  
  ggsave(plot = grid.arrange(pav23$plot[[1]]),
         filename = paste0("output/Qav_k2-3-", subset_i, ".pdf"),
         width = 6,
         height = 1.3)
  })

