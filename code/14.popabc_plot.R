# genetic diversity of pops
###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# July 2025
###.............................................................................
#GOAL: regression step and plot results from POPABC
#PROJECT: rbaluensis-trusmadi
###.............................................................................

library(stringr)
library(locfit)
library(dplyr)

#import Mark Beaumont's scripts. move the scripts to the working directories
source("code/functions/make_pd2.r")
source("code/functions/loc2plot_d.r")
source("code/functions/plot_params.R")

# 1. checking the simulations with a PCA
# 2. regression step to weight the simulations according to their distance
#     to the real data and get HPD.

### read files ####
# path to files
p_txt <- "data/intermediate/popabc/run1/run1.txt"
p_rej <- "data/intermediate/popabc/run1/run1.rej"
p_trg <- "data/intermediate/popabc/run1/run1.trg"

## TXT ## read summary txt with parameters and stats used
txt <- readLines(p_txt)
#   extract names of parameters
names_params <-
  txt[grep("parameters", txt)] |>
  str_extract("\\[(.*)\\]", group = 1) |>
  str_split_1("\\|")
#   extract names of pop stats
names_stats <-
  txt[grep("summstats", txt)] |>
  str_extract("\\[(.*)\\]", group = 1) |>
  str_split_1("\\|")
#   number of parameters
n_params <- length(names_params)
n_stats <- length(names_stats)

### REJECTION ### read rejection file
rej <-
  read.table(p_rej)
names(rej) <- c(names_params, names_stats)
rej <-
  filter(rej, Ne1 < NeA1) # filter cases where Kinabalu pop are larger
rej_stats <- rej[ , (n_params + 1):(n_params + n_stats)]
rej_params <- rej[ , 1:n_params]
 
### TRG - read the real data genetic statistics
target <- read.table(p_trg)
names(target) <- names(rej_stats)
target <- as.matrix(target)

### PRIORs
p_priors <- "data/intermediate/popabc/run1/run1.pri"
priors <- as.matrix(read.table(p_priors,
                               col.names = names_params))


#### 1. Checking the simulation with a PCA #### 
# Principal component analysis to check over the rejection file
pca_rej <- prcomp(rej_stats, center = T, scale = T)
real_data_PC <- predict(pca_rej, newdata = target)

pdf("output/popabc-pca.pdf",
    width = 4,
    height = 4)

plot(pca_rej$x, pch = 19, col = rgb(0, 0, 0, alpha = 0.07))

#plot the real data genetic statistics
points(real_data_PC[1],
       real_data_PC[2],
       col="red", cex=2, pch=16)
dev.off()

#### 2. Regression step ####

# Plot line for splitting time and save it in a .eps file
time_range <- c(1, 30000)
pdf("output/popabc_time.pdf",
    width = 5,
    height = 4)
hpd_time <-
  plot_params(c("t1", "t2"),
            rej = rej,
            param_range = time_range,
            xlab = "Years ago", lw = 2,
            prob = 0.05)
dev.off()

# Plot line for splitting time and save it in a .eps file
ne_range <- c(100, 100000)
# current pops sizes
pdf("output/popabc_popsize.pdf",
    width = 5,
    height = 4)
plot_params(c("Ne2", "Ne1", "Ne3"),
            rej = rej, lw = 2,
            param_range = ne_range,
            xlab = "Effective population size")
dev.off()

# ancestral
pdf("output/popabc_popsize-ancestral.pdf",
    width = 5,
    height = 4)
plot_params(c("NeA2", "NeA1"),
            rej = rej, lw = 2,
            param_range = ne_range,
            xlab = "Effective population size")
dev.off()

# all
pdf("output/popabc_popsize-ancestral.pdf",
    width = 6,
    height = 4)
hpd_ne <-
  plot_params(c("Ne2", "Ne1", "Ne3", "NeA2", "NeA1"),
            rej = rej, lw = 2,
            param_range = ne_range,
            xlab = "Effective population size")
dev.off()

# export HPD stats
hpd_stats <-
  do.call(c(hpd_ne, hpd_time), what = "rbind") |>
  as.data.frame() |>
  setNames(c("mode", "lowerHPD", "higherHPD"))
write.table(hpd_stats,
            file = "output/popabc-mode.txt")

