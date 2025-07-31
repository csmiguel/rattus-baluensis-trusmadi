###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com
# May 2024
###.............................................................................
#GOAL: compute and plot multilocus heterozygosity
#PROJECT: rbaluensis-trusmadi
###.............................................................................
library(inbreedR)
library(tidyverse)
library(ggrepel)
library(tidyGenR)

# genotypes
gen <-
  readRDS("output/genotypes_rbal_poly.rds") |>
  tidyGenR::gen_tidy2wide() |>
  tibble::column_to_rownames("sample")

# metadata
meta <-
  read.csv("data/raw/metadata.csv")
# generate input for inbreedR
gen_inbreedR <-
  #select(gen, -sample) |>
  tidyr::separate_wider_delim(gen,
                              cols = everything(),
                              delim = "/",
                              names_sep  = ".",
                              too_few = "align_start") %>%
  inbreedR::convert_raw()
# check data
inbreedR::check_data(gen_inbreedR)
# multilocus heterocigosity with metadata
mlh <-
  inbreedR::MLH(gen_inbreedR) |>
  setNames(rownames(gen)) %>%
  {data.frame(sample = names(.),
              pop = as.factor(meta$population[match(names(.), meta$sample)]),
              elevation = as.numeric(meta$elevation[match(names(.), meta$sample)]),
              missing_calls = apply(gen_inbreedR, 1, function(x) sum(is.na(x))),
              mlh = .)} |>
  arrange(pop, sample, mlh)

# plot mlh
p1 <-
  ggplot(mlh, aes(x = pop, y = mlh, label = sample, color = elevation)) +
  geom_boxplot(color = "grey60", outliers = F) +
  geom_jitter(aes(size = missing_calls), width = .2) +
  ggrepel::geom_text_repel(size = 2) +
  scale_color_gradient(low = "green", high = "brown") +
  ylab("inbreedR: multilocus heterozygosity") +
  xlab(NULL) +
  theme_classic()
ggsave(plot = p1,
       filename = "output/mlh.pdf",
       width = 5, height = 4)

# model: are there differences in individual heterozygosity between mountains
mlh2 <-
  mlh %>%
  dplyr::mutate(elevation = elevation - 2000) %>%
  dplyr::select(-sample, -missing_calls)
m1 <- lm(mlh ~ elevation * pop, data = mlh2)
m2_noInt <- lm(mlh ~ elevation + pop, data = mlh2)
m3_only_pop <- lm(mlh ~ pop, data = mlh2)
m4_only_elev <- lm(mlh ~ elevation, data = mlh2)
m5 <- lm(mlh ~ pop + elevation:pop, data = mlh2)

sink("output/mlh_model.txt")
summary(m1)
anova(m1, reduced_model)
sink()

# plot with elevation
p2 <-
  ggplot(data = mlh, aes(x = elevation, y = mlh, color = pop)) +
  stat_smooth(method = "lm", linewidth = 1.5) +
  geom_point() +
  ylab("Multilocus heterozigosity (MLH)") +
  xlab("Elevation (m)") +
  #scale_color_manual(name = "Population") +
  theme_classic() +
  theme(legend.title = element_blank())
ggsave(plot = p2,
       filename = "output/mlh2.pdf",
       width = 6, height = 4)
