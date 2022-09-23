#!/usr/bin/env Rscript

library(tidyverse)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)

vcf <- args[1]
metadata <- args[2]
plot_name <- args[3]

meta <- read_tsv(metadata)

variants <-
  read_tsv(vcf, comment = '##') %>%
  pivot_longer(-(1:9), names_to = 'accession') %>%
  mutate(GT = str_extract(value, '^([^:]+)') %>% str_count('1')) %>%
  select(accession, POS, ALT, GT) %>%
  unite('VID', POS, ALT) %>%
  inner_join(select(meta, accession, date),
             by = 'accession') %>%
  group_by(VID) %>%
  filter(sum(GT == 2) > 1,
         sum(GT == 2) > sum(GT == 1)) %>%
  group_by(accession) %>%
  filter(sum(GT == 1) < 3) %>%
  ungroup() %>%
  mutate(GT = case_when(GT == 0 ~ 0, GT == 2 ~ 1))

variant_matrix <-
  variants %>%
  select(-accession) %>%
  arrange(date) %>%
  pivot_wider(names_from = date, values_from = GT) %>%
  as.data.frame() %>%
  column_to_rownames('VID') %>%
  as.matrix()


png(plot_name, width=6, height=4, units="in", res=1200)

pheatmap(variant_matrix,
         color = viridisLite::cividis(2),
         cluster_cols = F,
         show_rownames = F,
         clustering_distance_rows  = 'manhattan',
         clustering_method = 'average',
         main = 'Victorian SARS-CoV-2 Varaints 2020',
         border_color = NA,
         legend_breaks = c(0,1),
         angle_col = 45)

dev.off()
