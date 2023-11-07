library(tidyverse)
library(DESeq2)
library(airway)


data("airway")

des <-
  airway |>
  # convert to DataSet for DESeq
  # analyse expression by cell line and dex treatment
  DESeqDataSet(design = ~ cell + dex) |>
  # Run DESeq2 normalization and regression models
  DESeq()

res <- results(des, tidy = TRUE)

res |>
  ggplot(aes(x = log2FoldChange, 
             y = - log10(padj))) +
  geom_point() +
  theme_bw(18)

ggsave('volcano.jpeg', width = 6, height = 5)


################################################################################

my.analysis <-
  res |>
  filter(padj <= 0.05, abs(log2FoldChange) >= 1) |>
  left_join(airway |> rowData() |> as_tibble(), c('row' = 'gene_id')) |>
  select(gene_id = row, gene_name, log2FoldChange, FDR = padj) |>
  arrange(FDR)

library(openxlsx)

my.analysis |>
  write.xlsx('my-analysis.xlsx')


# my.analysis |>
#   head() |>
#   knitr::kable('markdown')
# |gene_id         |gene_name | log2FoldChange| FDR|
# |:---------------|:---------|--------------:|---:|
# |ENSG00000152583 |SPARCL1   |      -4.574919|   < 1e-100|
# |ENSG00000165995 |CACNB2    |      -3.291062|   < 1e-100|
# |ENSG00000120129 |DUSP1     |      -2.947810|   < 1e-100|
# |ENSG00000101347 |SAMHD1    |      -3.766995|   < 1e-100|
# |ENSG00000189221 |MAOA      |      -3.353580|   < 1e-100|
# |ENSG00000211445 |GPX3      |      -3.730403|   < 1e-100|
