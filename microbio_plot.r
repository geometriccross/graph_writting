library(tidyverse)
library(phyloseq)
library(microViz)
library(qiime2R)


# QIIME2メタデータの読み込みとクリーニング（型定義行の削除とID列の修正）
clean_sample_data <- read_tsv("./q2objs/metadata.tsv", comment = "") %>%
    rename_with(~"SampleID", .cols = 1) %>%
    filter(!str_detect(SampleID, "^#q2:types|^categorical|^numeric")) %>%
    column_to_rownames("SampleID") %>%
    sample_data()


# Feature table等の読み込み、分類名の修正、メタデータの結合
phyloseq_object <- qza_to_phyloseq(
    features = "./q2objs/common_biology_free_table.qza",
    taxonomy = "./q2objs/common_biology_free_classification.qza",
    tree     = "./q2objs/biology_free_rooted_tree.qza"
) %>%
    tax_fix() %>%
    tax_fix(unknowns = c("uncultured")) %>%
    merge_phyloseq(clean_sample_data)
