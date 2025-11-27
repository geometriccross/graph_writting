library(tidyverse)
library(phyloseq)
library(microViz)
library(qiime2R)
library(gtools)


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


natural_sorted <- phyloseq_object %>%
    sample_data() %>%
    as.data.frame() %>%
    {
        rownames(.)[mixedorder(.$RawID)] # RawIDを対象とする
    }


composition_barplot <- phyloseq_object %>%
    comp_barplot(
        tax_level = "Genus",
        sample_order = natural_sorted,
        label = "RawID",
        n_taxa = 41, # brewerPlusは最大個まで色分けできる
        palette = distinct_palette(n = 41, pal = "brewerPlus", add = "gray"),
        bar_outline_colour = NA,
        bar_width = 0.9,
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


ggsave("barplot.svg", plot = composition_barplot, width = 16, height = 9, dpi = 600)
