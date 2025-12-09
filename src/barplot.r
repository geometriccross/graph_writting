library(tidyverse)
library(phyloseq)
library(microViz)
library(qiime2R)
library(gtools)

source("./src/q2obj.r")


write_plot <- function(ps, filename, width = 16, height = 9, dpi = 600) {
    plot_obj <- comp_barplot(
        ps,
        tax_level = "Genus",
        sample_order = natural_sorted,
        label = "RawID",
        n_taxa = 41, # brewerPlusは最大個まで色分けできる
        palette = distinct_palette(n = 41, pal = "brewerPlus", add = "gray"),
        bar_outline_colour = NA,
        bar_width = 0.9,
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    ggsave(filename, plot = plot_obj, width = width, height = height, dpi = dpi)
}


phyloseq_object <- load_q2obj() %>% tax_tweak()

natural_sorted <- phyloseq_object %>%
    sample_data() %>%
    as.data.frame() %>%
    {
        rownames(.)[mixedorder(.$RawID)] # RawIDを対象とする
    }


# 除外する正規表現パターンのリスト
exclude_patterns <- c(
    "Family", # 属だけで見たいため
    "[A-Z]{3}-\\d+", # UCG-003のような不明瞭な属
    "_group$",
    "1174-901-12"
) %>% paste(collapse = "|")

phyloseq_filtered <- phyloseq_object %>%
    subset_taxa(!grepl(exclude_patterns, Genus))


write_plot(phyloseq_filtered, "./barplot_genus_filtered.svg")
