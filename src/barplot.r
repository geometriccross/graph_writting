library(tidyverse)
library(phyloseq)
library(microViz)
library(qiime2R)
library(gtools)

source("./src/q2obj.r")


write_plot <- function(ps, filename, tax_order = sum, width = 16, height = 9, dpi = 600) {
    n_taxa_value <- if (is.function(tax_order)) 41 else length(tax_order)

    color_pal <- distinct_palette(
        n = n_taxa_value, pal = "brewerPlus", add = "gray"
    )

    plot_obj <- comp_barplot(
        ps,
        tax_level = "Genus",
        sample_order = natural_sorted,
        label = "RawID",
        n_taxa = n_taxa_value, # brewerPlusは最大個まで色分けできる
        tax_order = tax_order,
        palette = color_pal,
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

# 細胞内共生細菌のみを表示(それ以外はOtherとして集約)
intracellular_taxa <- c(
    "Bartonella",
    "Coxiella",
    "Rickettsia",
    "Anaplasma",
    "Ehrlichia",
)

write_plot(
    "./barplot_intracellular.svg",
    tax_order = c(intracellular_taxa, "Other")
)

symbiotic_taxa <- c(
    "Wolbachia",
    "Spiroplasma",
    "Rickettsiella",
    "Arsenophonus",
    "Cardinium",
    "Hamiltonella",
    "Regiella",
    "Sodalis",
    "Serratia",
    "Buchnera",
    "Blochmannia",
    "Baumannia",
    "Sulcia",
    "Nasuia"
)
