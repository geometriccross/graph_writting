library(conflicted)
library(tidyverse)
library(phyloseq)
library(microViz)
library(gtools)
library(dplyr)

source("./src/q2obj.r")


conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")


read_abundance <- function(df, key = NULL) {
    result <- df

    if (!(is.null(key) || key == "*")) {
        result <- result %>% filter(Genus == key)
    }

    result %>%
        select(-FeatureID, -Genus) %>%
        colSums() %>%
        enframe(name = "SampleID", value = "Abundance")
}


relative_abundance <- function(a, b) {
    a %>%
        left_join(b, by = "SampleID", suffix = c("_a", "_b")) %>%
        mutate(Ratio = Abundance_a / Abundance_b) %>%
        select(SampleID, Ratio)
}

get_ratio <- function(df, key = NULL) {
    relative_abundance(
        read_abundance(df, key),
        read_abundance(df, "*")
    ) %>%
        mutate(Genus = key)
}


# 複数の菌属のratio dataframeを統合する関数
combine_ratio_dfs <- function(ratio_list, genus_names) {
    combined <- map2_dfr(
        ratio_list,
        genus_names,
        ~ .x %>% mutate(Genus = .y)
    )
    return(combined)
}


ps <- load_q2obj()


# OTUテーブルを取得
otu <- ps %>%
    otu_table() %>%
    as.data.frame() %>%
    rownames_to_column("FeatureID")

# 分類学的情報を取得
tax <- ps %>%
    tax_table() %>%
    as.data.frame() %>%
    select("Genus") %>%
    rownames_to_column("FeatureID")

# 結合して表示
df <- otu %>%
    left_join(tax, by = "FeatureID")

ratio_df <- get_ratio(df, key = "Bartonella")
print(ratio_df)
