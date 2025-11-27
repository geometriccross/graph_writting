library(conflicted)
library(tidyverse)
library(phyloseq)
library(microViz)
library(gtools)
library(dplyr)

source("./src/q2obj.r")


conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")


read_abundance <- function(df, genus_name) {
    df %>%
        filter(Genus == genus_name) %>%
        select(-FeatureID, -Genus) %>%
        colSums()
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

