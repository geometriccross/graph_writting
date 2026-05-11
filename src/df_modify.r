library(tidyverse)
library(phyloseq)

# phyloseqオブジェクトからOTU + Taxonomyを結合したdata.frameを生成
ps_to_df <- function(ps) {
    otu <- ps %>%
        otu_table() %>%
        as.data.frame() %>%
        rownames_to_column("FeatureID")

    tax <- ps %>%
        tax_table() %>%
        as.data.frame() %>%
        select("Genus") %>%
        rownames_to_column("FeatureID")

    left_join(otu, tax, by = "FeatureID")
}

# 指定Genus（または全体）のサンプルごとAbundanceを取得
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

# 2つのAbundanceの比を計算
relative_abundance <- function(a, b) {
    a %>%
        left_join(b, by = "SampleID", suffix = c("_a", "_b")) %>%
        mutate(Ratio = Abundance_a / Abundance_b) %>%
        select(SampleID, Ratio)
}

# 指定Genusの相対存在比を取得
get_ratio <- function(df, key = NULL) {
    relative_abundance(
        read_abundance(df, key),
        read_abundance(df, "*")
    ) %>%
        mutate(Genus = key)
}

# 比率データから統計情報を計算 (n, mean, sd, se, detection_rate)
calculate_stats <- function(ratio_df) {
    ratio_df %>%
        group_by(Genus) %>%
        summarise(
            n = n(),
            mean_ratio = mean(Ratio, na.rm = TRUE),
            sd_ratio = sd(Ratio, na.rm = TRUE),
            se_ratio = sd_ratio / sqrt(n),
            detection_rate = sum(Ratio > 0, na.rm = TRUE) / n,
            .groups = "drop"
        )
}
