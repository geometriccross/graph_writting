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


# 統計情報を計算する関数（平均、標準偏差、標準誤差、サンプル数）
calculate_stats <- function(ratio_df) {
    ratio_df %>%
        group_by(Genus) %>%
        summarise(
            n = n(),
            mean_ratio = mean(Ratio, na.rm = TRUE),
            sd_ratio = sd(Ratio, na.rm = TRUE),
            se_ratio = sd_ratio / sqrt(n),
            .groups = "drop"
        )
}


# エラーバー付き散布図を作成する関数
plot_ratio_scatter <- function(ratio_df,
                               title = "Relative Abundance by Genus",
                               x_label = "Genus",
                               y_label = "Relative Abundance (Ratio)",
                               point_size = 3,
                               point_alpha = 0.5,
                               error_bar_width = 0.2) {
    # 統計情報を計算
    stats_df <- calculate_stats(ratio_df)

    # プロット作成
    p <- ggplot() +
        geom_jitter(
            data = ratio_df,
            aes(x = Genus, y = Ratio, color = Genus),
            width = 0.1,
            height = 0,
            alpha = point_alpha,
            size = point_size
        ) +
        # 平均値
        geom_point(
            data = stats_df,
            aes(x = Genus, y = mean_ratio, color = Genus),
            size = point_size * 2,
            shape = 20 # ●
        ) +
        geom_errorbar(
            data = stats_df,
            aes(
                x = Genus,
                ymin = mean_ratio - se_ratio,
                ymax = mean_ratio + se_ratio,
                color = Genus
            ),
            width = error_bar_width,
            linewidth = 1
        ) +
        labs(
            title = title,
            x = x_label,
            y = y_label
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 13, face = "bold"),
            legend.position = "right",
            legend.title = element_text(face = "bold"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_line(color = "grey95")
        )

    return(p)
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

plot_ratio_scatter(ratio_df, title = "Relative Abundance of Bartonella") %>%
    ggsave("relative_bartonella.svg", plot = .)
