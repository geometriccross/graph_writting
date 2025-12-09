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


# 統計情報を計算する関数（平均、標準偏差、標準誤差、サンプル数）
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


# エラーバー付き散布図を作成する関数
plot_ratio_scatter <- function(ratio_df,
                               title = "Relative Abundance by Genus",
                               x_label = "Genus",
                               y_label = "Relative Abundance (Ratio)",
                               point_size = 3,
                               point_alpha = 0.5,
                               error_bar_width = 0.2,
                               color_map = NULL) {
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
            aes(x = Genus, y = mean_ratio, color = "black"),
            size = point_size,
            shape = 3 # 十
        ) +
        geom_errorbar(
            data = stats_df,
            aes(
                x = Genus,
                ymin = mean_ratio - se_ratio,
                ymax = mean_ratio + se_ratio,
                color = "black"
            ),
            width = error_bar_width,
            linewidth = 0.5
        ) +
        labs(
            title = title,
            x = x_label,
            y = y_label
        ) +
        scale_color_manual(values = color_map) +
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

# ------------------------------------------------------------------------------------------------

# "Candidatus_Puchtella" = "#A6CEE3",
# "Bartonella" = "#1F78B4",
# "Prevotella" = "#B2DF8A",
# "Lactobacillus" = "#33A02C",
# "Blautia" = "#FB9A99",
# "Escherichia-Shigella" = "#E31A1C",
# "Faecalibacterium" = "#FDBF6F",
# "Streptococcus" = "#FF7F00",
# "Neisseriaceae Family" = "#CAB2D6",
# "Coxiella" = "#6A3D9A",
# "Butyricicoccaceae Family" = "#FFFF99",
# "Lachnospiraceae Family" = "#B15928",
# "Flavobacterium" = "#1ff8ff",
# "Pseudomonas" = "#1B9E77",
# "Succinivibrio" = "#D95F02",
# "Sphingomonas" = "#7570B3",
# "Ruminococcus" = "#E7298A",
# "Arsenophonus" = "#66A61E",
# "Corynebacterium" = "#E6AB02",
# "Bifidobacterium" = "#A6761D",
# "UCG-005" = "#666666",
# "Staphylococcus" = "#4b6a53",
# "Methylobacterium-Methylorubrum" = "#b249d5",
# "Alloprevotella" = "#7edc45",
# "Dorea" = "#5c47b8",
# "Lachnospiraceae_ND3007_group" = "#cfd251",
# "1174-901-12" = "#ff69b4",
# "Clostridia_UCG-014" = "#69c86c",
# "UCG-002" = "#cd3e50",
# "[Eubacterium]_coprostanoligenes_group" = "#83d5af",
# "[Eubacterium]_hallii_group" = "#da6130",
# "CAG-352" = "#5e79b2",
# "Anaerococcus" = "#c29545",
# "Nocardioides" = "#532a5a",
# "Erysipelotrichaceae_UCG-003" = "#5f7b35",
# "Subdoligranulum" = "#c497cf",
# "Peptoniphilus" = "#773a27",
# "Comamonadaceae Family" = "#7cb9cb",
# "Agathobacter" = "#594e50",
# "Anaerovibrio" = "#d3c4a8",
# "[Eubacterium]_ruminantium_group" = "#c17e7f",
# "Other" = "#D3D3D3"


target_map <- list(
    "Bartonella" = "#1F78B4",
    "Escherichia-Shigella" = "#E31A1C",
    "Pseudomonas" = "#1B9E77",
    "Streptococcus" = "#FF7F00",
    "Staphylococcus" = "#4b6a53",
    "Corynebacterium" = "#E6AB02",
    "Candidatus_Puchtella" = "#A6CEE3",
    "Coxiella" = "#6A3D9A"
)

names(target_map) %>%
    map_dfr(~ get_ratio(df, key = .x)) %>%
    mutate(Genus = factor(Genus, levels = names(target_map))) %>%
    plot_ratio_scatter(title = "占有率の比較", color_map = target_map) %>%
    ggsave("relative_abundance_comparison.svg", plot = .)

names(target_map) %>%
    map_dfr(~ get_ratio(df, key = .x) %>%
        calculate_stats()) %>%
    write_csv("relative_abundance_stats.csv")
