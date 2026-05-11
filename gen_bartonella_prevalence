#!/usr/bin/env Rscript

show_help <- function() {
    cat("Usage: Rscript prevalence.r <qza_folder> <metadata_path> [output_dir]\n\n")
    cat("Arguments:\n")
    cat("  qza_folder      QIIME2 .qzaファイルが格納されたフォルダのパス\n")
    cat("  metadata_path   メタデータファイルのパス (TSV, sample_name列とspecies列が必要)\n")
    cat("  output_dir      出力先ディレクトリ (デフォルト: カレントディレクトリ)\n\n")
    cat("Options:\n")
    cat("  --help, -h      このヘルプメッセージを表示して終了する\n\n")
    cat("Outputs:\n")
    cat("  bartonella_prevalence_by_species.csv  種ごとのBartonella保有率\n")
    cat("  bartonella_prevalence_by_species.svg  種ごとのBartonella保有率 棒グラフ\n")
    quit(status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && (args[1] == "--help" || args[1] == "-h")) {
    show_help()
}

library(conflicted)
library(tidyverse)
library(phyloseq)
library(dplyr)

source("./src/q2obj.r")

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")


# Bartonella相対存在比をサンプルごとに計算
# scatter.r の get_ratio(df, "Bartonella") と同等
bartonella_ratio <- function(df) {
    bartonella_abundance <- df %>%
        filter(Genus == "Bartonella") %>%
        select(-FeatureID, -Genus) %>%
        colSums() %>%
        enframe(name = "SampleID", value = "Abundance_a")

    total_abundance <- df %>%
        select(-FeatureID, -Genus) %>%
        colSums() %>%
        enframe(name = "SampleID", value = "Abundance_b")

    bartonella_abundance %>%
        left_join(total_abundance, by = "SampleID") %>%
        mutate(Ratio = Abundance_a / Abundance_b) %>%
        select(SampleID, Ratio)
}


# speciesごとの保有率を計算
calculate_prevalence <- function(ratio_df, sample_info) {
    ratio_df %>%
        left_join(sample_info, by = "SampleID") %>%
        filter(!is.na(species)) %>%
        group_by(species) %>%
        summarise(
            n = n(),
            detected = sum(Ratio > 0),
            prevalence = detected / n,
            mean_ratio = mean(Ratio, na.rm = TRUE),
            sd_ratio = sd(Ratio, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(desc(prevalence))
}


# 種ごとの保有率棒グラフを作成
plot_prevalence <- function(prevalence_df, title = "Bartonella保有率（種別）") {
    # 表示用にアンダースコアをスペースに置換
    plot_df <- prevalence_df %>%
        mutate(
            species_label = gsub("_", " ", species),
            species_label = factor(species_label, levels = gsub("_", " ", species))
        )

    ggplot(plot_df, aes(x = species_label, y = prevalence, fill = species_label)) +
        geom_col(width = 0.7) +
        geom_text(
            aes(label = sprintf("n=%d", n)),
            vjust = -0.5, size = 3.5
        ) +
        geom_text(
            aes(label = sprintf("%.1f%%", prevalence * 100)),
            vjust = -1.8, size = 3.5, fontface = "bold"
        ) +
        scale_y_continuous(
            limits = c(0, 1.15),
            labels = scales::percent_format()
        ) +
        labs(
            title = title,
            x = "Species",
            y = "Bartonella prevalence"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            legend.position = "none"
        )
}


# --- Main ---

if (length(args) < 2) {
    stop("Usage: Rscript prevalence.r <qza_folder> <metadata_path> [output_dir]\nUse --help for more information.")
}
qza_folder <- args[1]
metadata_path <- args[2]
output_dir <- if (length(args) >= 3) args[3] else "."

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# phyloseqオブジェクト構築 (OTU + Taxonomy)
ps <- load_q2obj(qza_folder) %>%
    tax_fix() %>%
    tax_fix(unknowns = c("uncultured"))

# メタデータ読み込み
otu_samples <- sample_names(ps)
sample_info <- read_tsv(metadata_path, comment = "", show_col_types = FALSE) %>%
    rename_with(~"SampleID", .cols = 1) %>%
    filter(!str_detect(.data[["SampleID"]], "^#q2:types|^categorical|^numeric")) %>%
    filter(SampleID %in% otu_samples) %>%
    select(SampleID, species)

# OTU + Tax 結合
otu <- ps %>%
    otu_table() %>%
    as.data.frame() %>%
    rownames_to_column("FeatureID")

tax <- ps %>%
    tax_table() %>%
    as.data.frame() %>%
    select("Genus") %>%
    rownames_to_column("FeatureID")

df <- otu %>% left_join(tax, by = "FeatureID")

# 種ごとのBartonella保有率を計算
prevalence_result <- bartonella_ratio(df) %>%
    calculate_prevalence(sample_info)

# CSV出力
write_csv(prevalence_result, file.path(output_dir, "bartonella_prevalence_by_species.csv"))

# 棒グラフ出力
prevalence_result %>%
    plot_prevalence() %>%
    ggsave(
        file.path(output_dir, "bartonella_prevalence_by_species.svg"),
        plot = .,
        width = 10, height = 6, dpi = 600
    )

cat("Done.\n")
