library(tidyverse)
library(phyloseq)
library(microViz)
library(qiime2R)

MAX_DEPTH <- 3

# フォルダ内からパターンに一致するファイルを再帰的に探索
find_file <- function(folder, pattern) {
    files <- list.files(
        folder,
        pattern = pattern,
        recursive = TRUE,
        full.names = TRUE
    )
    # max depth でフィルタ（folder 自体を depth 0 とする）
    files <- files[map_int(files, \(f) length(strsplit(sub(paste0("^", folder, "/?"), "", f), "/")[[1]])) <= MAX_DEPTH]
    if (length(files) == 0) {
        stop(sprintf("No file matching '%s' found in %s (max depth %d)", pattern, folder, MAX_DEPTH))
    }
    if (length(files) > 1) {
        warning(sprintf("Multiple files matched '%s'; using the first: %s", pattern, files[1]))
    }
    files[1]
}

# QIIME2メタデータの読み込みとクリーニング（型定義行の削除とID列の修正）
clean_sample_data <- function(folder) {
    metadata_path <- find_file(folder, "metadata\\.tsv$")
    read_tsv(metadata_path, comment = "") %>%
        rename_with(~"SampleID", .cols = 1) %>%
        filter(!str_detect(SampleID, "^#q2:types|^categorical|^numeric")) %>%
        column_to_rownames("SampleID") %>%
        sample_data()
}


# Feature table等の読み込み、分類名の修正、メタデータの結合
load_q2obj <- function(folder) {
    features_path <- find_file(folder, "common_biology_free_table.qza$")
    taxonomy_path <- find_file(folder, "common_biology_free_classification.qza$")
    tree_path     <- find_file(folder, "biology_free_tree.qza$")
    qza_to_phyloseq(
        features = features_path,
        taxonomy = taxonomy_path,
        tree     = tree_path
    )
}

tax_tweak <- function(phyloseq_object, folder) {
    phyloseq_object %>%
        tax_fix() %>%
        tax_fix(unknowns = c("uncultured")) %>%
        merge_phyloseq(clean_sample_data(folder))
}
