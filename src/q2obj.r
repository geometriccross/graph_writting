library(tidyverse)
library(phyloseq)
library(microViz)
library(qiime2R)

MAX_DEPTH <- 3

# .qza ファイルの QIIME2 型（type）を取得
# .qza は ZIP アーカイブで、<uuid>/metadata.yaml に型情報が記録されている
get_q2_type <- function(qza_path) {
    tryCatch(
        {
            tmp <- tempfile()
            dir.create(tmp)
            on.exit(unlink(tmp, recursive = TRUE))

            # read_qza と同じ手法
            unpacked <- unzip(qza_path, list = TRUE)
            uuid_dir <- gsub("/.*", "", unpacked$Name[1])
            yaml_path <- file.path(tmp, uuid_dir, "metadata.yaml")
            unzip(qza_path, files = file.path(uuid_dir, "metadata.yaml"), exdir = tmp)

            lines <- readLines(yaml_path)
            match <- str_match(paste(lines, collapse = "\n"), "type:\s*(.+?)\s*$")
            if (is.na(match[1, 2])) stop("type field not found in metadata.yaml")
            match[1, 2]
        },
        error = function(e) stop(sprintf("Failed to read q2 type from %s: %s", qza_path, e$message))
    )
}

# .qza ファイルの型が期待値と一致するか検証
assert_q2_type <- function(qza_path, expected_type) {
    actual <- get_q2_type(qza_path)
    if (actual != expected_type) {
        stop(sprintf(
            "QIIME2 type mismatch for %s\n  expected: %s\n  actual:   %s",
            basename(qza_path), expected_type, actual
        ))
    }
    invisible(actual)
}

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

    # QIIME2 型の検証
    assert_q2_type(features_path, "FeatureTable[Frequency]")
    assert_q2_type(taxonomy_path, "FeatureData[Taxonomy]")
    assert_q2_type(tree_path,     "Phylogeny[Rooted]")

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
