library(tidyverse)
library(phyloseq)
library(microViz)
library(qiime2R)


pseq <- qza_to_phyloseq(
    features = "./q2objs/common_biology_free_table.qza",
    taxonomy = "./q2objs/common_biology_free_classification.qza",
    metadata = "./q2objs/metadata.tsv",
    tree = "./q2objs/biology_free_rooted_tree.qza"
) %>%
    tax_fix() %>%
    tax_fix(unknowns = c("uncultured"))

metadata <- read_tsv("./q2objs/metadata.tsv", comment = "") %>%
    rename(SampleID = `#SampleID`) %>%
    column_to_rownames("SampleID")

sample_data(pseq) <- sample_data(metadata)

print(colnames(sample_data(pseq)))

p <- pseq %>%
    ps_arrange("RawID") %>%
    comp_barplot(
        tax_level = "Genus",
        n_taxa = 32,
        bar_outline_colour = NA,
        sample_order = "asis",
        label = "RawID"
    ) +
    coord_flip()


ggsave("barplot.svg", plot = p, width = 10, height = 8, dpi = 300)
