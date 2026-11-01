suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: export_imvigor210.R cds.RData signatures.tsv output_dir")
load(args[[1]])
signatures <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
genes <- unique(c("ITGB1", trimws(unlist(strsplit(signatures$genes, ",", fixed = TRUE)))))
symbols <- as.character(rowData(cds)$symbol)
if (all(is.na(symbols)) || !any(symbols == "ITGB1", na.rm = TRUE)) symbols <- as.character(rowData(cds)$Symbol)
normalized <- counts(cds, normalized = TRUE)
ranked <- apply(normalized, 2, rank, ties.method = "average") / nrow(normalized)
rownames(ranked) <- symbols
feature_table <- data.frame(sample_id = colnames(normalized), check.names = FALSE)
feature_table$ITGB1 <- colMeans(log2(normalized[symbols == "ITGB1", , drop = FALSE] + 1))
for (i in seq_len(nrow(signatures))) {
  signature_genes <- setdiff(trimws(unlist(strsplit(signatures$genes[[i]], ",", fixed = TRUE))), "ITGB1")
  selected <- which(symbols %in% signature_genes)
  feature_table[[signatures$signature[[i]]]] <- if (length(selected) >= 2) colMeans(ranked[selected, , drop = FALSE]) else NA_real_
}
keep <- !is.na(symbols) & symbols %in% genes
matrix <- normalized[keep, , drop = FALSE]
rownames(matrix) <- symbols[keep]
matrix <- rowsum(matrix, group = rownames(matrix), reorder = FALSE) /
  as.numeric(table(factor(rownames(matrix), levels = unique(rownames(matrix)))))
matrix <- log2(matrix + 1)
dir.create(args[[3]], recursive = TRUE, showWarnings = FALSE)
write.table(matrix, gzfile(file.path(args[[3]], "expression_log2.tsv.gz")), sep = "\t", quote = FALSE, col.names = NA)
write.table(feature_table, file.path(args[[3]], "analysis_features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
metadata <- as.data.frame(colData(cds)); metadata$sample_id <- rownames(metadata)
write.table(metadata, file.path(args[[3]], "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
