# Writes biotapy's R golden files and R-only test fixtures. Run only in tests/r/Dockerfile's image, from the
# repo root: docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/work biotapy-golden
suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
})

write_golden <- function(df, path) {
  # write.csv keeps 15 significant digits; zlib's gzip header has no timestamp, so reruns are identical.
  con <- gzfile(path, "w")
  write.csv(df, con, row.names = FALSE)
  close(con)
  if (file.size(path) >= 1e6) stop(path, " is ", file.size(path), " bytes; rules.md R6.6 caps data files at 1 MB")
}

samples_as_rows <- function(physeq) {
  m <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) t(m) else m
}

dense_frame <- function(physeq) {
  m <- samples_as_rows(physeq)
  data.frame(sample_id = rownames(m), m, check.names = FALSE, row.names = NULL)
}

for (dir in c("tests/golden/global_patterns", "tests/data/phyloseq", "tests/data/dada2")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

## GlobalPatterns golden files (derived numbers only)
data(GlobalPatterns)
rel <- samples_as_rows(transform_sample_counts(GlobalPatterns, function(x) x / sum(x)))
nz <- which(rel != 0, arr.ind = TRUE)
nz <- nz[order(nz[, "col"], nz[, "row"]), , drop = FALSE]
write_golden(
  data.frame(sample_id = rownames(rel)[nz[, "row"]], taxon_id = colnames(rel)[nz[, "col"]], value = rel[nz]),
  "tests/golden/global_patterns/relative.csv.gz"
)
write_golden(dense_frame(tax_glom(GlobalPatterns, "Phylum")), "tests/golden/global_patterns/tax_glom_phylum.csv.gz")
write_golden(dense_frame(tax_glom(GlobalPatterns, "Genus")), "tests/golden/global_patterns/tax_glom_genus.csv.gz")

## Synthetic phyloseq fixtures: biotapy's toy() numbers, no third-party data
counts <- rbind(
  c(10, 5, 20, 30, 0, 2, 1, 0), c(8, 7, 25, 22, 3, 0, 0, 1), c(12, 4, 18, 35, 1, 5, 2, 0),
  c(2, 1, 5, 10, 0, 40, 15, 3), c(0, 2, 3, 12, 2, 38, 20, 5), c(1, 0, 4, 8, 1, 45, 12, 2)
)
dimnames(counts) <- list(paste0("s", 1:6), paste0("f", 1:8))
firm <- c("Bacteria", "Firmicutes", "Clostridia")
bact <- c("Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales")
prot <- c("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae")
taxonomy <- rbind(
  c(firm, "Lachnospirales", "Lachnospiraceae", "Blautia"), c(firm, "Lachnospirales", "Lachnospiraceae", "Roseburia"),
  c(firm, "Oscillospirales", "Ruminococcaceae", "Faecalibacterium"), c(bact, "Bacteroidaceae", "Bacteroides"),
  c(bact, "Bacteroidaceae", "Bacteroides"), c(bact, "Prevotellaceae", "Prevotella"),
  c(prot, "Escherichia"), c(prot, NA)
)
dimnames(taxonomy) <- list(paste0("f", 1:8), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
newick <- "(((f1:0.1,f2:0.12)n4:0.05,f3:0.2)n1:0.1,((f4:0.02,f5:0.03)n5:0.05,f6:0.15)n2:0.1,(f7:0.1,f8:0.12)n3:0.2)root;"
groups <- sample_data(data.frame(group = rep(c("A", "B"), each = 3), row.names = paste0("s", 1:6)))

toy <- phyloseq(otu_table(t(counts), taxa_are_rows = TRUE), tax_table(taxonomy), groups, phy_tree(ape::read.tree(text = newick)))
toy_b <- prune_samples(c("s1", "s2"), toy)
saveRDS(toy, "tests/data/phyloseq/toy.rds")
save(toy, file = "tests/data/phyloseq/toy.RData")
save(toy, toy_b, file = "tests/data/phyloseq/two_objects.RData")
saveRDS(phyloseq(otu_table(counts, taxa_are_rows = FALSE), groups), "tests/data/phyloseq/samples_as_rows.rds")
sequences <- DNAStringSet(setNames(c("ACGT", "ACGA", "ACGC", "ACGG", "TCGT", "TCGA", "TCGC", "TCGG"), paste0("f", 1:8)))
saveRDS(merge_phyloseq(toy, sequences), "tests/data/phyloseq/with_refseq.rds")
saveRDS(prune_samples("s1", toy), "tests/data/phyloseq/single_sample.rds")
zeroed <- counts
zeroed["s6", ] <- 0
saveRDS(phyloseq(otu_table(t(zeroed), taxa_are_rows = TRUE), tax_table(taxonomy)), "tests/data/phyloseq/zero_sample.rds")

## DADA2 .rds fixtures: the same values as the CSV strings in tests/io/test_dada2.py
seqtab <- rbind(S1 = c(10L, 0L, 5L, 0L), S2 = c(0L, 0L, 0L, 0L), S3 = c(3L, 7L, 1L, 0L))
colnames(seqtab) <- c("ACGTACGT", "TTGACCAA", "GGGCCCAA", "CCCCAAAA")
saveRDS(seqtab, "tests/data/dada2/seqtab.rds")
taxa <- rbind(ACGTACGT = c("Bacteria", "Firmicutes", "Blautia"), TTGACCAA = c("Bacteria", "Bacteroidota", NA))
colnames(taxa) <- c("Kingdom", "Phylum", "Genus")
saveRDS(taxa, "tests/data/dada2/taxa.rds")

writeLines(c(
  R.version.string,
  paste("Bioconductor", as.character(BiocManager::version())),
  paste0("phyloseq ", packageVersion("phyloseq"))
), "tests/golden/VERSIONS.txt")
