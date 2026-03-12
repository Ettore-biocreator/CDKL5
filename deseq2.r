rm(list = ls())

library(DESeq2)
library(dplyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
coldata_file  <- "PATH/TO/samplesheet.txt"
counts_file   <- "PATH/TO/rawcounts.txt"
annot_file    <- "PATH/TO/annotation.txt"
output_file   <- "PATH/TO/output/DEG_results.txt"
# ──────────────────────────────────────────────────────────────────────────────

# ── Parameters ────────────────────────────────────────────────────────────────
sample_filter      <- "3m_H"        # Pattern to select sample subset
col_range          <- 1:54          # Columns to retain from raw count matrix
min_counts         <- 10            # Minimum counts threshold for filtering
smallest_group     <- 4             # Minimum number of samples above threshold
contrast_group_exp <- "EXP6"        # Experimental group label (in 'Condition')
contrast_group_ctr <- "CTR6"        # Control group label (in 'Condition')
# ──────────────────────────────────────────────────────────────────────────────

# ── Load data ─────────────────────────────────────────────────────────────────
coldata <- read.table(coldata_file, fill = TRUE, sep = "\t", header = TRUE)
data    <- read.table(counts_file,  header = TRUE, check.names = FALSE,
                      fill = TRUE, row.names = 1)
data    <- data[, col_range]

stopifnot(all(coldata$Sample.Name == colnames(data)))

# ── Subset samples ────────────────────────────────────────────────────────────
data1    <- data[, grep(sample_filter, colnames(data))]
coldata1 <- coldata[grep(sample_filter, colnames(data)), ]

# ── DESeq2 ───────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = data1,
  colData   = coldata1,
  design    = ~ Condition
)

# Low-count filtering
keep <- rowSums(counts(dds) >= min_counts) >= smallest_group
dds  <- dds[keep, ]

dds <- DESeq(dds)

# ── Results ───────────────────────────────────────────────────────────────────
res        <- results(dds, contrast = c("Condition", contrast_group_exp, contrast_group_ctr))
resOrdered <- res[order(res$padj), ]

# ── Annotation ────────────────────────────────────────────────────────────────
resOrdered$Geneid <- rownames(resOrdered)

annot <- read.table(annot_file, header = TRUE, sep = "\t")

resOrdered <- as.data.frame(resOrdered) %>%
  left_join(annot, by = "Geneid") %>%
  select(c(7, 8, 2, 5, 6)) %>%
  mutate(Gene.name = ifelse(Gene.name == "" | is.na(Gene.name), Geneid, Gene.name))

# ── Export ────────────────────────────────────────────────────────────────────
write.table(resOrdered, file = output_file, quote = FALSE, row.names = FALSE, sep = "\t")
