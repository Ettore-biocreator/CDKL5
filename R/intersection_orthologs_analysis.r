library(dplyr)
library(tools)
library(ggplot2)
library(VennDiagram)
library(grid)
library(futile.logger)
library(openxlsx)

# ── Paths ─────────────────────────────────────────────────────────────────────
# Transcriptomics
mouse_deg_folder  <- "PATH/TO/transcriptomics/mouse_DEG"
zebra_deg_folder  <- "PATH/TO/transcriptomics/zebrafish_DEG"
cell_deg_folder   <- "PATH/TO/transcriptomics/cell_line_DEG"
human_deg_file    <- "PATH/TO/transcriptomics/human/DEG_Human_Patient.txt"

# Proteomics
mouse_dep_folder  <- "PATH/TO/proteomics/mouse_DEP"
zebra_dep_file    <- "PATH/TO/proteomics/zebrafish/DEP_Zebrafish.txt"
hela_dep_file     <- "PATH/TO/proteomics/cell_line/proteomics_results_HeLa.xlsx"
hek_dep_file      <- "PATH/TO/proteomics/cell_line/proteomics_results_HEK293.xlsx"

# Orthology
orthology_file    <- "PATH/TO/orthology/ORTHOLOGY-ALLIANCE_COMBINED.tsv"

# BioMart synonym tables (used for gene symbol harmonisation)
biomart_mouse_file <- "PATH/TO/biomart/biomart_mouse.txt"  # cols: Gene.name, Gene.Synonym
biomart_zebra_file <- "PATH/TO/biomart/biomart_zebrafish.txt"
# ──────────────────────────────────────────────────────────────────────────────

# ── Parameters ────────────────────────────────────────────────────────────────
padj_mouse   <- 0.1
padj_zebra   <- 0.05
padj_human   <- 0.05
padj_cell    <- 0.1
lfc_cell     <- 1.0    # |log2FC| threshold for cell line DEGs

# Column names for gene symbols (vary by dataset)
gene_col_mouse  <- "Gene.name"
gene_col_zebra  <- "gene_name"
gene_col_cell   <- "external_gene_name"

# HeLa / HEK sheet index and significance column
hela_sheet      <- 4
hela_sig_col    <- "Significantly.deregulated.HeLa_KO_HeLa_Rescue"
hek_sheet       <- 4
hek_sig_col     <- "Significantly.deregulated.KO_IBAQ_Res_IBAQ"

# Venn diagram
venn_categories <- c("Cell Line (Human)", "Mouse", "Zebrafish")
venn_colors     <- c("#73EC8B", "#FD8B51", "#FFD700")
# ──────────────────────────────────────────────────────────────────────────────


# ── Helper: load all .txt files from a folder into a named list ───────────────
load_txt_folder <- function(folder, sep = "\t") {
  files <- list.files(folder, pattern = "\\.txt$", full.names = TRUE)
  data  <- lapply(files, read.table, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  names(data) <- tools::file_path_sans_ext(basename(files))
  data
}

# ── Helper: synonym harmonisation ────────────────────────────────────────────
harmonise_genes <- function(genes, biomart_df) {
  syn_map <- setNames(biomart_df$Gene.name, biomart_df$Gene.Synonym)
  sapply(genes, function(g) {
    if (g %in% biomart_df$Gene.name) g
    else if (g %in% names(syn_map))  syn_map[[g]]
    else g
  })
}

# ── Helper: filter DEG list ───────────────────────────────────────────────────
filter_deg <- function(data_list, gene_col, padj_thr, lfc_thr = NULL) {
  lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]] %>%
      filter(complete.cases(padj), complete.cases(.data[[gene_col]]), padj < padj_thr)
    if (!is.null(lfc_thr))
      df <- df %>% filter(abs(log2FoldChange) > lfc_thr)
    df
  }) |> setNames(names(data_list))
}

# ── Helper: orthology mapping loop ───────────────────────────────────────────
map_orthologs <- function(data_list, alliance_df, gene_col) {
  lapply(seq_along(data_list), function(i) {
    data_list[[i]] %>%
      left_join(alliance_df, by = setNames("Gene1Symbol", gene_col)) %>%
      pull(Gene2Symbol) %>%
      na.omit() %>%
      as.character()
  }) |> setNames(names(data_list))
}


# ── Load transcriptomics ──────────────────────────────────────────────────────
data_list_Mouse <- load_txt_folder(mouse_deg_folder, sep = " ")   # Mouse uses space separator
data_list_Zebra <- load_txt_folder(zebra_deg_folder)
data_list_Cell  <- load_txt_folder(cell_deg_folder)

human_g <- read.delim2(human_deg_file)
human_g_filt <- human_g %>%
  mutate(padj = as.numeric(padj), log2FoldChange = as.numeric(log2FoldChange)) %>%
  filter(padj < padj_human)

# ── Load proteomics ───────────────────────────────────────────────────────────
data_list_Mouse_proteo <- load_txt_folder(mouse_dep_folder)

proteo_zebra     <- read.delim2(zebra_dep_file)
proteo_zebra_fil <- proteo_zebra %>% filter(complete.cases(Significant))

protein_hela   <- read.xlsx(hela_dep_file, sheet = hela_sheet)
protein_hela_f <- protein_hela[protein_hela[[hela_sig_col]] == "+", ]

protein_hek    <- read.xlsx(hek_dep_file, sheet = hek_sheet)
protein_hek_f  <- protein_hek[protein_hek[[hek_sig_col]] == "+", ]

# ── Load BioMart tables ───────────────────────────────────────────────────────
biomart_mouse <- read.table(biomart_mouse_file, header = TRUE, sep = "\t")
biomart_zebra <- read.table(biomart_zebra_file, header = TRUE, sep = "\t")

# ── Filter DEGs ───────────────────────────────────────────────────────────────
list_mouse_filt  <- filter_deg(data_list_Mouse, gene_col_mouse, padj_mouse)
list_zebra_filt  <- filter_deg(data_list_Zebra, gene_col_zebra, padj_zebra)
list_cell_filt   <- filter_deg(data_list_Cell,  gene_col_cell,  padj_cell, lfc_thr = lfc_cell)

list_mouse_proteo_filt <- lapply(data_list_Mouse_proteo, function(df) {
  df %>% filter(complete.cases(Significant))
})

# ── Load orthology (Alliance) ─────────────────────────────────────────────────
Alliance <- read.delim2(orthology_file, sep = "\t", skip = 15)

Alliance_step  <- Alliance %>% filter(Gene2SpeciesName == "Homo sapiens")
Alliance_zebra <- Alliance_step %>% filter(Gene1SpeciesName == "Danio rerio")  %>% select(Gene1Symbol, Gene2Symbol)
Alliance_Mouse <- Alliance_step %>% filter(Gene1SpeciesName == "Mus musculus") %>% select(Gene1Symbol, Gene2Symbol)

# ── Map to human orthologs ────────────────────────────────────────────────────
process_mouse <- map_orthologs(list_mouse_filt,       Alliance_Mouse, gene_col_mouse)
process_zebra <- map_orthologs(list_zebra_filt,       Alliance_zebra, gene_col_zebra)

# Mouse proteomics: harmonise gene symbols before ortholog mapping
list_mouse_proteo_harm <- lapply(list_mouse_proteo_filt, function(df) {
  df$Genes <- harmonise_genes(df$Genes, biomart_mouse)
  df
})
process_mouse_pro <- map_orthologs(list_mouse_proteo_harm, Alliance_Mouse, "Genes")

# Zebrafish proteomics
proteo_zebra_fil$Genes <- harmonise_genes(proteo_zebra_fil$Genes, biomart_zebra)
all_zebra_protein <- proteo_zebra_fil %>%
  left_join(Alliance_zebra, by = c("Genes" = "Gene1Symbol")) %>%
  pull(Gene2Symbol) %>%
  na.omit() %>%
  as.character()

# ── Build final gene sets ─────────────────────────────────────────────────────
all_mouse_gene    <- unique(unlist(process_mouse))
all_mouse_protein <- unique(unlist(process_mouse_pro))
all_zebra_gene    <- unique(unlist(process_zebra))

cell_rna <- unique(c(
  list_cell_filt$HEK_DEG[[gene_col_cell]],
  list_cell_filt$HELA21_DEG[[gene_col_cell]]
))
cell_pro <- unique(c(protein_hek_f$`T:.PG.Genes`, protein_hela_f$`T:.PG.Genes`))

A <- unique(c(cell_rna, cell_pro))
B <- unique(c(all_mouse_gene, all_mouse_protein))
C <- unique(c(all_zebra_gene, all_zebra_protein))

# ── Venn Diagram ──────────────────────────────────────────────────────────────
venn.plot <- draw.triple.venn(
  area1    = length(A),
  area2    = length(B),
  area3    = length(C),
  n12      = length(intersect(A, B)),
  n23      = length(intersect(B, C)),
  n13      = length(intersect(A, C)),
  n123     = length(Reduce(intersect, list(A, B, C))),
  category = venn_categories,
  fill     = venn_colors,
  cat.cex  = 1.2,
  cex      = 1.2
)
