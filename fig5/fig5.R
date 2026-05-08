## upset_args
req_pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "ggplot2", "UpSetR")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
}
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(UpSetR)

path <- "fig5\\plasmid.xlsx"
raw  <- readxl::read_excel(path)

names(raw) <- trimws(names(raw))
normalize_name <- function(x) gsub("[^a-z]", "", tolower(trimws(x)))
nm <- normalize_name(names(raw))
seq_keys <- normalize_name(c("Plasmid Sequence"))

arg_keys <- normalize_name(c("Resistance gene"))

seq_idx <- which(nm %in% seq_keys)
arg_idx <- which(nm %in% arg_keys)

seq_col <- names(raw)[seq_idx[1]]
arg_col <- names(raw)[arg_idx[1]]

df <- raw %>%
  transmute(
    SEQUENCE = as.character(.data[[seq_col]]),
    ARG_raw  = as.character(.data[[arg_col]])
  ) %>%
  filter(!is.na(SEQUENCE), SEQUENCE != "", !is.na(ARG_raw), ARG_raw != "")

split_regex <- "\\s*(,|;|/|\\||\\+)\\s*"
arg_long <- df %>%
  separate_rows(ARG_raw, sep = split_regex) %>%
  mutate(
    gene = str_trim(ARG_raw),
    gene = str_replace_all(gene, "\\s+", ""),           
    gene = str_replace_all(gene, "^\\((.*)\\)$", "\\1") 
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(SEQUENCE, gene)                              

top_genes <- 15  
genes_keep <- arg_long %>%
  count(gene, sort = TRUE, name = "freq") %>%
  slice_head(n = top_genes) %>%
  pull(gene)

arg_long_top <- arg_long %>% filter(gene %in% genes_keep)

## UpSetR
mat <- arg_long_top %>%
  mutate(value = 1L) %>%
  tidyr::pivot_wider(
    id_cols = SEQUENCE,
    names_from = gene,
    values_from = value,
    values_fill = 0L
  )

mat_df <- as.data.frame(mat)
rownames(mat_df) <- mat_df$SEQUENCE
mat_df$SEQUENCE <- NULL


for (cn in names(mat_df)) {
  mat_df[[cn]] <- as.integer(mat_df[[cn]] > 0)
}

col_lightblue   <- "#a6cee3"   
col_lightpurple <- "#b39ddb"   

while (!is.null(dev.list())) dev.off()
try(dev.new(width = 10, height = 8), silent = TRUE)

old_theme <- ggplot2::theme_get()
ggplot2::theme_set(ggplot2::theme_gray())
on.exit(ggplot2::theme_set(old_theme), add = TRUE)

top_intersections <- 30  

suppressWarnings(
  UpSetR::upset(
    mat_df,
    nsets = ncol(mat_df),
    nintersects = top_intersections,
    order.by = "freq",
    mb.ratio = c(0.65, 0.35),
    text.scale = c(1.2, 1.2, 1.0, 1.0, 1.1, 1.0),
    main.bar.color = col_lightblue,
    sets.bar.color = col_lightblue,
    matrix.color   = col_lightpurple,
    point.size = 3.2,
    line.size  = 0.7
  )
)

## 弦图
req_pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "circlize", "randomcoloR")
for (p in req_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(readxl); library(dplyr); library(tidyr); library(stringr)
library(circlize); library(randomcoloR)

path <- "fig5\\plasmid.xlsx"
raw <- readxl::read_excel(path)

normalize_name <- function(x) gsub("[^a-z]", "", tolower(x))
nm <- normalize_name(names(raw))
rep_idx <- which(nm %in% c("plasmidreplicon", "replicon", "plasmidtype"))
org_idx <- which(nm == "organism")
stopifnot(length(rep_idx) >= 1, length(org_idx) >= 1)
rep_col <- names(raw)[rep_idx[1]]
org_col <- names(raw)[org_idx[1]]

df <- raw %>%
  select(replicon = all_of(rep_col), organism = all_of(org_col)) %>%
  mutate(replicon = as.character(replicon),
         organism = as.character(organism)) %>%
  filter(!is.na(replicon), !is.na(organism),
         str_trim(replicon) != "", str_trim(organism) != "")

split_regex <- "[,;/|]+"
df_long <- df %>%
  mutate(replicon = str_replace_all(replicon, "\\s+", " ")) %>%
  separate_rows(replicon, sep = split_regex) %>%
  mutate(
    replicon = str_trim(replicon),
    replicon = str_replace_all(replicon, "^inc", "Inc"),
    organism = str_trim(organism)
  ) %>%
  filter(replicon != "", organism != "")

filter_by      <- "both"      
top_replicon   <- 8           
top_organism   <- 8           

df_long_filt <- df_long

if (filter_by %in% c("replicon", "both")) {
  keep_rep <- df_long %>%
    count(replicon, sort = TRUE, name = "freq") %>%
    slice_head(n = top_replicon) %>%
    pull(replicon)
  df_long_filt <- df_long_filt %>% filter(replicon %in% keep_rep)
}

if (filter_by %in% c("organism", "both")) {
  keep_org <- df_long %>%
    count(organism, sort = TRUE, name = "freq") %>%
    slice_head(n = top_organism) %>%
    pull(organism)
  df_long_filt <- df_long_filt %>% filter(organism %in% keep_org)
}

stopifnot(nrow(df_long_filt) > 0)

pair_counts <- df_long_filt %>%
  count(replicon, organism, name = "n") %>%
  arrange(desc(n))

min_count <- 1     
max_edges <- 40    
pair_filtered <- pair_counts %>% filter(n >= min_count) %>% slice_head(n = max_edges)
if (nrow(pair_filtered) == 0) pair_filtered <- pair_counts

pair_plot <- pair_filtered %>%
  transmute(from = paste0("Rep: ", replicon),
            to   = paste0("Org: ", organism),
            value = n)

rep_nodes <- unique(pair_plot$from)
org_nodes <- unique(pair_plot$to)
plot_order <- c(rep_nodes, org_nodes)

set.seed(2025)
grid_col <- c(
  setNames(randomcoloR::distinctColorPalette(length(rep_nodes)), rep_nodes),
  setNames(randomcoloR::distinctColorPalette(length(org_nodes)), org_nodes)
)

while (!is.null(dev.list())) dev.off()
try(dev.new(width = 10, height = 10), silent = TRUE)
if (.Platform$OS.type == "windows") try(windows(10,10), silent = TRUE)

circos.clear()

gap <- rep(2, length(plot_order))
gap[length(rep_nodes)]  <- 8
gap[length(plot_order)] <- 8
circos.par(gap.after = gap, start.degree = 90, track.margin = c(0.01, 0.01))

chordDiagram(
  x = pair_plot,
  order = plot_order,
  grid.col = grid_col,
  transparency = 0.25,
  directional = 0,
  annotationTrack = "grid",                 
  preAllocateTracks = list(track.height = 0.09)
)

circos.trackPlotRegion(
  track.index = get.current.track.index(),
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    label  <- gsub("^Rep: |^Org: ", "", sector)  
    
    circos.text(
      x = get.cell.meta.data("xcenter"),
      y = get.cell.meta.data("ylim")[2] + mm_y(2.5),  
      labels = label,
      facing = "clockwise", niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.9, font = 2
    )
  },
  bg.border = NA
)

title(sprintf("Chord: Replicon (%s) \u2194 Organism (%s)",
              if (exists("filter_by") && filter_by %in% c("replicon","both")) paste0("Top-", top_replicon) else "All",
              if (exists("filter_by") && filter_by %in% c("organism","both")) paste0("Top-", top_organism) else "All"),
      cex.main = 1.1)


## TA heatmap
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

library(tidyverse)
library(pheatmap)
library(stringr)
library(readxl)
library(writexl)

setwd("fig5")

blast_res <- read.table("plasmid_TA_blastp.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
colnames(blast_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")

annot_df <- read.csv("TADB_annotation.csv", header = TRUE, stringsAsFactors = FALSE)
if(ncol(annot_df) >= 3) colnames(annot_df)[1:3] <- c("sseqid", "Accession", "GeneName")

if(file.exists("plasmid.xlsx")) {
  metadata <- read_excel("plasmid.xlsx", sheet = 1)
} else if(file.exists("plasmid.xlsx - Sheet1.csv")) {
  metadata <- read.csv("plasmid.xlsx - Sheet1.csv", stringsAsFactors = FALSE)

id_col <- which(str_detect(colnames(metadata), regex("Sequence|ID|Name", ignore_case = TRUE)))[1]
rep_col <- which(str_detect(colnames(metadata), regex("Replicon|Inc", ignore_case = TRUE)))[1]

if(is.na(id_col)) id_col <- 1
if(is.na(rep_col)) rep_col <- 3

colnames(metadata)[id_col] <- "PlasmidID"
colnames(metadata)[rep_col] <- "Replicon"

blast_filtered <- blast_res %>% filter(pident > 80, qcovs > 80)

blast_annotated <- left_join(blast_filtered, annot_df, by = "sseqid")

ta_keywords <- list(
  "VapC" = "VapBC", "VapB" = "VapBC", "VapD" = "VapBC", "VapX" = "VapBC",
  "RelE" = "RelBE", "RelB" = "RelBE", "ParE" = "RelBE", "RelJ" = "RelBE", "RelF" = "RelBE",
  "MazF" = "MazEF", "MazE" = "MazEF", "PemK" = "MazEF", "PemI" = "MazEF", "ChpB" = "MazEF", "ChpS" = "MazEF",
  "HigB" = "HigBA", "HigA" = "HigBA",
  "HipA" = "HipAB", "HipB" = "HipAB",
  "CcdB" = "CcdAB", "CcdA" = "CcdAB",
  "ParD" = "ParDE", 
  "Phd"  = "Phd/Doc", "Doc"  = "Phd/Doc", "YefM" = "Phd/Doc", "death-on-curing" = "Phd/Doc",
  "MqsR" = "MqsRA", "MqsA" = "MqsRA",
  "HicA" = "HicAB", "HicB" = "HicAB",
  "PrlF" = "PrlF/YhaV", "YhaV" = "PrlF/YhaV", "Schmidt" = "PrlF/YhaV", 
  "Hha"  = "Hha/TomB",  "TomB" = "Hha/TomB",  "YgfX" = "Hha/TomB",
  "GhoT" = "GhoST", "GhoS" = "GhoST", # Type V
  "CbtA" = "CbtA/CbeA", "CbeA" = "CbtA/CbeA",
  "YacA" = "YacAB", "YacB" = "YacAB",
  "TacA" = "TacA",
  "Hok" = "Hok/Sok", "Sok" = "Hok/Sok", "Gef" = "Hok/Sok", "Mok" = "Hok/Sok",
  "Ldr" = "Ldr", "LdrA" = "Ldr", "LdrB" = "Ldr", "LdrC" = "Ldr", "LdrD" = "Ldr",
  "TisB" = "TisB/IstR", "IstR" = "TisB/IstR",
  "Ibs"  = "Ibs/Sib",   "Sib"  = "Ibs/Sib",
  "ShoB" = "ShoB", 
  "SymE" = "SymE",
  "Fst"  = "Fst/Ltx",
  "YeeU" = "YeeUV", "YeeV" = "YeeUV",
  "PasT" = "PasTI", "PasI" = "PasTI", "RatA" = "PasTI",
  "BrnT" = "BrnTA", "BrnA" = "BrnTA",
  "PezT" = "PezAT", "PezA" = "PezAT",
  "DarT" = "DarTG", "DarG" = "DarTG",
  "AbiE" = "Abi system", "AbiG" = "Abi system", 
  "DinJ" = "RelBE", "YafQ" = "RelBE", "YafO" = "RelBE", "YafN" = "RelBE", 
  "GNAT" = "GNAT-family", 
  "SdhE" = "SdhE",        
  "CopG" = "CopG/RHH",    
  "ribbon-helix-helix" = "CopG/RHH"
)

get_ta_family <- function(desc_string) {
  if (is.na(desc_string) || desc_string == "") return("Other_TA")
  for (keyword in names(ta_keywords)) {
    if (str_detect(desc_string, regex(keyword, ignore_case = TRUE))) {
      return(ta_keywords[[keyword]])
    }
  }
  return("Other_TA")
}

blast_annotated$TA_Family <- sapply(blast_annotated$GeneName, get_ta_family)

real_ids <- metadata$PlasmidID

blast_annotated$PlasmidID <- sapply(blast_annotated$qseqid, function(q_id) {
  hit <- real_ids[str_detect(q_id, fixed(real_ids))]
  if(length(hit) > 0) return(hit[1]) else return(NA)
})

blast_final <- blast_annotated %>% filter(!is.na(PlasmidID))

ta_counts <- blast_final %>%
  count(PlasmidID, TA_Family) %>%
  pivot_wider(names_from = TA_Family, values_from = n, values_fill = 0)

full_data <- metadata %>%
  select(PlasmidID, Replicon) %>% 
  left_join(ta_counts, by = "PlasmidID")

full_data[is.na(full_data)] <- 0
write_xlsx(full_data, "Plasmid_TA_Result.xlsx")

plot_data <- full_data %>%
  group_by(Replicon) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  column_to_rownames("Replicon")

if("Other_TA" %in% colnames(plot_data)) {
  plot_data <- plot_data %>% select(-Other_TA)
}

plot_data <- plot_data[, colSums(plot_data) > 0, drop = FALSE]
plot_data_save <- plot_data %>%
  rownames_to_column(var = "Replicon")

if(nrow(plot_data) > 0 && ncol(plot_data) > 0) {
  pheatmap(plot_data, 
           cluster_cols = TRUE, 
           cluster_rows = TRUE, 
           display_numbers = FALSE, 
           number_format = "%.1f",
           color = colorRampPalette(c("white", "#E64B35"))(50), 
           border_color = "grey95",
           fontsize_row = 8,  
           fontsize_col = 10, 
           angle_col = 45)