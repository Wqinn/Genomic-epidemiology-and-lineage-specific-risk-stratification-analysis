## snp matrix
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)

infile <- "fig3/snp.xlsx"
df <- read_xlsx(infile, guess_max = 50000)

strain <- trimws(as.character(df[[1]]))
src    <- trimws(as.character(df[[2]]))
names(src) <- strain

mat_df <- df[, -(1:2)]
mat <- as.matrix(sapply(mat_df, function(x) as.numeric(gsub(",", "", as.character(x)))))
rownames(mat) <- strain
colnames(mat) <- trimws(colnames(mat_df))

common <- intersect(rownames(mat), colnames(mat))
mat <- mat[common, common, drop = FALSE]
src <- src[common]

mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
diag(mat) <- 0

hc  <- hclust(as.dist(mat), method = "average")
ord <- hc$order
mat2 <- mat[ord, ord, drop = FALSE]
src2 <- src[ord]

mx <- max(mat2, na.rm = TRUE)
col_fun <- colorRamp2(c(0, mx), c("white", "#000080"))

source_cols <- c(
  "IPF_pig"                   = "#1f77b4",
  "BYF_pig"                   = "#aec7e8",
  "slaughterhouse"            = "#ff7f0e",
  "slaughterhouse wastewater" = "#2ca02c",
  "market wastewater"         = "#98df8a",
  "pork"                      = "#ffbb78",
  "transport vehicles"        = "#9467bd",
  "pig workers"               = "#d62728",
  "diarrhea patients"         = "#ff9896",
  "healthy human"             = "#c5b0d5"
)

missing_src <- setdiff(unique(src2), names(source_cols))
if (length(missing_src) > 0) {
  source_cols <- c(source_cols, setNames(rep("#DDDDDD", length(missing_src)), missing_src))
}

ha_row <- rowAnnotation(
  Source = src2,
  col = list(Source = source_cols),
  annotation_name_side = "top"
)

ha_col <- HeatmapAnnotation(
  Source = src2,
  col = list(Source = source_cols),
  show_annotation_name = FALSE
)

ht <- Heatmap(
  mat2,
  name = "SNP",
  col = col_fun,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = ha_row,
  top_annotation  = ha_col,
  use_raster = TRUE, raster_quality = 2,
  rect_gp = gpar(col = "grey90", lwd = 0.2),
  heatmap_legend_param = list(title = "SNP", at = pretty(c(0, mx), n = 4))
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")


## network
library(readxl)
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)

infile <- "fig3/snp.xlsx"
df <- read_xlsx(infile, guess_max = 50000)

strain <- trimws(as.character(df[[1]]))
src    <- trimws(as.character(df[[2]]))
names(src) <- strain

mat_df <- df[, -(1:2)]
mat <- as.matrix(sapply(mat_df, function(x) as.numeric(gsub(",", "", as.character(x)))))
rownames(mat) <- strain
colnames(mat) <- trimws(colnames(mat_df))

common <- intersect(rownames(mat), colnames(mat))
mat <- mat[common, common, drop = FALSE]
src <- src[common]

mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
diag(mat) <- 0

thr <- 10
idx <- which(mat <= thr & mat > 0 & upper.tri(mat), arr.ind = TRUE)

pairs <- data.frame(
  i = rownames(mat)[idx[,1]],
  j = colnames(mat)[idx[,2]],
  snp = mat[idx],
  si = src[rownames(mat)[idx[,1]]],
  sj = src[colnames(mat)[idx[,2]]],
  stringsAsFactors = FALSE
)

within_counts <- pairs %>%
  filter(si == sj) %>%
  count(Source = si, name = "within_clonal_pairs")

between_counts <- pairs %>%
  filter(si != sj) %>%
  transmute(
    from = pmin(si, sj),
    to   = pmax(si, sj)
  ) %>%
  count(from, to, name = "between_clonal_pairs")

source_cols <- c(
  "IPF_pig"                   = "#1f77b4",
  "BYF_pig"                   = "#aec7e8",
  "slaughterhouse"            = "#ff7f0e",
  "slaughterhouse wastewater" = "#2ca02c",
  "market wastewater"         = "#98df8a",
  "pork"                      = "#ffbb78",
  "transport vehicles"        = "#9467bd",
  "pig workers"               = "#d62728",
  "diarrhea patients"         = "#ff9896",
  "healthy human"             = "#c5b0d5"
)

sources_all <- sort(unique(src))
verts <- data.frame(
  name = sources_all,
  within_clonal_pairs = 0,
  stringsAsFactors = FALSE
)

m <- match(within_counts$Source, verts$name)
verts$within_clonal_pairs[m] <- within_counts$within_clonal_pairs

miss <- setdiff(verts$name, names(source_cols))
if (length(miss) > 0) source_cols <- c(source_cols, setNames(rep("#DDDDDD", length(miss)), miss))

if (nrow(between_counts) == 0) {
  stop("未检测到")
}

g <- graph_from_data_frame(between_counts, directed = FALSE, vertices = verts)

set.seed(1)

p <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = between_clonal_pairs),
                 color = "grey60", alpha = 0.8) +
  geom_node_point(aes(size = within_clonal_pairs, color = name),
                  stroke = 0.3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = source_cols, guide = "none") +
  scale_edge_width(range = c(0.5, 3.5)) +   
  scale_size(range = c(5, 14)) +          
  theme_void()

print(p)
