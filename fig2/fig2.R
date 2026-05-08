## Antimicrobial susceptibility profiles
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)

file_path <- "fig2/AST_result.xlsx"
ast <- read_xlsx(file_path)
unique(ast$Source)

ast <- ast %>%
  mutate(
    source_group = case_when(
      Source %in% c("IPF_pig", "BYF_pig") ~ "Farm",
      Source %in% c("slaughterhouse", "slaughterhouse wastewater", "transport vehicles") ~ "Slaughter",
      Source %in% c("pork", "market wastewater") ~ "Retail",
      Source %in% c("pig workers", "diarrhea patients", "healthy human") ~ "Human",
      TRUE ~ NA_character_
    )
  )

table(ast$source_group)
abx_cols <- c("DOX", "TGC", "AMS", "CAZ", "CTX", "FOX",
              "CFZ", "CN", "CIP", "IPM", "SXT", "FFC",
              "C", "AZM", "AK")

res_by_source <- ast %>%
  filter(source_group %in% c("Farm", "Slaughter", "Retail", "Human")) %>%
  group_by(source_group) %>%
  summarise(
    n_isolates = n(),
    across(all_of(abx_cols), ~ mean(. == 1, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = all_of(abx_cols),
               names_to = "drug",
               values_to = "R_prop") %>%
  mutate(
    R_pct = R_prop * 100
  )

res_by_source
drug_order <- res_by_source %>%
  group_by(drug) %>%
  summarise(R_mean = mean(R_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(R_mean) %>%
  pull(drug)

res_by_source <- res_by_source %>%
  mutate(
    drug = factor(drug, levels = drug_order),
    source_group = factor(source_group,
                          levels = c("Farm", "Slaughter", "Retail", "Human"))
  )

source_cols <- c(
  "Farm"      = "#EF7A6D",
  "Slaughter" = "#B1CE46",
  "Retail"    = "#F1D77E",
  "Human"     = "#9DC3E7"
)

p_ast_facets <- ggplot(res_by_source,
                       aes(x = drug, y = R_pct, fill = source_group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", R_pct),
                y = R_pct + 3),
            size = 3) +
  facet_wrap(~ source_group, nrow = 1) +   
  scale_fill_manual(values = source_cols, guide = "none") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.10))  
  ) +
  coord_flip() +
  xlab("") +
  ylab("Resistance (%)") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )
p_ast_facets

## prevalence of ARGs
pkgs <- c("readxl","dplyr","tidyr","stringr","ggplot2","forcats","scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(scales)

file_path <- "fig2/Antibiotics resistance genes.xlsx"
dat <- read_xlsx(file_path)

source_col <- c("Source","sample","Sample","source")[c("Source","sample","Sample","source") %in% names(dat)][1]
if (is.na(source_col)) stop("无Source列")

dat <- dat %>%
  mutate(
    Source_raw  = as.character(.data[[source_col]]),
    Source_norm = str_to_lower(str_trim(Source_raw)),
    Source_norm = case_when(
      str_detect(Source_norm, "^transport") ~ "transport vehicles",
      TRUE ~ Source_norm
    ),
    source_group = case_when(
      Source_norm %in% c("ipf_pig","byf_pig") ~ "Farm",
      Source_norm %in% c("slaughterhouse","slaughterhouse wastewater","transport vehicles") ~ "Slaughter",
      Source_norm %in% c("pork","market wastewater") ~ "Retail",
      Source_norm %in% c("pig workers","diarrhea patients","healthy human") ~ "Human",
      TRUE ~ NA_character_
    ),
    source_group = factor(source_group, levels = c("Farm","Slaughter","Retail","Human"))
  )

meta_cols <- intersect(
  c(source_col, "Source_raw","Source_norm","source_group",
    "Strain","Isolate","ID","Genome","Accession","Species",
    "NUM_FOUND","num_found","n_found"),
  names(dat)
)
gene_cols <- setdiff(names(dat), meta_cols)

is_present <- function(x){
  x_chr <- toupper(trimws(as.character(x)))
  x_chr %in% c("1","TRUE","T","YES","Y","PRESENT","POS")
}

dat2 <- dat %>%
  filter(!is.na(source_group)) %>%
  mutate(across(all_of(gene_cols), is_present))

gene_prev <- dat2 %>%
  pivot_longer(cols = all_of(gene_cols), names_to = "gene", values_to = "present") %>%
  group_by(source_group, gene) %>%
  summarise(
    n_isolates = n(),
    pct = mean(present, na.rm = TRUE) * 100,
    .groups = "drop"
  )

min_any_pct <- 1   
gene_keep <- gene_prev %>%
  group_by(gene) %>%
  summarise(max_pct = max(pct, na.rm = TRUE), .groups = "drop") %>%
  filter(max_pct >= min_any_pct) %>%
  pull(gene)

gene_prev_f <- gene_prev %>% filter(gene %in% gene_keep)

classify_gene <- function(g){
  g <- as.character(g)
  
  if (str_detect(g, "^(aac\\(|aadA|ant\\(|aph\\(|armA$|rmt)")) return("Aminoglycosides")
  if (str_detect(g, "^bla")) return("Beta-lactams")
  if (str_detect(g, "^fosA")) return("Fosfomycin")
  if (str_detect(g, "^(lnu\\()")) return("Lincosamides")
  if (str_detect(g, "^(erm\\(|mph\\(|mef\\(|msr\\(|ere\\()")) return("Macrolides")
  if (str_detect(g, "^(cat$|catA|catB|cmlA|floR|cfr)")) return("Chloramphenicols")
  if (str_detect(g, "^(qnr|qepA|OqxA$|OqxB$|aac\\(6'\\)-Ib-cr$|crpP$|tmex|TOprJ)")) return("Quinolones")
  if (str_detect(g, "^ARR")) return("Rifamycins")
  if (str_detect(g, "^sul")) return("Sulfonamides")
  if (str_detect(g, "^tet")) return("Tetracyclines")
  if (str_detect(g, "^dfr")) return("Trimethoprims")
  if (str_detect(g, "^mcr")) return("Polymyxins")
  
  return("Other")
}

class_levels <- c(
  "Aminoglycosides","Beta-lactams","Fosfomycin","Lincosamides","Macrolides",
  "Chloramphenicols","Quinolones","Rifamycins","Sulfonamides",
  "Tetracyclines","Trimethoprims","Polymyxins","Other"
)

gene_class_df <- gene_prev_f %>%
  group_by(gene) %>%
  summarise(max_pct = max(pct, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    class = vapply(gene, classify_gene, character(1)),
    class = factor(class, levels = class_levels)
  )

gene_order <- gene_class_df %>%
  arrange(class, desc(max_pct), gene) %>%
  pull(gene)

gene_idx_df <- tibble(gene = gene_order, gene_idx = seq_along(gene_order)) %>%
  left_join(gene_class_df %>% select(gene, class), by = "gene")

plot_df <- gene_prev_f %>%
  left_join(gene_idx_df, by = "gene")

bracket_df <- gene_idx_df %>%
  group_by(class) %>%
  summarise(
    xmin = min(gene_idx),
    xmax = max(gene_idx),
    xmid = (xmin + xmax)/2,
    .groups = "drop"
  ) %>%
  filter(!is.na(class)) %>%
  arrange(class) %>%
  mutate(
    idx = row_number(),
    y0 = -4,
    y1 = -8,
    ytext = ifelse(idx %% 2 == 0, -22, -16)
  )

source_cols <- c(
  "Farm"      = "#EF7A6D",
  "Slaughter" = "#B1CE46",
  "Retail"    = "#F1D77E",
  "Human"     = "#9DC3E7"
)

p <- ggplot(plot_df, aes(x = gene_idx, y = pct, fill = source_group)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  facet_wrap(~ source_group, ncol = 1) +  
  scale_fill_manual(values = source_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(ylim = c(-30, 100), clip = "off") +
  scale_x_continuous(
    breaks = seq_along(gene_order),
    labels = gene_order,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(x = NULL, y = "Percent of isolates (%)") +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    plot.margin = margin(5.5, 5.5, 110, 5.5)
  ) +
  
geom_segment(
  data = bracket_df,
  aes(x = xmin - 0.4, xend = xmax + 0.4, y = y0, yend = y0),
  inherit.aes = FALSE, linewidth = 0.5
) +
  geom_segment(
    data = bracket_df,
    aes(x = xmin - 0.4, xend = xmin - 0.4, y = y0, yend = y1),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_segment(
    data = bracket_df,
    aes(x = xmax + 0.4, xend = xmax + 0.4, y = y0, yend = y1),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_text(
    data = bracket_df,
    aes(x = xmid, y = ytext, label = as.character(class)),
    inherit.aes = FALSE,
    angle = 35, vjust = 1, size = 4, fontface = "bold"
  )

print(p)


