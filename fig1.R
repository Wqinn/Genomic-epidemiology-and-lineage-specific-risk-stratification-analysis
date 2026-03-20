## tigecycline resistance positive rate（%）
pkgs <- c("ggplot2", "dplyr", "tibble", "binom")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggplot2)
library(dplyr)
library(tibble)
library(binom)

df <- tribble(
  ~sample,                     ~sample_num, ~tetX4_ecoli,
  "IPF_pig",                          443,        383,
  "BYF_pig",                          408,        291,
  "slaughterhouse",                   315,        202,
  "pork",                             309,        215,
  "diarrhea patients",                335,         25,
  "healthy human",                     97,          2,
  "pig workers",                      115,         32,
  "slaughterhouse wastewater",         80,         52,
  "market wastewater",                100,         80,
  "transport vehicles",                14,          6
)

df <- df %>%
  mutate(
    group = case_when(
      sample %in% c("IPF_pig", "BYF_pig") ~ "Farm",
      sample %in% c("slaughterhouse", "slaughterhouse wastewater", "transport vehicles") ~ "Slaughter",
      sample %in% c("pork", "market wastewater") ~ "Retail",
      sample %in% c("pig workers", "diarrhea patients", "healthy human") ~ "Human",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("Farm", "Slaughter", "Retail", "Human"))
  )

ci <- binom::binom.confint(x = df$tetX4_ecoli, n = df$sample_num, methods = "wilson")

df <- df %>%
  mutate(
    positive = ci$mean * 100,
    lower = ci$lower * 100,
    upper = ci$upper * 100
  )

sample_order <- c(
  "IPF_pig", "BYF_pig",
  "slaughterhouse", "slaughterhouse wastewater", "transport vehicles",
  "pork", "market wastewater",
  "pig workers", "diarrhea patients", "healthy human"
)
df$sample <- factor(df$sample, levels = sample_order)

group_colors <- c(
  "Farm"      = "#EF7A6D",
  "Slaughter" = "#B1CE46",
  "Retail"    = "#F1D77E",
  "Human"     = "#9DC3E7"
)

p <- ggplot(df, aes(x = sample, y = positive, fill = group)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.6) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(
    name = "tigecycline resistance positive rate (%)",
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = NULL) +
  guides(fill = "none") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.spacing.x = unit(0.8, "lines"),
    panel.grid.major.y = element_line(linewidth = 0.3, linetype = "solid"),
    panel.grid.minor.y = element_line(linewidth = 0.2, linetype = "dotted"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
print(p)


## tet(X4)-positive Ecoli(%)
pkgs <- c("ggplot2", "dplyr", "tibble", "binom")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(ggplot2)
library(dplyr)
library(tibble)
library(binom)

df <- tribble(
  ~sample,                     ~sample_num, ~tetX4_ecoli,
  "IPF_pig",                          443,        318,
  "BYF_pig",                          408,        179,
  "slaughterhouse",                   315,        111,
  "pork",                             309,         95,
  "diarrhea patients",                335,          4,
  "healthy human",                     97,          1,
  "pig workers",                      115,         17,
  "slaughterhouse wastewater",         80,         37,
  "market wastewater",                100,         26,
  "transport vehicles",                14,          2
)

df <- df %>%
  mutate(
    group = case_when(
      sample %in% c("IPF_pig", "BYF_pig") ~ "Farm",
      sample %in% c("slaughterhouse", "slaughterhouse wastewater", "transport vehicles") ~ "Slaughter",
      sample %in% c("pork", "market wastewater") ~ "Retail",
      sample %in% c("pig workers", "diarrhea patients", "healthy human") ~ "Human",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("Farm", "Slaughter", "Retail", "Human"))
  )

ci <- binom::binom.confint(x = df$tetX4_ecoli, n = df$sample_num, methods = "wilson")

df <- df %>%
  mutate(
    positive = ci$mean * 100,
    lower = ci$lower * 100,
    upper = ci$upper * 100
  )

sample_order <- c(
  "IPF_pig", "BYF_pig",
  "slaughterhouse", "slaughterhouse wastewater", "transport vehicles",
  "pork", "market wastewater",
  "pig workers", "diarrhea patients", "healthy human"
)
df$sample <- factor(df$sample, levels = sample_order)

group_colors <- c(
  "Farm"      = "#EF7A6D",
  "Slaughter" = "#B1CE46",
  "Retail"    = "#F1D77E",
  "Human"     = "#9DC3E7"
)

p <- ggplot(df, aes(x = sample, y = positive, fill = group)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.6) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(
    name = "tet(X4)-positive E. coli (%)",
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = NULL) +
  guides(fill = "none") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.spacing.x = unit(0.8, "lines"),
    panel.grid.major.y = element_line(linewidth = 0.3, linetype = "solid"),
    panel.grid.minor.y = element_line(linewidth = 0.2, linetype = "dotted"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
print(p)
