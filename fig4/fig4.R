library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(scales)
library(patchwork)
library(ggrepel)

base_dir <- "fig4"
file_arg <- file.path(base_dir, "ARG.xlsx")
file_vf  <- file.path(base_dir, "VF.xlsx")

if(!file.exists(file_arg) | !file.exists(file_vf)) stop("文件未找到")

arg_raw <- read_xlsx(file_arg)
vf_raw  <- read_xlsx(file_vf)

arg <- arg_raw %>% rename(ARG_count = NUM_FOUND)
vf  <- vf_raw  %>% rename(VF_count  = NUM_FOUND)

valid_strains <- intersect(arg$Strain, vf$Strain)

clean_meta <- function(df) {
  df %>%
    filter(Strain %in% valid_strains) %>%
    filter(!is.na(Source) & Source != "NA" & Source != "") %>%
    filter(tolower(Source) != "unknown") %>%
    filter(!is.na(MLST) & MLST != "NA" & MLST != "-" & MLST != 0)
}

arg <- clean_meta(arg)
vf  <- clean_meta(vf)

common_ids <- intersect(arg$Strain, vf$Strain)
arg <- arg %>% filter(Strain %in% common_ids) %>% arrange(Strain)
vf  <- vf  %>% filter(Strain %in% common_ids) %>% arrange(Strain)

cat("有效菌株数: ", length(common_ids), "\n")

meta_cols <- c("Strain","Collection date","Country","Source","Phylogroup","MLST")
arg_cols <- setdiff(names(arg), c(meta_cols, "ARG_count"))
vf_cols  <- setdiff(names(vf),  c(meta_cols, "VF_count"))

arg[arg_cols] <- lapply(arg[arg_cols], as.numeric)
vf[vf_cols]   <- lapply(vf[vf_cols],  as.numeric)

all_sources <- sort(unique(c(arg$Source, vf$Source)))
macaroon_cols <- c("#8ECFC9", "#FFBE7A", "#FA7F6F", "#82B0D2", "#BEB8DC", "#FF9AA2", 
                   "#E7DAD2", "#2878B5", "#B2EBF2", "#999999", "#F9E79F", "#AED6F1", 
                   "#F5B7B1", "#D2B4DE")
src_colors <- setNames(rep(macaroon_cols, length.out = length(all_sources)), all_sources)

arg$Source <- factor(arg$Source, levels = all_sources)
vf$Source  <- factor(vf$Source,  levels = all_sources)

## 风险逻辑定义
has_regex <- function(df, pat, pool) {
  cols <- intersect(pool, names(df))
  hit <- cols[str_detect(cols, regex(pat, ignore_case=TRUE))]
  if(length(hit)==0) return(rep(FALSE, nrow(df)))
  rowSums(as.matrix(df[, hit, drop=FALSE])==1, na.rm=TRUE) > 0
}
has_any <- function(df, cols) {
  cols <- intersect(cols, names(df))
  if(length(cols)==0) return(rep(FALSE, nrow(df)))
  rowSums(as.matrix(df[, cols, drop=FALSE])==1, na.rm=TRUE) > 0
}

# ARG
arg$has_carb <- has_regex(arg, "^bla(NDM|KPC|IMP|VIM)", arg_cols) | has_regex(arg, "^blaOXA-(48|181|232|485|488)$", arg_cols)
arg$has_mcr  <- has_regex(arg, "^mcr-", arg_cols)
arg$has_tmex <- has_regex(arg, "^tmex[CD]", arg_cols) | has_regex(arg, "^TOprJ", arg_cols)
arg$risk_critical <- arg$has_carb | arg$has_mcr | arg$has_tmex

arg$has_esbl_crit <- has_regex(arg, "^blaCTX-M", arg_cols) | has_regex(arg, "^blaSHV-12", arg_cols) | has_regex(arg, "^blaVEB", arg_cols)
arg$has_ampc <- has_regex(arg, "^bla(CMY|DHA|ACT|MIR)-", arg_cols)
arg$risk_high_esbl <- (arg$has_esbl_crit | arg$has_ampc) & !arg$risk_critical

q_arg75 <- quantile(arg$ARG_count, 0.75, na.rm=TRUE)
arg$high_res_load <- arg$ARG_count > q_arg75

arg$high_res <- arg$risk_critical | arg$risk_high_esbl | arg$high_res_load

# 标记类型
arg <- arg %>% mutate(
  Res_Mech_Type = case_when(
    risk_critical ~ "Critical (Carb/MCR/Tmex)",
    risk_high_esbl ~ "High Risk (ESBL/AmpC)",
    high_res_load ~ "High Load (>Q75)",
    TRUE ~ "Low/Baseline"
  )
)

# VF
vf$mk_pap <- has_regex(vf, "^pap[ACEFGHIJKX]", vf_cols)
vf$mk_sfa <- has_regex(vf, "^(sfa|foc)", vf_cols)
vf$mk_afa <- has_regex(vf, "^(afa|dra)", vf_cols)
vf$mk_iuc <- has_regex(vf, "^(iuc|iut)", vf_cols)
vf$mk_kps <- has_regex(vf, "^kps[MD]", vf_cols) | has_any(vf, c("kpsMII", "neuC"))
vf$patho_ExPEC <- (vf$mk_pap + vf$mk_sfa + vf$mk_afa + vf$mk_iuc + vf$mk_kps) >= 2

vf$has_hly <- has_regex(vf, "^hly[ABCD]", vf_cols)
vf$has_cnf <- has_regex(vf, "^cnf", vf_cols)
vf$patho_UPEC <- vf$patho_ExPEC & (vf$mk_pap | vf$has_hly | vf$has_cnf)
vf$has_stx <- has_regex(vf, "^stx", vf_cols)
vf$has_eae <- has_any(vf, c("eae", "eaeA", "eaeH", "eaeX"))
vf$patho_EHEC <- vf$has_stx & vf$has_eae
vf$patho_STEC <- vf$has_stx & !vf$has_eae
vf$patho_EPEC <- vf$has_eae & !vf$has_stx
vf$has_ibeA <- has_any(vf, "ibeA")
vf$has_K1 <- if("ECOK1" %in% vf_cols) vf$ECOK1==1 else (has_any(vf, "neuC") | has_regex(vf, "^kps", vf_cols))
vf$patho_NMEC <- vf$has_K1 & vf$has_ibeA
vf$patho_ETEC <- has_regex(vf, "^(elt|est|sta)", vf_cols)
vf$patho_EAEC <- has_any(vf, c("aggR", "aatA")) | has_regex(vf, "^aai", vf_cols)
vf$patho_EIEC <- has_any(vf, c("ipaH", "invE", "invG"))

vf <- vf %>% mutate(
  Pathotype_primary = case_when(
    patho_EHEC ~ "EHEC", patho_STEC ~ "STEC", patho_EIEC ~ "EIEC",
    patho_EAEC ~ "EAEC", patho_ETEC ~ "ETEC", patho_EPEC ~ "EPEC",
    patho_NMEC ~ "NMEC", patho_UPEC ~ "UPEC", patho_ExPEC ~ "ExPEC",
    TRUE ~ "Commensal"
  )
)

q_vf75 <- quantile(vf$VF_count, 0.75, na.rm=TRUE)
vf$high_vf_load <- vf$VF_count > q_vf75
vf$high_vf <- (vf$Pathotype_primary != "Commensal") | vf$high_vf_load

## Risk Category
dat_all <- arg %>%
  select(Strain, MLST, Source, Phylogroup, ARG_count, high_res, 
         risk_critical, risk_high_esbl, high_res_load, Res_Mech_Type,
         has_carb, has_mcr, has_tmex, has_esbl_crit, has_ampc) %>%
  inner_join(
    vf %>% select(Strain, VF_count, high_vf, high_vf_load, Pathotype_primary,
                  patho_ExPEC, patho_UPEC, patho_NMEC, patho_STEC, 
                  patho_EHEC, patho_EPEC, patho_ETEC, patho_EAEC, patho_EIEC),
    by = "Strain"
  )

dat_all <- dat_all %>%
  mutate(
    Risk_Category = case_when(
      high_res & high_vf ~ "Hybrid (High-Res + High-VF)",
      high_res & !high_vf ~ "High-Res Only",
      !high_res & high_vf ~ "High-VF Only",
      TRUE ~ "Low/Baseline"
    )
  )

write.csv(dat_all, file.path(base_dir, "Strain_Risk_Classification.csv"), row.names = FALSE)

## ST风险气泡图
st_stat <- dat_all %>%
  group_by(MLST) %>%
  summarise(
    n = n(),
    Med_ARG = median(ARG_count),
    Med_VF  = median(VF_count),
    Pct_HighRes = mean(high_res) * 100,
    Pct_HighVF  = mean(high_vf) * 100
  ) %>%
  filter(n >= 3) %>%
  arrange(desc(n)) %>%
  mutate(
    Risk_Class_ST = case_when(
      Pct_HighRes >= 50 & Pct_HighVF >= 50 ~ "Hybrid (Res+Vir)",
      Pct_HighRes >= 50 ~ "MDR Dominant",
      Pct_HighVF >= 50  ~ "Virulent Dominant",
      TRUE ~ "Low Risk"
    )
  )

st_colors <- c("Hybrid (Res+Vir)"="#D53E4F", "MDR Dominant"="#FDAE61", 
               "Virulent Dominant"="#3288BD", "Low Risk"="#E6F598")

p_bubble <- ggplot(st_stat, aes(x = Med_ARG, y = Med_VF)) +
  geom_hline(yintercept = q_vf75, linetype="dashed", color="grey70") +
  geom_vline(xintercept = q_arg75, linetype="dashed", color="grey70") +
  geom_point(aes(size = n, fill = Risk_Class_ST), shape=21, alpha=0.9) +
  geom_text_repel(data=subset(st_stat, Risk_Class_ST=="Hybrid (Res+Vir)" | n>50),
                  aes(label=paste0("ST", MLST)), size=4, max.overlaps=30) +
  scale_fill_manual(values=st_colors) +
  scale_size_continuous(range=c(2, 20)) +
  labs(title="Lineage Risk Stratification", subtitle=paste0("Ref: Q75 ARG=",q_arg75, ", VF=",q_vf75),
       x="Median ARG Count", y="Median VF Count", fill="ST Risk Profile") +
  theme_bw()

print(p_bubble)


## ST组合图
top3_st_ids <- st_stat$MLST[1:3]
sub_dat <- dat_all %>% 
  filter(MLST %in% top3_st_ids) %>% 
  mutate(MLST = factor(MLST, levels = top3_st_ids))

pie_dat <- sub_dat %>%
  count(MLST, Source) %>%
  group_by(MLST) %>%
  mutate(prop = n/sum(n))
p_pie <- ggplot(pie_dat, aes(x="", y=prop, fill=Source)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y") +
  facet_wrap(~MLST, nrow=1) +
  scale_fill_manual(values=src_colors) +
  geom_text(aes(label = ifelse(prop > 0, percent(prop, accuracy=1), "")), 
            position = position_stack(vjust = 0.5), size=3) +
  theme_void() +
  labs(title="A. Source Distribution") +
  theme(plot.title=element_text(face="bold"))
print(p_pie)

risk_bar_cols <- c(
  "Hybrid (High-Res + High-VF)" = "#C82423", 
  "High-Res Only"               = "#F46D43", 
  "High-VF Only"                = "#4575B4", 
  "Low/Baseline"                = "#E0F3F8"
)
bar_dat <- sub_dat %>%
  mutate(Risk_Category = factor(Risk_Category, levels=names(risk_bar_cols))) %>%
  count(MLST, Risk_Category) %>%
  group_by(MLST) %>%
  mutate(prop = n/sum(n))
p_bar <- ggplot(bar_dat, aes(x = MLST, y = prop, fill = Risk_Category)) +
  geom_col(width=0.7) +
  geom_text(aes(label = ifelse(prop > 0.02, percent(prop, accuracy=1), "")), 
            position = position_stack(vjust = 0.5), size=3) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = risk_bar_cols, name="Isolate Profile") +
  labs(y="Proportion", title="B. Risk Profile Composition") +
  theme_bw() +
  theme(plot.title=element_text(face="bold"), axis.title.x=element_blank())
print(p_bar)

heat_dat <- sub_dat %>%
  count(MLST, Source, Risk_Category) 
p_heat <- ggplot(heat_dat, aes(x = Source, y = Risk_Category, fill = n)) +
  geom_tile(color="white") +
  geom_text(aes(label = n, color = ifelse(n > max(n, na.rm=T)*0.6, "white", "black")), size=3, show.legend=FALSE) +
  scale_color_identity() + 
  facet_wrap(~MLST, nrow=1, scales="free_x") +
  scale_fill_gradient(low="grey90", high="#C82423", name="Count") +
  labs(title="C. Source-Risk Heatmap", x=NULL, y=NULL) +
  theme_bw() +
  theme(plot.title=element_text(face="bold"), 
        axis.text.x=element_text(angle=45, hjust=1),
        strip.background=element_rect(fill="grey95"))
print(p_heat)

layout <- "
AAAA
BBCC
"
final_plot <- p_pie / (p_bar | p_heat) + 
  plot_layout(heights=c(1, 1.4), design=layout) +
  plot_annotation(title="Detailed Characterization of Top 3 Lineages", 
                  theme=theme(plot.title=element_text(size=16, face="bold")))

print(final_plot)