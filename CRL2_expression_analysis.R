library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tibble)

# Human data (GSE116250)
data_GSE116250 <- read_tsv("GSE116250_rpkm.txt.gz", show_col_types = FALSE)

if (colnames(data_GSE116250)[2] != "Common_name") {
  colnames(data_GSE116250)[2] <- "Common_name"
}

human_sample_info <- data.frame(
  Sample_ID = colnames(data_GSE116250)[3:ncol(data_GSE116250)], 
  stringsAsFactors = FALSE
) %>%
  mutate(Group = case_when(
    grepl("^NF", Sample_ID) ~ "Control",     
    grepl("^DCM", Sample_ID) ~ "DCM",  
    TRUE ~ "Other"
  )) %>%
  filter(Group %in% c("Control", "DCM")) %>%
  mutate(Group = factor(Group, levels = c("Control", "DCM")))

human_target_genes <- c("ZER1", "CUL2", "TCEB1", "TCEB2", "RBX1", "ZYG11B")

human_plot_data_long <- data_GSE116250 %>%
  filter(Common_name %in% human_target_genes) %>%
  dplyr::select(Common_name, all_of(human_sample_info$Sample_ID)) %>% 
  pivot_longer(cols = -Common_name, names_to = "Sample_ID", values_to = "RPKM") %>%
  dplyr::rename(GeneName = Common_name) %>%
  left_join(human_sample_info, by = "Sample_ID") %>%
  mutate(Log2RPKM = log2(RPKM + 1))

p_human_genes <- ggplot(human_plot_data_long, aes(x = Group, y = Log2RPKM)) + 
  geom_boxplot(aes(fill = Group), width = 0.6, outlier.shape = NA, color = "black") +  
  geom_jitter(width = 0.2, size = 1.5, color = "black", alpha = 0.4) +
  stat_compare_means(comparisons = list(c("Control", "DCM")), method = "wilcox.test", 
                     label = "p.signif", size = 5, vjust = 0.5) +
  facet_wrap(~GeneName, scales = "free_y", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("Control" = "white", "DCM" = "gray60")) +
  theme_classic() + 
  labs(title = "Human: CRL2 Complex Genes", y = "Log2(RPKM + 1)", x = "") +
  theme(legend.position = "none", strip.text = element_text(face = "bold.italic"))

human_score_matrix <- human_plot_data_long %>%
  dplyr::select(GeneName, Sample_ID, Log2RPKM) %>%
  pivot_wider(names_from = Sample_ID, values_from = Log2RPKM) %>%
  column_to_rownames("GeneName")

human_z_score_matrix <- scale(t(human_score_matrix))
human_data_score_plot <- data.frame(
  Sample_ID = rownames(human_z_score_matrix),
  Signature_Score = rowMeans(human_z_score_matrix, na.rm = TRUE)
) %>% 
  left_join(human_sample_info, by = "Sample_ID")

p_human_score <- ggplot(human_data_score_plot, aes(x = Group, y = Signature_Score)) +
  geom_boxplot(aes(fill = Group), width = 0.5, outlier.shape = NA, color = "black") +  
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.4) +
  stat_compare_means(comparisons = list(c("Control", "DCM")), method = "wilcox.test", 
                     label = "p.signif", size = 8, vjust = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("Control" = "white", "DCM" = "gray60")) +
  theme_classic() + 
  labs(title = "Human: CRL2 Signature Score", y = "Z-score", x = "") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

human_corr_crl2 <- c("ZER1",......)
human_corr_markers <- c("NPPA", "NPPB", "MYH7")

human_raw_corr <- data_GSE116250 %>%
  filter(Common_name %in% c(human_corr_crl2, human_corr_markers)) %>%
  dplyr::select(Common_name, all_of(human_sample_info$Sample_ID))

human_mat_corr <- as.matrix(human_raw_corr[, -1])
rownames(human_mat_corr) <- human_raw_corr$Common_name
human_df_corr <- as.data.frame(t(log2(human_mat_corr + 1)))

human_valid_crl2 <- intersect(human_corr_crl2, colnames(human_df_corr))
human_valid_markers <- intersect(human_corr_markers, colnames(human_df_corr))

human_df_corr$CRL2_Score <- rowMeans(scale(human_df_corr[, human_valid_crl2]), na.rm = TRUE)
human_df_corr$Hypertrophy_Score <- rowMeans(scale(human_df_corr[, human_valid_markers]), na.rm = TRUE)

p_human_corr <- ggplot(human_df_corr, aes(x = CRL2_Score, y = Hypertrophy_Score)) +
  geom_point(fill = "gray60", shape = 21, color = "black", size = 3.5, stroke = 0.8, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", fill = "gray40", alpha = 0.2) +
  stat_cor(method = "pearson", label.x.npc = "center", label.y.npc = "top", size = 6) +
  theme_classic() +
  labs(title = "Human: CRL2-Hypertrophy Correlation", x = "CRL2 Score (Z-score)", 
       y = "Hypertrophy Score (Z-score)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

# Mouse data (GSE203083)
fpkm_data_203083 <- read_tsv("GSE203083_Raw_Gene_FPKM_Matrix.txt.gz", show_col_types = FALSE)
mouse_cols <- c("Sham_1", "Sham_2", "Sham_3", "TAC_1", "TAC_2", "TAC_3")

mouse_sample_info <- data.frame(
  Sample_ID = mouse_cols,
  Group = factor(c(rep("Sham", 3), rep("TAC", 3)), levels = c("Sham", "TAC")),
  stringsAsFactors = FALSE
)

mouse_expr_matrix <- fpkm_data_203083 %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  dplyr::select(gene_symbol, all_of(mouse_cols)) %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("gene_symbol")

mouse_target_genes <- c("Zer1", "Cul2", "Eloc", "Elob", "Rbx1", "Zyg11b")
mouse_valid_genes <- intersect(rownames(mouse_expr_matrix), mouse_target_genes)

mouse_plot_data_long <- as.data.frame(t(log2(mouse_expr_matrix[mouse_valid_genes, ] + 1))) %>%
  rownames_to_column("Sample_ID") %>%
  pivot_longer(cols = -Sample_ID, names_to = "Gene", values_to = "Log2FPKM") %>%
  left_join(mouse_sample_info, by = "Sample_ID")

p_mouse_genes <- ggplot(mouse_plot_data_long, aes(x = Group, y = Log2FPKM)) +
  geom_boxplot(aes(fill = Group), width = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.5) +
  stat_compare_means(comparisons = list(c("Sham", "TAC")), method = "t.test", 
                     label = "p.signif", size = 5, vjust = 0.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = 3) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("Sham" = "white", "TAC" = "gray60")) +
  theme_classic() +
  labs(title = "Mouse: CRL2 Complex Genes", y = "Log2(FPKM + 1)", x = "") +
  theme(legend.position = "none", strip.text = element_text(face = "bold.italic"))

mouse_z_mat <- scale(t(log2(mouse_expr_matrix[mouse_valid_genes, ] + 1)))
mouse_data_score <- data.frame(
  Sample_ID = rownames(mouse_z_mat),
  Signature_Score = rowMeans(mouse_z_mat, na.rm = TRUE)
) %>% 
  left_join(mouse_sample_info, by = "Sample_ID")

p_mouse_score <- ggplot(mouse_data_score, aes(x = Group, y = Signature_Score)) +
  geom_boxplot(aes(fill = Group), width = 0.5, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.5) +
  stat_compare_means(comparisons = list(c("Sham", "TAC")), method = "t.test", 
                     label = "p.signif", size = 8, vjust = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("Sham" = "white", "TAC" = "gray60")) +
  theme_classic() +
  labs(title = "Mouse: CRL2 Signature Score", y = "Z-score", x = "") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

mouse_corr_crl2 <- c("Zer1", ......)
mouse_corr_markers <- c("Nppa", "Nppb", "Myh7")

mouse_valid_crl2 <- intersect(mouse_corr_crl2, rownames(mouse_expr_matrix))
mouse_valid_markers <- intersect(mouse_corr_markers, rownames(mouse_expr_matrix))

mouse_genes_for_corr <- c(mouse_valid_crl2, mouse_valid_markers)
mouse_corr_df <- as.data.frame(t(log2(mouse_expr_matrix[mouse_genes_for_corr, ] + 1)))

if(length(mouse_valid_crl2) > 1) {
  mouse_corr_df$CRL2_Score <- rowMeans(scale(mouse_corr_df[, mouse_valid_crl2]), na.rm = TRUE)
} else {
  mouse_corr_df$CRL2_Score <- as.numeric(scale(mouse_corr_df[, mouse_valid_crl2]))
}

if(length(mouse_valid_markers) > 1) {
  mouse_corr_df$Hypertrophy_Score <- rowMeans(scale(mouse_corr_df[, mouse_valid_markers]), na.rm = TRUE)
} else {
  mouse_corr_df$Hypertrophy_Score <- as.numeric(scale(mouse_corr_df[, mouse_valid_markers]))
}

p_mouse_corr <- ggplot(mouse_corr_df, aes(x = CRL2_Score, y = Hypertrophy_Score)) +
  geom_point(fill = "gray60", shape = 21, color = "black", size = 4, stroke = 0.8, alpha = 0.9) +
  geom_smooth(method = "lm", color = "black", fill = "gray40", alpha = 0.2) +
  stat_cor(method = "pearson", label.x.npc = "center", label.y.npc = "top", size = 6) +
  theme_classic() +
  labs(title = "Mouse: CRL2-Hypertrophy Correlation", x = "CRL2 Score (Z-score)", 
       y = "Hypertrophy Score (Z-score)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

get_circle_pts <- function(cx, cy, r) {
  theta <- seq(0, 2 * pi, length.out = 100)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

circ_dcm <- get_circle_pts(-0.6, 0, 1.2)
circ_tac <- get_circle_pts(0.6, 0, 1.2)

p_venn <- ggplot() +
  geom_polygon(data = circ_dcm, aes(x, y), fill = "white", color = "black", size = 1, alpha = 0.5) +
  geom_polygon(data = circ_tac, aes(x, y), fill = "gray80", color = "black", size = 1, alpha = 0.5) +
  annotate("text", x = -1.0, y = 0, label = "CUL2\nELOC", size = 5, fontface = "italic") +
  annotate("text", x = 1.0, y = 0, label = "RBX1", size = 5, fontface = "italic") +
  annotate("text", x = 0, y = 0, label = "ZER1", size = 6, fontface = "bold.italic", color = "red") +
  annotate("text", x = -1.2, y = 1.4, label = "DCM (Human)", size = 6, fontface = "bold") +
  annotate("text", x = 1.2, y = 1.4, label = "TAC (Mouse)", size = 6, fontface = "bold") +
  coord_fixed() + theme_void() + 
  labs(title = "Differentially Expressed CRL2 Genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 20)))

print(p_human_genes)
print(p_human_score)
print(p_human_corr)
print(p_mouse_genes)
print(p_mouse_score)
print(p_mouse_corr)
print(p_venn)