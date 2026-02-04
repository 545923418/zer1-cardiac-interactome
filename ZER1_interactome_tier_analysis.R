

library(readxl)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(openxlsx)

# Load data
ip_ms_data <- read_excel("path/to/IP-MS_data.xlsx")
strict_genes_df <- read_excel("path/to/trans_output.xlsx", sheet = "Strict_Group")
target_gene_list <- strict_genes_df$`Gene Symbol`

print(paste("Loaded genes:", length(target_gene_list)))

# Prepare expression matrix
subset_expression <- data_GSE116250 %>% filter(Common_name %in% target_gene_list)

plot_matrix <- subset_expression %>% dplyr::select(-1, -2) %>% as.matrix()
rownames(plot_matrix) <- subset_expression$Common_name
plot_matrix_log <- log2(plot_matrix + 1)

# Extract NF samples only
nf_samples <- sample_info_GSE116250$Sample_ID[sample_info_GSE116250$Group == "Control"]
nf_matrix_raw <- plot_matrix[, nf_samples]
nf_matrix_log <- plot_matrix_log[, nf_samples]

# Filter by median and breadth
nf_medians <- apply(nf_matrix_raw, 1, median)
nf_breadth <- apply(nf_matrix_raw, 1, function(x) sum(x >= 1) / length(x))

keep_index <- (nf_medians > 2) & (nf_breadth >= 0.70)
expressed_matrix_log <- nf_matrix_log[keep_index, ]

print(paste("Filtered genes:", nrow(expressed_matrix_log)))

# Prepare filter plot data
filter_plot_df <- data.frame(
  Gene = rownames(nf_matrix_raw),
  Median_FPKM = nf_medians,
  Breadth = nf_breadth,
  Status = ifelse((nf_medians > 2) & (nf_breadth >= 0.7), "Passed", "Filtered"),
  stringsAsFactors = FALSE
)

count_kept <- sum(filter_plot_df$Status == "Passed")
count_dropped <- sum(filter_plot_df$Status == "Filtered")

highlight_list <- c("DVL2", "NFKB1", "PRKCA", "CAMKK2")
highlight_df <- subset(filter_plot_df, Gene %in% highlight_list)
kept_color <- "#CD3333"

# Plot filtering results
p_filter <- ggplot(filter_plot_df, aes(x = Median_FPKM, y = Breadth)) +
  geom_point(aes(color = Status), alpha = 0.5, size = 1.5) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_point(data = highlight_df, aes(x = Median_FPKM, y = Breadth),
             shape = 21, fill = kept_color, color = "black", size = 2, stroke = 0.6) +
  geom_label_repel(data = highlight_df, aes(x = Median_FPKM, y = Breadth, label = Gene),
                   color = "black", fill = "white", fontface = "bold", size = 3.5,
                   label.size = 0.25, box.padding = 0.6, max.overlaps = Inf) +
  scale_color_manual(values = c("Filtered" = "gray30", "Passed" = kept_color)) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 2, 5, 10, 50, 100)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1.05)) +
  theme_bw() +
  labs(title = paste0("Gene Filtering (NF only): Passed ", count_kept, " / Filtered ", count_dropped),
       subtitle = "Median > 2 & Breadth >= 0.7",
       x = "Median FPKM (log1p)", y = "Breadth")

print(p_filter)

# Pathway matching
target_genes <- rownames(expressed_matrix_log)

info_df <- AnnotationDbi::select(org.Hs.eg.db, keys = target_genes,
                                 columns = c("GENENAME", "GO"), keytype = "SYMBOL")

all_go_ids <- unique(info_df$GO[!is.na(info_df$GO)])
go_terms_df <- AnnotationDbi::select(GO.db, keys = all_go_ids, columns = "TERM", keytype = "GOID")

keywords <- list(
  "ECM_Collagen" = "extracellular matrix|extracellular structure|collagen|cell-matrix adhesion|wound healing|tissue remodeling",
  "Muscle_Structure" = "cardiac muscle hypertrophy|muscle tissue development|striated muscle|sarcomere|myofibril|actin filament|actin cytoskeleton",
  "Mech_Response" = "mechanical stimulus|response to pressure|response to stretch|shear stress",
  "Stress_Survival" = "apoptotic|oxidative stress|autophagy"
)

match_results <- data.frame(Gene = target_genes, Pathway_Label = "Other", stringsAsFactors = FALSE)

for(i in 1:nrow(match_results)) {
  g <- match_results$Gene[i]
  sub_info <- info_df[info_df$SYMBOL == g, ]
  all_text <- c(unique(sub_info$GENENAME), 
                go_terms_df$TERM[go_terms_df$GOID %in% unique(sub_info$GO)])
  all_text <- all_text[!is.na(all_text)]
  
  found_cats <- c()
  for(cat_name in names(keywords)) {
    if(any(grepl(keywords[[cat_name]], all_text, ignore.case = TRUE))) {
      found_cats <- c(found_cats, cat_name)
    }
  }
  
  if(length(found_cats) > 0) {
    match_results$Pathway_Label[i] <- paste(unique(found_cats), collapse = "; ")
  }
}

hits <- match_results[match_results$Pathway_Label != "Other", ]
print(paste("Pathway hits:", nrow(hits)))

write.csv(hits, "output/pathway_matches.csv", row.names = FALSE)

# Heatmap
if(nrow(hits) > 0) {
  filtered_matrix <- expressed_matrix_log[hits$Gene, ]
  
  row_annot <- data.frame(Pathway = hits$Pathway_Label, row.names = hits$Gene)
  annotation_col <- data.frame(Group = rep("Control", length(nf_samples)), row.names = nf_samples)
  
  pheatmap(filtered_matrix, annotation_col = annotation_col, annotation_row = row_annot,
           show_rownames = nrow(filtered_matrix) < 100, show_colnames = FALSE,
           scale = "row", cluster_cols = TRUE, cluster_rows = TRUE, fontsize_row = 8,
           main = "Heatmap: High Expr Genes in Pathways (NF only)", border_color = NA)
}

# Merge with IP-MS data
use_gene_col <- "Gene Symbol"
use_val_col <- "Log2FC"

plot_data <- match_results %>%
  inner_join(ip_ms_data, by = setNames(use_gene_col, "Gene")) %>%
  mutate(Tier_Group = ifelse(Pathway_Label != "Other", "Tier 1", "Tier 2"),
         Y_Value = .data[[use_val_col]])

print(paste("Merged genes:", nrow(plot_data)))

highlight_df_ipms <- plot_data %>% filter(Gene %in% highlight_list)
n_tier1 <- sum(plot_data$Tier_Group == "Tier 1")
n_tier2 <- sum(plot_data$Tier_Group == "Tier 2")

# Tier comparison plot
p_tier <- ggplot(plot_data, aes(x = Tier_Group, y = Y_Value)) +
  geom_violin(aes(fill = Tier_Group), trim = FALSE, alpha = 0.4, color = NA) +
  geom_jitter(aes(color = Tier_Group), width = 0.2, size = 1.5, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
  geom_point(data = highlight_df_ipms, aes(x = Tier_Group, y = Y_Value),
             color = "black", fill = "#CD3333", shape = 21, size = 3.5, stroke = 1) +
  geom_label_repel(data = highlight_df_ipms, aes(label = Gene), nudge_x = 0.4,
                   color = "black", fill = "white", fontface = "bold", box.padding = 0.5, max.overlaps = Inf) +
  scale_fill_manual(values = c("Tier 1" = "#CD3333", "Tier 2" = "gray70")) +
  scale_color_manual(values = c("Tier 1" = "#CD3333", "Tier 2" = "gray50")) +
  scale_x_discrete(labels = c(paste0("Tier 1\n(Core, n=", n_tier1, ")"),
                              paste0("Tier 2\n(Others, n=", n_tier2, ")"))) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 11, face = "bold"), panel.grid.major.x = element_blank()) +
  labs(title = "Validation: Core Pathway Genes show Higher IP-MS Enrichment",
       subtitle = "Tier 1 vs Tier 2", y = use_val_col, x = "")

print(p_tier)

# Tier1 ranking barplot
tier1_data <- plot_data %>% filter(Tier_Group == "Tier 1") %>% arrange(desc(Y_Value))
target_genes_bar <- c("DVL2", "NFKB1", "PRKCA", "CAMKK2")

bar_plot_data <- bind_rows(head(tier1_data, 20),
                           tier1_data %>% filter(Gene %in% target_genes_bar)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(Highlight = ifelse(Gene %in% target_genes_bar, "Target", "Top Ranked"))

p_rank <- ggplot(bar_plot_data, aes(x = reorder(Gene, Y_Value), y = Y_Value)) +
  geom_col(aes(fill = Highlight), width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("Target" = "#CD3333", "Top Ranked" = "gray85")) +
  geom_text(aes(label = round(Y_Value, 1)), hjust = -0.2, size = 3, color = "black") +
  theme_bw() +
  theme(legend.position = "none", axis.text.y = element_text(size = 10),
        panel.grid.major.y = element_blank()) +
  labs(title = "Top Enrichment Ranks in Tier 1", subtitle = "Highlighted genes validated by CoIP",
       x = "", y = use_val_col) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

print(p_rank)

# Statistical test
wilcox_result <- wilcox.test(Y_Value ~ Tier_Group, data = plot_data, alternative = "less")
print(paste("Wilcoxon p-value:", format(wilcox_result$p.value, scientific = TRUE)))
print(paste("Tier 1 median:", round(median(plot_data$Y_Value[plot_data$Tier_Group == "Tier 1"]), 2)))
print(paste("Tier 2 median:", round(median(plot_data$Y_Value[plot_data$Tier_Group == "Tier 2"]), 2)))

# Extract detailed pathway info
tier1_genes <- match_results %>% filter(Pathway_Label != "Other") %>% pull(Gene)
detailed_pathway_info <- data.frame(Gene = character(), Pathway_Category = character(),
                                    Specific_GO_Terms = character(), stringsAsFactors = FALSE)

for(g in tier1_genes) {
  sub_info <- info_df[info_df$SYMBOL == g, ]
  all_text <- c(unique(sub_info$GENENAME),
                go_terms_df$TERM[go_terms_df$GOID %in% unique(sub_info$GO)])
  all_text <- all_text[!is.na(all_text)]
  
  matched_categories <- c()
  matched_terms <- list()
  
  for(cat_name in names(keywords)) {
    pattern <- keywords[[cat_name]]
    hits <- grep(pattern, all_text, ignore.case = TRUE, value = TRUE)
    if(length(hits) > 0) {
      matched_categories <- c(matched_categories, cat_name)
      matched_terms[[cat_name]] <- unique(hits)
    }
  }
  
  if(length(matched_categories) > 0) {
    category_str <- paste(unique(matched_categories), collapse = "; ")
    all_specific_terms <- unique(unlist(matched_terms))
    specific_str <- if(length(all_specific_terms) > 3) {
      paste(c(all_specific_terms[1:3], "..."), collapse = " | ")
    } else {
      paste(all_specific_terms, collapse = " | ")
    }
    
    detailed_pathway_info <- rbind(detailed_pathway_info,
                                   data.frame(Gene = g, Pathway_Category = category_str,
                                              Specific_GO_Terms = specific_str, stringsAsFactors = FALSE))
  }
}

print(paste("Detailed pathway info extracted:", nrow(detailed_pathway_info)))

# Build supplementary tables
table_sx1 <- ip_ms_clean %>%
  dplyr::rename(Gene_symbol = all_of(use_gene_col), log2FC_ZER1_vs_IgG = all_of(use_val_col)) %>%
  left_join(seq_info_clean_full, by = c("Gene_symbol" = "Gene Symbol")) %>%
  left_join(heart_metrics_clean, by = c("Gene_symbol" = "Gene")) %>%
  left_join(go_tier_clean, by = c("Gene_symbol" = "Gene")) %>%
  left_join(detailed_pathway_info, by = c("Gene_symbol" = "Gene")) %>%
  mutate(Selected_for_CoIP = ifelse(Gene_symbol %in% c("DVL2", "NFKB1", "PRKCA", "CAMKK2"), "Y", "N"),
         CoIP_result = case_when(Gene_symbol == "DVL2" ~ "Positive",
                                 Gene_symbol %in% c("NFKB1", "PRKCA", "CAMKK2") ~ "Negative",
                                 TRUE ~ NA_character_),
         Tier = ifelse(is.na(Tier), "Filtered (Low Expr)", Tier)) %>%
  dplyr::select(Gene_symbol, any_of(c("Accession", "UniProt_ID", "Description")),
                log2FC_ZER1_vs_IgG, any_of(c("Adj. P-Value", "adj.P.Val", "P Value", "p_value")),
                AA, Met_cleavage_predicted, Nt1_in_GASTC, NF_median_TPM, NF_frac_TPM_ge_1,
                HeartExpressed_pass, Tier, Pathway_Category, Specific_GO_Terms,
                Selected_for_CoIP, CoIP_result)

table_sx2 <- table_sx1 %>%
  filter(Tier == "Tier 1") %>%
  arrange(desc(log2FC_ZER1_vs_IgG)) %>%
  mutate(Rank_within_Tier1 = row_number()) %>%
  dplyr::select(Rank_within_Tier1, Gene_symbol, log2FC_ZER1_vs_IgG,
                any_of(c("Adj. P-Value", "adj.P.Val")), AA, Met_cleavage_predicted,
                NF_median_TPM, NF_frac_TPM_ge_1, HeartExpressed_pass,
                Pathway_Category, Specific_GO_Terms, Selected_for_CoIP, any_of("CoIP_result"))

print(paste("Table Sx1:", nrow(table_sx1), "rows | Table Sx2:", nrow(table_sx2), "rows"))

# Export Excel
output_file <- "output/Supplementary_Tables_Final.xlsx"

if(file.exists(output_file)) {
  file.remove(output_file)
  Sys.sleep(1)
}

wb <- createWorkbook()
addWorksheet(wb, "Table Sx1_Full List")
writeData(wb, sheet = 1, table_sx1)
addWorksheet(wb, "Table Sx2_Tier1 Ranked")
writeData(wb, sheet = 2, table_sx2)
setColWidths(wb, sheet = 1, cols = 1:ncol(table_sx1), widths = "auto")
setColWidths(wb, sheet = 2, cols = 1:ncol(table_sx2), widths = "auto")
saveWorkbook(wb, output_file, overwrite = TRUE)

print(paste("Tables exported:", output_file))
print("Analysis complete")



