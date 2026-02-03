library(readxl)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(stringr)
library(writexl)

# Read filtered gene list from Python output
excel_path <- "trans_output.xlsx"

if(file.exists(excel_path)) {
  strict_genes_df <- read_excel(excel_path, sheet = "Strict_Group")
  target_gene_list <- strict_genes_df$`Gene Symbol`
  print(paste("Loaded Strict genes:", length(target_gene_list)))
} else {
  stop("File not found. Check path.")
}

# Prepare expression data (Strict genes only)
subset_expression <- data_GSE116250 %>%
  filter(Common_name %in% target_gene_list)

plot_matrix <- subset_expression %>%
  dplyr::select(-1, -2) %>%
  as.matrix()
rownames(plot_matrix) <- subset_expression$Common_name

plot_matrix_log <- log2(plot_matrix + 1)

# Extract Control samples only
nf_samples <- sample_info_GSE116250$Sample_ID[sample_info_GSE116250$Group == "Control"]
nf_matrix_raw <- plot_matrix[, nf_samples]
nf_matrix_log <- plot_matrix_log[, nf_samples]

print(paste("Analysis limited to Control samples:", length(nf_samples)))

# Apply filtering criteria
nf_medians <- apply(nf_matrix_raw, 1, median)
cond_median <- nf_medians > 2

nf_breadth <- apply(nf_matrix_raw, 1, function(x) sum(x >= 1) / length(x))
cond_breadth <- nf_breadth >= 0.70

keep_index <- cond_median & cond_breadth
expressed_matrix_log <- nf_matrix_log[keep_index, ]

print(paste("Genes passing filter (Median>2 & Breadth>70%):", nrow(expressed_matrix_log)))

# Prepare filter plot data
filter_plot_df <- data.frame(
  Gene = rownames(nf_matrix_raw),
  Median = nf_medians,
  Breadth = nf_breadth,
  Status = ifelse(keep_index, "Kept (Reserved)", "Filtered (Dropped)")
)

count_kept <- sum(keep_index)
count_dropped <- sum(!keep_index)

highlight_list <- c("DVL2", "NFKB1", "PRKCA", "CAMKK2")
highlight_df <- subset(filter_plot_df, Gene %in% highlight_list)
kept_color <- "#CD3333"

p_filter <- ggplot(filter_plot_df, aes(x = Median, y = Breadth)) +
  geom_point(aes(color = Status), alpha = 0.5, size = 1.5) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", size = 0.8) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "black", size = 0.8) +
  geom_point(data = highlight_df, aes(x = Median, y = Breadth),
             color = "black", fill = kept_color, shape = 21, size = 2, stroke = 0.6) +
  geom_label_repel(data = highlight_df, aes(label = Gene),
                   color = "black", fill = "white", fontface = "bold",
                   size = 3.5, label.size = 0.25, box.padding = 0.6, max.overlaps = Inf) +
  scale_color_manual(values = c("Filtered (Dropped)" = "gray85", "Kept (Reserved)" = kept_color)) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 2, 5, 10, 50, 100)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1.05)) +
  theme_bw() +
  labs(title = paste0("Gene Filtering: Kept ", count_kept, " / Dropped ", count_dropped),
       x = "Median FPKM (log1p)", y = "Breadth")

print(p_filter)

# GO pathway matching
target_genes <- rownames(expressed_matrix_log)

info_df <- AnnotationDbi::select(org.Hs.eg.db, keys = target_genes,
                                 columns = c("GENENAME", "GO"), keytype = "SYMBOL")

all_go_ids <- unique(info_df$GO[!is.na(info_df$GO)])
go_terms_df <- AnnotationDbi::select(GO.db, keys = all_go_ids, columns = c("TERM"), keytype = "GOID")

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
  all_text <- c(unique(sub_info$GENENAME), go_terms_df$TERM[go_terms_df$GOID %in% unique(sub_info$GO)])
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

check_genes <- c("DVL2", "NFKB1", "PRKCA", "CAMKK2")
for(g in check_genes) {
  if(g %in% hits$Gene) {
    print(paste("✓", g, "Hit:", hits$Pathway_Label[hits$Gene == g]))
  } else if (g %in% target_genes) {
    print(paste("○", g, "In list but marked as Other"))
  } else {
    print(paste("✗", g, "Filtered out"))
  }
}

write.csv(hits, "Target_Pathway_Matches.csv", row.names = FALSE)

# Heatmap
keep_indices <- match_results$Pathway_Label != "Other"
core_genes <- match_results$Gene[keep_indices]

if(length(core_genes) > 0) {
  filtered_matrix <- expressed_matrix_log[core_genes, ]
  row_annot <- data.frame(Pathway = match_results$Pathway_Label[keep_indices])
  rownames(row_annot) <- core_genes
  
  annotation_col <- data.frame(Group = sample_info_GSE116250$Group)
  rownames(annotation_col) <- sample_info_GSE116250$Sample_ID
  annotation_col <- annotation_col[nf_samples, , drop = FALSE]
  
  pheatmap(filtered_matrix, annotation_col = annotation_col, annotation_row = row_annot,
           show_rownames = nrow(filtered_matrix) < 100, show_colnames = FALSE,
           scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
           fontsize_row = 8, main = "Pathway Genes Heatmap", border_color = NA)
}

# IP-MS integration
ip_ms_data <- read_excel("IPMS_data.xlsx")

use_gene_col <- "Gene Symbol"
use_val_col <- "Log2FC"

plot_data <- match_results %>%
  inner_join(ip_ms_data, by = setNames(use_gene_col, "Gene")) %>%
  mutate(Tier_Group = ifelse(Pathway_Label != "Other", "Tier 1", "Tier 2"),
         Y_Value = .data[[use_val_col]])

if(nrow(plot_data) == 0) {
  stop("Merge failed. Check gene name format consistency.")
}

print(paste("Merged:", nrow(plot_data), "genes"))

highlight_df_tier <- plot_data %>% filter(Gene %in% highlight_list)
n_tier1 <- sum(plot_data$Tier_Group == "Tier 1")
n_tier2 <- sum(plot_data$Tier_Group == "Tier 2")

p_tier <- ggplot(plot_data, aes(x = Tier_Group, y = Y_Value)) +
  geom_violin(aes(fill = Tier_Group), trim = FALSE, alpha = 0.4, color = NA) +
  geom_jitter(aes(color = Tier_Group), width = 0.2, size = 1.5, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
  geom_point(data = highlight_df_tier, aes(x = Tier_Group, y = Y_Value),
             color = "black", fill = "#CD3333", shape = 21, size = 3.5, stroke = 1) +
  geom_label_repel(data = highlight_df_tier, aes(label = Gene),
                   nudge_x = 0.4, color = "black", fill = "white", fontface = "bold", box.padding = 0.5) +
  scale_fill_manual(values = c("Tier 1" = "#CD3333", "Tier 2" = "gray70")) +
  scale_color_manual(values = c("Tier 1" = "#CD3333", "Tier 2" = "gray50")) +
  scale_x_discrete(labels = c(paste0("Tier 1\n(Core, n=", n_tier1, ")"),
                              paste0("Tier 2\n(Others, n=", n_tier2, ")"))) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank()) +
  labs(title = "Core Pathway Genes show Higher IP-MS Enrichment",
       subtitle = "Tier 1 vs Tier 2", y = use_val_col, x = "")

print(p_tier)

# Tier 1 ranking
tier1_data <- plot_data %>% filter(Tier_Group == "Tier 1") %>% arrange(desc(Y_Value))
top20_df <- head(tier1_data, 20)
targets_df <- tier1_data %>% filter(Gene %in% highlight_list)
bar_plot_data <- bind_rows(top20_df, targets_df) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(Highlight = ifelse(Gene %in% highlight_list, "Target", "Top Ranked"))

p_rank <- ggplot(bar_plot_data, aes(x = reorder(Gene, Y_Value), y = Y_Value)) +
  geom_col(aes(fill = Highlight), width = 0.7, color = "black", size = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("Target" = "#CD3333", "Top Ranked" = "gray85")) +
  geom_text(aes(label = round(Y_Value, 1)), hjust = -0.2, size = 3, color = "black") +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major.y = element_blank()) +
  labs(title = "Top Enrichment Ranks in Tier 1", subtitle = "Highlighted genes validated by CoIP",
       x = "", y = use_val_col) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

print(p_rank)

# Generate supplementary tables
ip_ms_clean <- ip_ms_data %>% distinct(!!sym(use_gene_col), .keep_all = TRUE)

info_df_all <- AnnotationDbi::select(org.Hs.eg.db, keys = ip_ms_clean[[use_gene_col]],
                                     columns = c("GENENAME", "GO"), keytype = "SYMBOL")

all_go_ids_full <- unique(info_df_all$GO[!is.na(info_df_all$GO)])
go_terms_df_full <- AnnotationDbi::select(GO.db, keys = all_go_ids_full, columns = c("TERM"), keytype = "GOID")

match_results_full <- data.frame(Gene = ip_ms_clean[[use_gene_col]], Tier_Group = "Tier 2",
                                 Broad_Category = "Other", Specific_Terms = "", stringsAsFactors = FALSE)

for(i in 1:nrow(match_results_full)) {
  g <- match_results_full$Gene[i]
  sub_info <- info_df_all[info_df_all$SYMBOL == g, ]
  current_go_ids <- sub_info$GO[!is.na(sub_info$GO)]
  current_terms <- go_terms_df_full$TERM[go_terms_df_full$GOID %in% current_go_ids]
  gene_desc <- unique(sub_info$GENENAME)
  
  hit_categories <- c()
  matched_terms_all <- c()
  
  for(cat_name in names(keywords)) {
    pattern <- keywords[[cat_name]]
    hits_go <- grep(pattern, current_terms, value = TRUE, ignore.case = TRUE)
    hits_desc <- grep(pattern, gene_desc, value = TRUE, ignore.case = TRUE)
    
    if(length(hits_go) > 0 || length(hits_desc) > 0) {
      hit_categories <- c(hit_categories, cat_name)
      matched_terms_all <- c(matched_terms_all, hits_go)
    }
  }
  
  if(length(hit_categories) > 0) {
    match_results_full$Tier_Group[i] <- "Tier 1"
    match_results_full$Broad_Category[i] <- paste(unique(hit_categories), collapse = "; ")
    unique_terms <- unique(matched_terms_all)
    if(length(unique_terms) > 0) {
      match_results_full$Specific_Terms[i] <- paste(head(unique_terms, 3), collapse = "; ")
    } else {
      match_results_full$Specific_Terms[i] <- "Inferred from Gene Name/Description"
    }
  }
}

go_tier_clean <- match_results_full %>%
  dplyr::select(Gene, Tier_Group, Broad_Category, Specific_Terms) %>%
  distinct(Gene, .keep_all = TRUE)

heart_metrics_clean <- data.frame(Gene = names(nf_medians), NF_median_TPM = as.numeric(nf_medians),
                                  NF_frac_TPM_ge_1 = as.numeric(nf_breadth)) %>%
  group_by(Gene) %>%
  summarise(NF_median_TPM = max(NF_median_TPM), NF_frac_TPM_ge_1 = max(NF_frac_TPM_ge_1), .groups = "drop") %>%
  mutate(HeartExpressed_pass = (NF_median_TPM > 2) & (NF_frac_TPM_ge_1 >= 0.7))

if("Processed_Sequence" %in% colnames(strict_genes_df)) {
  seq_info_clean <- strict_genes_df %>%
    dplyr::select(`Gene Symbol`, Processed_Sequence) %>%
    distinct(`Gene Symbol`, .keep_all = TRUE) %>%
    mutate(Mature_Nt1 = substr(Processed_Sequence, 1, 1),
           Nt1_in_GASTC = Mature_Nt1 %in% c("G", "A", "S", "T", "C"))
} else {
  seq_info_clean <- data.frame(`Gene Symbol` = unique(ip_ms_clean[[use_gene_col]]), Mature_Nt1 = NA, check.names = FALSE)
}

table_sx1 <- ip_ms_clean %>%
  dplyr::rename(Gene_symbol = all_of(use_gene_col), log2FC_ZER1_vs_IgG = all_of(use_val_col)) %>%
  left_join(seq_info_clean, by = c("Gene_symbol" = "Gene Symbol")) %>%
  left_join(heart_metrics_clean, by = c("Gene_symbol" = "Gene")) %>%
  left_join(go_tier_clean, by = c("Gene_symbol" = "Gene")) %>%
  mutate(Selected_for_CoIP = ifelse(Gene_symbol %in% c("DVL2", "NFKB1", "PRKCA", "CAMKK2"), "Y", "N"),
         CoIP_result = case_when(Gene_symbol == "DVL2" ~ "Positive",
                                 Gene_symbol %in% c("NFKB1", "PRKCA", "CAMKK2") ~ "Negative",
                                 TRUE ~ NA_character_),
         Tier_Group = ifelse(is.na(Tier_Group), "Filtered", Tier_Group),
         Broad_Category = ifelse(is.na(Broad_Category), "", Broad_Category),
         Specific_Terms = ifelse(is.na(Specific_Terms), "", Specific_Terms)) %>%
  dplyr::select(Gene_symbol, any_of(c("Accession", "UniProt_ID", "Description")), log2FC_ZER1_vs_IgG,
                any_of(c("Adj. P-Value", "adj.P.Val")), Mature_Nt1, NF_median_TPM,
                Tier_Group, Broad_Category, Specific_Terms, Selected_for_CoIP, CoIP_result)

print(paste("Table Sx1 rows:", nrow(table_sx1)))

table_sx2 <- table_sx1 %>%
  filter(Tier_Group == "Tier 1") %>%
  arrange(desc(log2FC_ZER1_vs_IgG)) %>%
  mutate(Rank_within_Tier1 = row_number()) %>%
  dplyr::select(Rank_within_Tier1, Gene_symbol, log2FC_ZER1_vs_IgG,
                any_of(c("Adj. P-Value", "adj.P.Val")), Mature_Nt1, NF_median_TPM,
                Pathway_Category = Broad_Category, Matched_Specific_GO_Terms = Specific_Terms,
                Selected_for_CoIP)

print(paste("Table Sx2 rows:", nrow(table_sx2)))

output_file <- "Supplementary_Tables_Sx1_Sx2.xlsx"
writexl::write_xlsx(list("Table Sx1" = table_sx1, "Table Sx2_Tier1" = table_sx2), path = output_file)
print(paste("Tables exported to:", output_file))

