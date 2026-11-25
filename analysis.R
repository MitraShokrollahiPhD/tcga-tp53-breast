library(cBioPortalData)
library(MultiAssayExperiment)
library(tidyverse)

brca <- cBioDataPack("brca_tcga_pan_can_atlas_2018", ask = FALSE)
names(experiments(brca))

mut <- brca[["mutations"]]
class(mut)
dim(mut)
names(assays(mut))
head(colnames(mut), 3)

expr <- brca[["mrna_seq_v2_rsem"]]
dim(expr)
head(colnames(expr), 3)

head(colnames(colData(brca)), 40)
library(cBioPortalData)
library(MultiAssayExperiment)
library(tidyverse)

brca <- cBioDataPack("brca_tcga_pan_can_atlas_2018", ask = FALSE)
mut <- brca[["mutations"]]
mut_gr <- unlist(as(mut, "GRangesList"))

length(mut_gr)
colnames(mcols(mut_gr))
head(names(mut_gr), 3)

expr <- brca[["mrna_seq_v2_rsem"]]
dim(expr)
head(colnames(expr), 3)
head(colnames(colData(brca)), 40)
mut <- brca[["mutations"]]
mut_gr <- unlist(as(mut, "GRangesList"))

length(mut_gr)
colnames(mcols(mut_gr))
head(names(mut_gr), 3)

expr <- brca[["mrna_seq_v2_rsem"]]
dim(expr)
head(colnames(expr), 3)
head(colnames(colData(brca)), 40)
mut_df <- data.frame(
  sample = sub("\\..*$", "", names(mut_gr)),
  gene   = sub("^[^.]*\\.", "", names(mut_gr)),
  vclass = as.character(mcols(mut_gr)$Variant_Classification)
)

nonsilent <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
               "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")

tp53_mut_samples <- mut_df %>%
  filter(gene == "TP53", vclass %in% nonsilent) %>%
  pull(sample) %>%
  unique()

length(tp53_mut_samples)
table(mut_df$vclass[mut_df$gene == "TP53"])
sum(colnames(expr) %in% colnames(mut))
head(rownames(expr), 5)
head(rownames(colData(brca)), 3)
tail(colnames(colData(brca)), 30)
table(colData(brca)$SUBTYPE)

cat("TP53 mutant samples:", length(tp53_mut_samples), "\n")
cat("Samples with expression AND mutation data:", sum(colnames(expr) %in% colnames(mut)), "\n")
cat("First 5 gene names in expression table:", head(rownames(expr), 5), "\n")
cat("TP53 mutation types:\n")
print(table(mut_df$vclass[mut_df$gene == "TP53"]))
cat("Survival-related columns:\n")
print(grep("OS_|DFS_|PFS_|DSS_|SURVIVAL", colnames(colData(brca)), value = TRUE))

# Getting the expression numbers as a plain matrix
expr_mat <- if (is.matrix(expr)) expr else assay(expr)

# Removing rows with blank gene names and duplicated names
expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & rownames(expr_mat) != "", ]
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]

# only tumors that have expression and mutation data both
keep_samples <- intersect(colnames(expr_mat), colnames(mut))
expr_mat <- expr_mat[, keep_samples]

# Log-transform
expr_log <- log2(expr_mat + 1)

# sample table: TP53 group + breast cancer subtype
clin <- as.data.frame(colData(brca))
sample_info <- data.frame(sample = keep_samples) %>%
  mutate(tp53 = ifelse(sample %in% tp53_mut_samples, "Mutant", "WT")) %>%
  left_join(dplyr::select(clin, SAMPLE_ID, SUBTYPE),
            by = c("sample" = "SAMPLE_ID"))

cat("Expression matrix (genes x tumors):", dim(expr_log), "\n")
cat("TP53 in expression table:", "TP53" %in% rownames(expr_log), "\n")
print(table(sample_info$tp53))
print(table(sample_info$tp53, sample_info$SUBTYPE, useNA = "ifany"))

# Dropping tumors with no subtype
info <- sample_info %>% filter(!is.na(SUBTYPE))
info$tp53 <- factor(info$tp53, levels = c("WT", "Mutant"))
info$SUBTYPE <- factor(info$SUBTYPE)
expr_use <- expr_log[, info$sample]

# Removing low expressed genes
keep <- rowMeans(expr_use > 1) >= 0.2
expr_use <- expr_use[keep, ]

# TP53 status, adjusting for subtype
design <- model.matrix(~ tp53 + SUBTYPE, data = info)
fit <- lmFit(expr_use, design)
fit <- eBayes(fit, trend = TRUE)
res <- topTable(fit, coef = "tp53Mutant", number = Inf, sort.by = "P")
res$gene <- rownames(res)

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/tp53_limma_results.csv", row.names = FALSE)

# Volcano plot
res$status <- case_when(
  res$adj.P.Val < 0.05 & res$logFC >  1 ~ "Up in TP53 mutant",
  res$adj.P.Val < 0.05 & res$logFC < -1 ~ "Down in TP53 mutant",
  TRUE ~ "Not significant"
)
volcano <- ggplot(res, aes(logFC, -log10(adj.P.Val), colour = status)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c("Up in TP53 mutant" = "firebrick",
                                 "Down in TP53 mutant" = "steelblue",
                                 "Not significant" = "grey70")) +
  labs(title = "TP53 mutant vs wild-type (adjusted for subtype)",
       x = "log2 fold change", y = "-log10 adjusted p-value") +
  theme_minimal()
print(volcano)
ggsave("results/volcano_tp53.png", volcano, width = 7, height = 5, dpi = 300)

cat("Genes tested:", nrow(res), "\n")
print(table(res$status))
print(head(res[, c("gene", "logFC", "adj.P.Val")], 10))

library(limma)
library(STRINGdb)
library(survival)
library(survminer)

# Setting up STRING with a cache folder
dir.create("string_cache", showWarnings = FALSE)
string_db <- STRINGdb$new(version = "12.0", species = 9606,
                          score_threshold = 400, network_type = "full",
                          input_directory = "string_cache")

# The top 100 genes and mapping them to STRING IDs
top_genes <- res %>%
  dplyr::slice_head(n = 100) %>%
  dplyr::select(gene, logFC, adj.P.Val)
mapped <- string_db$map(top_genes, "gene", removeUnmappedRows = TRUE)
cat("Genes mapped to STRING:", nrow(mapped), "of 100\n")

# Using all tested genes as the background for enrichment
all_mapped <- string_db$map(res[, "gene", drop = FALSE], "gene",
                            removeUnmappedRows = TRUE)

# survival table
surv_df <- clin %>%
  dplyr::select(SAMPLE_ID, OS_STATUS, OS_MONTHS) %>%
  dplyr::inner_join(info, by = c("SAMPLE_ID" = "sample")) %>%
  dplyr::mutate(
    time  = as.numeric(OS_MONTHS),
    event = ifelse(grepl("^1", OS_STATUS), 1, 0)
  ) %>%
  dplyr::filter(!is.na(time), time > 0)

cat("Tumors in survival analysis:", nrow(surv_df), "\n")
print(table(surv_df$OS_STATUS))

# Kaplan-Meier curve (unadjusted)
fit_km <- survfit(Surv(time, event) ~ tp53, data = surv_df)
km_plot <- ggsurvplot(fit_km, data = surv_df, pval = TRUE, risk.table = TRUE,
                      palette = c("steelblue", "firebrick"),
                      legend.labs = c("WT", "TP53 mutant"),
                      xlab = "Months", ylab = "Overall survival")
print(km_plot)
png("results/km_tp53.png", width = 2000, height = 1800, res = 300)
print(km_plot)
dev.off()

# Cox model adjusted for subtype
cox <- coxph(Surv(time, event) ~ tp53 + SUBTYPE, data = surv_df)
print(summary(cox))
string_db$set_background(all_mapped$STRING_id)

# Plot and save the network (summary text turned off)
string_db$plot_network(mapped$STRING_id, add_summary = FALSE)
png("results/string_network_tp53.png", width = 2000, height = 2000, res = 300)
string_db$plot_network(mapped$STRING_id, add_summary = FALSE)
dev.off()

# Functional enrichment
enrich <- string_db$get_enrichment(mapped$STRING_id)
enrich <- enrich[order(enrich$fdr), ]
write.csv(enrich, "results/string_enrichment.csv", row.names = FALSE)
print(head(enrich[, c("category", "description", "number_of_genes", "fdr")], 10))

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")