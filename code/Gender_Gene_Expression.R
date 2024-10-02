library(DESeq2)
library(limma)
library(edgeR)
library(dplyr)
library(stringr)


file_path <- "Combined_Expression_Gender.csv"
df <- read.csv(file_path)

df$Group <- trimws(df$Group)
df$Gender <- trimws(df$Gender)

df$Group <- factor(df$Group, levels = c("Control", "Disease"))
df$Gender <- factor(df$Gender)
df$Sample <- as.factor(df$Sample)

exprs <- as.matrix(df[, -(1:3)])  # Remove Sample, Gender, and Group columns
rownames(exprs) <- df$Sample
exprs <- t(exprs)

count_data <- round(exprs)  # DESeq2 requires integer counts

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = df, design = ~ Group)

dds <- DESeq(dds)

res_deseq <- results(dds)

res_deseq <- res_deseq[order(res_deseq$padj), ]

res_sig_deseq <- res_deseq[which(res_deseq$padj < 0.05), ]

y <- DGEList(counts = exprs)
y <- calcNormFactors(y)
norm_exprs <- cpm(y, log = TRUE)  # Log-transformed counts per million (CPM)

plotMDS(norm_exprs, labels = df$Group, col = as.numeric(df$Gender))

design <- model.matrix(~ 0 + Group:Gender, data = df)
colnames(design) <- make.names(colnames(design))

contrast.matrix <- makeContrasts(
  MaleControl_vs_MaleDisease = GroupDisease.GenderMale - GroupControl.GenderMale,
  FemaleControl_vs_FemaleDisease = GroupDisease.GenderFemale - GroupControl.GenderFemale,
  FemaleDisease_vs_MaleDisease = GroupDisease.GenderFemale - GroupDisease.GenderMale,
  FemaleControl_vs_MaleControl = GroupControl.GenderFemale - GroupControl.GenderMale,
  GenderInteraction = (GroupDisease.GenderMale - GroupControl.GenderMale) - (GroupDisease.GenderFemale - GroupControl.GenderFemale),
  levels = design
)

fit <- lmFit(norm_exprs, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract the top tables for each comparison
results_MaleControl_vs_MaleDisease <- topTable(fit2, coef = "MaleControl_vs_MaleDisease", adjust.method = "BH", number = Inf)
results_FemaleControl_vs_FemaleDisease <- topTable(fit2, coef = "FemaleControl_vs_FemaleDisease", adjust.method = "BH", number = Inf)
results_FemaleDisease_vs_MaleDisease <- topTable(fit2, coef = "FemaleDisease_vs_MaleDisease", adjust.method = "BH", number = Inf)
results_FemaleControl_vs_MaleControl <- topTable(fit2, coef = "FemaleControl_vs_MaleControl", adjust.method = "BH", number = Inf)
results_GenderInteraction <- topTable(fit2, coef = "GenderInteraction", adjust.method = "BH", number = Inf)

# Step 10: Filter significant results (adjusted p-value < 0.05) for limma/voom results
significant_MaleControl_vs_MaleDisease <- results_MaleControl_vs_MaleDisease %>% filter(adj.P.Val < 0.05)
significant_FemaleControl_vs_FemaleDisease <- results_FemaleControl_vs_FemaleDisease %>% filter(adj.P.Val < 0.05)
significant_FemaleDisease_vs_MaleDisease <- results_FemaleDisease_vs_MaleDisease %>% filter(adj.P.Val < 0.05)
significant_FemaleControl_vs_MaleControl <- results_FemaleControl_vs_MaleControl %>% filter(adj.P.Val < 0.05)
significant_GenderInteraction <- results_GenderInteraction %>% filter(adj.P.Val < 0.05)

