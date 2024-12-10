#Limma Differential Analysis 

#Script to run Differential Analysis using GSVA Hallmark Scores instead of gene expression. 
# Drugs ran: Olabarib
# Drugs to run: Carboplatin, Alisertib

#Libraries
library(limma)

#Load Data
counts <- read.csv("GSVA_Hallmark_13Osteosarcoma_CellLines.csv", row.names = 1)
counts1 <- t(counts)

#Design matrix
ola_des <- read.csv("input/Olaparib.csv")
ola_des1 <- ola_des[, c(1, 3)]
row.names(ola_des1) <- ola_des1$Cell.Line
ola_des1$Cell.Line <- NULL

counts2 <- counts1[, row.names(ola_des1)]

ola_des1$Res.vs.Sens <- factor(ola_des1$Res.vs.Sens)

design <- model.matrix(~ola_des1$Res.vs.Sens)
colnames(design)[1] <- "Intercept"  # For the intercept column
colnames(design)[2] <- "Res_vs_SensS"  # For the condition contrast column

contrast_matrix <- makeContrasts(Res_vs_SensS, levels = design)

fit <- lmFit(counts2, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, adjust.method = "BH", number = Inf)
print(results)

# Filter for significant results (e.g., adjusted p-value < 0.05)
significant_gene_sets <- results[results$adj.P.Val < 0.05, ]
print(significant_gene_sets)

# Volcano plot
volcanoplot(fit2, coef=1, main="Volcano plot for Treatment vs Control", xlab="LogFC", ylab="-log10(p-value)")

# Extract logFC for significant gene sets
sig_gene_sets <- rownames(significant_gene_sets)
heatmap_data <- gsva_results[sig_gene_sets, ]

# Create a heatmap of the significant gene sets
heatmap(heatmap_data, scale="row", main="Heatmap of Differentially Expressed Gene Sets")
