#load dependencies

library(tidyverse)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(limma)

#Data Preprocssing for cancer dataset----

lines <- readLines("C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//GSE134359_series_matrix.txt")

start_index <- which(lines == "!series_matrix_table_begin") + 1
end_index <- which(lines == "!series_matrix_table_end") - 1

data_lines <- lines[start_index:end_index]

writeLines(data_lines, "C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//GSE134359_series_matrix.processed.txt")


raw_counts <- read.delim("C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//GSE134359_series_matrix.processed.txt", header = TRUE, row.names = 1, sep = "\t")


normal_tissue <- colnames(raw_counts[, which(colnames(raw_counts)=="GSM3943700"):which(colnames(raw_counts) == "GSM3943711")])

cancer <- colnames(raw_counts[, which(colnames(raw_counts)=="GSM3943712"):which(colnames(raw_counts) == "GSM3943741")])


subset <- c(normal_tissue, cancer)

raw_count_subtypes <- raw_counts[, subset]



#Create metadata

metaData <- data.frame(
  condition = factor(c(rep("Normal_tissue", 12), rep("Cancer", 30)))
  

)

rownames(metaData) <- colnames(raw_count_subtypes)



#Check if rows of metadata is in the same order as the columns of the raw counts
all(rownames(metaData)==colnames(raw_count_subtypes))

all_integers <- all(apply(raw_count_subtypes, 2, function(x) all(x == as.integer(x))))

if (all_integers) {
  cat("All values in raw_counts_subset are integers.\n")
} else {
  cat("Some values in raw_counts_subset are not integers.\n")
}


# Convert the raw counts matrix to integers
raw_counts_subtypes <- round(raw_count_subtypes)
head(raw_counts_subtypes)




# Normalize between arrays (quantile normalization)
exprs_norm <- normalizeBetweenArrays(raw_count_subtypes, method = "quantile")

# Build the design matrix
design <- model.matrix(~0 + condition, data = metaData)
colnames(design) <- levels(metaData$condition)



fit <- lmFit(exprs_norm, design)

# Compare Cancer vs Normal
contrast_matrix <- makeContrasts(Cancer - Normal_tissue, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# View top differentially expressed genes
topTable(fit2, adjust = "fdr", number = 20)


library(ggplot2)

# Extract top table
deg_table <- topTable(fit2, number = Inf, adjust = "fdr")

# Add a column to flag significance
deg_table$significant <- with(deg_table, adj.P.Val < 0.05 & abs(logFC) > 1)

# Basic volcano plot
ggplot(deg_table, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("grey", "firebrick")) +
  labs(title = "Volcano Plot: Cancer vs Normal",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Install annotation package (only once)
if (!require("hta20transcriptcluster.db")) {
  BiocManager::install("hta20transcriptcluster.db")
}

library(hta20transcriptcluster.db)

# Your table of results (topTable)
deg_table$Symbol <- mapIds(hta20transcriptcluster.db,
                           keys = rownames(deg_table),
                           column = "SYMBOL",
                           keytype = "PROBEID",
                           multiVals = "first")

# Assuming deg_table is your full results table with logFC and adj.P.Val
upregulated <- rownames(deg_table[deg_table$logFC > 1 & deg_table$adj.P.Val < 0.05, ])
downregulated <- rownames(deg_table[deg_table$logFC < -1 & deg_table$adj.P.Val < 0.05, ])

install.packages("VennDiagram")


library(VennDiagram)

venn.plot <- venn.diagram(
  x = list("Upregulated" = upregulated, "Downregulated" = downregulated),
  filename = NULL, 
  fill = c("firebrick", "steelblue"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0,
  main = "DEG Overlap: Cancer vs Normal"
)

grid.newpage()
grid.draw(venn.plot)



if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}


library(clusterProfiler)
library(org.Hs.eg.db)

# Filter for significant genes (adjust as needed)
sig_genes <- deg_table$Symbol[deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) > 1]

# Filter upregulated genes (padj < 0.05 and logFC > 1)
upregulated <- deg_table[deg_table$adj.P.Val < 0.05 & deg_table$logFC > 1, ]

# Sort by descending logFC
upregulated_sorted <- upregulated[order(-upregulated$logFC), ]

# Extract the top 500

# Remove NAs
sig_genes <- na.omit(top500_genes)

top500_genes <- head(upregulated_sorted[,c("Symbol", "logFC")], 500)

#write into csv
write.csv(sig_genes[, c('Symbol', 'logFC'), ], "C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//deg_genes.csv", row.names = FALSE)
 


write.table(top500_genes, file = "C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//deg_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(sig_genes, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene          = entrez_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",  # Biological Process
                pAdjustMethod = "fdr",
                qvalueCutoff  = 0.05)

head(ego)

barplot(ego, showCategory = 10, title = "GO Biological Processes")
dotplot(ego, showCategory = 10, title = "GO Enrichment: Cancer DEGs")


ego_MF <- enrichGO(gene          = entrez_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "MF",  # Biological Process
                pAdjustMethod = "fdr",
                qvalueCutoff  = 0.05)

barplot(ego_MF, showCategory = 10, title = "Molecular Functions")


#Data preprocessing of diabetes data

lines_diabetes <- readLines("C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//GSE20966_series_matrix.txt")

start_index_db <- which(lines_diabetes == "!series_matrix_table_begin") + 1
end_index_db <- which(lines_diabetes == "!series_matrix_table_end") - 1

data_lines_db <- lines_diabetes[start_index_db:end_index_db]

writeLines(data_lines_db, "C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//GSE20966_series_matrix.processed.txt")


raw_counts_db <- read.delim("C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//GSE20966_series_matrix.processed.txt", header = TRUE, row.names = 1, sep = "\t")

head(raw_counts_db)

#Create metadata

metaData_db <- data.frame(
  condition = factor(c(rep("Nondiabetic", 10), rep("t2diabetic", 10)))
  
  
)

rownames(metaData_db) <- colnames(raw_counts_db)



#Check if rows of metadata is in the same order as the columns of the raw counts
all(rownames(metaData_db)==colnames(raw_counts_db))

all_integers_db <- all(apply(raw_counts_db, 2, function(x) all(x == as.integer(x))))

if (all_integers_db) {
  cat("All values in raw_counts_subset are integers.\n")
} else {
  cat("Some values in raw_counts_subset are not integers.\n")
}


# Convert the raw counts matrix to integers
raw_counts_db <- round(raw_counts_db)
head(raw_counts_db)

# Normalize between arrays (quantile normalization)
exprs_norm_db <- normalizeBetweenArrays(raw_counts_db, method = "quantile")

# Build the design matrix
design <- model.matrix(~0 + condition, data = metaData_db)
colnames(design) <- levels(metaData_db$condition)



fit_db <- lmFit(exprs_norm_db, design)

# Compare t2d vs Nont2d
contrast_matrix_db <- makeContrasts(t2diabetic - Nondiabetic, levels = design)
fit2_db <- contrasts.fit(fit_db, contrast_matrix_db)
fit2_db <- eBayes(fit2_db)

# View top differentially expressed genes
topTable(fit2_db, adjust = "fdr", number = 20)

# Extract top table
deg_table_db <- topTable(fit2_db, number = Inf, adjust = "fdr")

# Add a column to flag significance
deg_table_db$significant <- with(deg_table_db, adj.P.Val < 0.05 & abs(logFC) > 1)

# Basic volcano plot
ggplot(deg_table_db, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("grey", "firebrick")) +
  labs(title = "Volcano Plot: T2D vs NonT2D",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")


library(hgu133x3p.db)


deg_table_db$Symbol <- mapIds(hta20transcriptcluster.db,
                           keys = rownames(deg_table_db),
                           column = "SYMBOL",
                           keytype = "PROBEID",
                           multiVals = "first")


head(keys(hta20transcriptcluster.db, keytype = "PROBEID"))

if (!requireNamespace("hgu133x3p.db", quietly = TRUE)) {
  BiocManager::install("hgu133x3p.db")
}
library(hgu133x3p.db)

version




#Export gene list for PPI analysis
write.table(gene_list, file = "C://Users//ghagi//Desktop//Differential Gene Expression Cancer and Diabetes//deg_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

