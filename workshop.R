# Installation of R Packages for Transcriptomic Analysis

## Dependencies

# Install 'BiocManager' if it's not already installed
# BiocManager is used for installing and managing Bioconductor packages.
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Check if 'vsn' package is available, if not, install it using BiocManager
# 'vsn' is used for variance stabilization and normalization of high-throughput data.
if (!requireNamespace("vsn", quietly = TRUE)) {
  BiocManager::install("vsn")
}

# Check if 'apeglm' package is available, if not, install it using BiocManager
# 'apeglm' is used for adaptive shrinkage of effect sizes in DESeq2 analyses.
if (!requireNamespace("apeglm", quietly = TRUE)) {
  BiocManager::install("apeglm")
}

# Check if 'ggplot2' package is available, if not, install it
# 'ggplot2' is used for data visualization and creating complex plots from data
# in a data frame.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Check if 'pheatmap' package is available, if not, install it
# 'pheatmap' is used for creating heatmaps, including clustering and annotation features.
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# Check if 'circlize' package is available, if not, install it
# 'circlize' is used for creating circular visualizations, useful for genomic data.
if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}

# Check if 'RColorBrewer' package is available, if not, install it
# 'RColorBrewer' provides color palettes for visualizing data.
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

## Dependences from Bioconductor

# Install and load BiocManager (if not already installed)
# BiocManager is used for installing and managing Bioconductor packages.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install and load the 'enrichplot' package
# 'enrichplot' is used for visualization of functional enrichment results.
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}

# Install and load the 'DOSE' package
# 'DOSE' is used for disease ontology semantic and enrichment analysis.
if (!requireNamespace("DOSE", quietly = TRUE)) {
  BiocManager::install("DOSE", force = TRUE)
}

# Install and load the 'clusterProfiler' package
# 'clusterProfiler' is used for statistical analysis and visualization of
# functional profiles for genes and gene clusters.
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

# Install and load the 'org.Dr.eg.db' package
# 'org.Dr.eg.db' provides genome-wide annotation for Zebrafish, including mappings
# between gene IDs and biological data.
if (!requireNamespace("org.Dr.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Dr.eg.db")
}

# Install and load the 'ComplexHeatmap' package
# 'ComplexHeatmap' is used for creating complex and highly customizable heatmaps
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}

# Check if 'DESeq2' package is available, if not, install it using BiocManager
# 'DESeq2' is used for differential gene expression analysis based on the negative
# binomial distribution.
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

# Import necessary libraries
# DESeq2 is used for differential expression analysis
# ggplot2 is used for creating various plots and visualizations
# pheatmap is used for creating heatmaps to visualize data
# pasilla provides an example RNA-seq dataset for demonstration
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(pasilla)

########################################################################################################
###################### DIFFERENTIAL EXPRESSION GENES (DEG) ANALYSIS USING DESeq2 #######################
########################################################################################################

# Define the RNA-seq data directory
# This specifies the location of the directory containing RNA-seq count data in TSV format
directory <- "/Users/deandradesilvae/Documents/edson/lab_meet/workshop/vcorimb/tsv/"


# Read the sample table from "sample.txt" and convert 'condition' to factor
# 'sampleTable' contains information about the samples, such as conditions
# The 'condition' column is converted to a factor to be used in the DESeq2 analysis
sampleTable <- read.table("/Users/deandradesilvae/Documents/edson/lab_meet/workshop/vcorimb/tsv/sample.txt", header = TRUE)
sampleTable$condition <- factor(sampleTable$condition)

sampleTable

# Create a DESeqDataSet from HTSeqCount data
# DESeqDataSetFromHTSeqCount creates a DESeq2 dataset from HTSeq count files
# The 'design' formula specifies the experimental design, in this case, the condition of the samples
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

# Perform the DESeq analysis
# DESeq performs differential analysis, including normalization and estimation of dispersions
ddsAB <- DESeq(ddsHTSeq)

# Get results for the comparison between conditions (LG vs. FP) and save to a file
# 'results' extracts the results of the differential expression analysis for the specified contrast
# The results are saved to a text file for further examination
resAB <- results(ddsAB, contrast = c("condition", "LG", "FP"))

head(resAB)

write.table(resAB, file = "resultados_diff.txt")

## Fold Change explanation

head(resAB)

# Transformations for data visualization
# Variance stabilizing transformation (VST) and regularized log transformation (Rlog) are
# applied for visualization
# These transformations make the data more suitable for plotting by stabilizing variance
vst <- vst(ddsAB, blind=FALSE)
rld <- rlog(ddsAB, blind=FALSE)

# Plot a PCA using the transformed data
# PCA (Principal Component Analysis) is used to visualize the overall patterns and relationships
# between samples
plotPCA(vst, intgroup=c("condition", "sizeFactor"), ntop = 3400)

# Perform shrinkage of the results matrix for MA analysis
# lfcShrink applies shrinkage to the log fold changes to improve visualization and interpretation
# Apeglm method is used here with a threshold for fold change
resApeT <- lfcShrink(ddsAB, coef=2, type="apeglm", lfcThreshold=1)

# Visualize an MA plot
# MA plot is used to visualize the relationship between fold change and average expression level
# Points are plotted with color indicating their significance
plotMA(resApeT, ylim=c(-3,3), cex=.8)

abline(h=c(-1,1), col="dodgerblue", lwd=2)

## Volcano Plot in Gene Expression Analysis:


########################################################################################################
################ USING THE VOLCANO PLOT TO REPRESENT DIFFERENTIAL GENE EXPRESSION (DGE) ################
########################################################################################################

# Volcano Plot
# Read the differential expression results from the "resultados_diff.txt" file
# The data is read with columns separated by spaces, and the decimal point as '.'
data_diff <- read.table("resultados_diff.txt", header = TRUE, sep = " ", dec = ".")
data_diff$gene_id <- rownames(data_diff)
colnames(data_diff) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_id")

head(data_diff)

# Plot the Volcano Plot
# The x-axis represents the log2 fold change, and the y-axis represents the -log10 adjusted p-value
plot(data_diff$log2FoldChange, -log10(data_diff$padj),
     pch=20,  # Point character
     xlab="Log2FC",  # X-axis label
     ylab="-log10(FDR)",  # Y-axis label
     main="12hai",  # Title of the plot
     cex=1,  # Size of the points
     xlim=c(-20,20),  # **Limits for the x-axis**
     ylim=c(0,110),  # Limits for the y-axis
     col = "#d9d9d9",  # Color for the points
     abline(h=c(1.3), v=c(-1.5,1.5), lty=3, lwd=3))  # Add horizontal and vertical lines

# Select and plot up-regulated points
# Up-regulated genes are those with an adjusted p-value < 0.05 and log2 fold change > 1.5
id <- subset(data_diff, padj < 0.05 & log2FoldChange > 1.5)

head(id)

write.table(id, file = "up_regulated.txt", sep = "\t", row.names = FALSE)  # Save up-regulated genes to a file
points(id$log2FoldChange, -log10(id$padj), pch=20, col="#bf523e", cex=1.2)  # Plot up-regulated points with a different color

# Select and plot down-regulated points
# Down-regulated genes are those with an adjusted p-value < 0.05 and log2 fold change < -1.5
id2 <- subset(data_diff, padj < 0.05 & log2FoldChange < -1.5)

head(id2)

write.table(id2, file = "down_regulated.txt", sep = "\t", row.names = FALSE)  # Save down-regulated genes to a file
points(id2$log2FoldChange, -log10(id2$padj), pch=20, col="#4d4dff", cex=1.2)  # Plot down-regulated points with a different color

# Volcano Plot
# Read the differential expression results from the "resultados_diff.txt" file
# The data is read with columns separated by spaces, and the decimal point as '.'
data_diff <- read.table("resultados_diff.txt", header = TRUE, sep = " ", dec = ".")
data_diff$gene_id <- rownames(data_diff)
colnames(data_diff) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_id")


########################################################################################################
################# APPLICATION OF THE HEATMAP  TO REPRESENT DIFFERENTIAL GENE EXPRESSION  ###############
########################################################################################################

# Load necessary libraries
library("RColorBrewer")  # RColorBrewer is used for color palettes
library("circlize")      # circlize provides circular visualization tools     
library(ComplexHeatmap)  # ComplexHeatmap is used for creating complex heatmaps  

# Filter results to include genes with an adjusted p-value < 0.001
# Convert the filtered results to a matrix for further analysis
best <- as.matrix(subset(resAB, padj < 0.001))

# Set the criterion for sorting (e.g., adjusted p-value)
criterio <- "padj"

# Sort genes based on the chosen criterion
genes_ordenados <- resAB[order(resAB[[criterio]]), ]

# Select the top 100 genes based on the sorted results
top_100_genes <- genes_ordenados[1:100, ]

# Display the top 100 genes for verification
head(top_100_genes)

# Extract normalized expression data from DESeq2 object
expressao_data_normal <- counts(ddsAB, normalized=FALSE)


# Subset the normalized expression data to include only the top 100 genes
expressao_data_normal <- subset(counts(ddsAB, normalized=FALSE), 
                                rownames(counts(ddsAB, normalized=FALSE)) %in% 
                                  rownames(top_100_genes))

# Scale the expression data for visualization
ztexp <- scale(t(expressao_data_normal))
ztexp <- t(ztexp)

# Check the range of scaled expression values
max(ztexp)
min(ztexp)

# Define a color ramp for the heatmap
mycol <- colorRamp2(c(-1, 1, 2), c("white", "yellow", "red"))

# Create a heatmap with the top 100 genes
Heatmap(ztexp, k = 2, col = mycol,
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 5))

## An example:

########################################################################################################
#################################### GENE ONTOLOGY TERM ENRICHMENT ANALYSIS#############################
########################################################################################################


#####################################################
################# Hypergeometric Test ############### 
#####################################################

### Example R Code for Hypergeometric Test

# Example values for the hypergeometric test
# Assume we are interested in the GO term "biosynthesis of anthocyanins"
# with ID GO:0046283
# Total number of genes in the genome
N <- 20000  

# Total number of genes annotated with the GO term "biosynthesis of anthocyanins"
M <- 500    

# Number of genes in the target set (e.g., differentially expressed genes)
n <- 100    

# Number of genes annotated with the GO term in the target set
k <- 20     

# Perform hypergeometric test
p_value <- phyper(k-1, M, N-M, n, lower.tail=FALSE)
p_value


### Gene Ontology (GO) enrichment analysis in blueberry

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("AnnotationForge",
#                       "clusterProfiler",
#                       "enrichplot"))


#####################################################
########  Performing the test as plot graphic ####### 
#####################################################

#library(GO.db)
#library(AnnotationDbi)
library(AnnotationForge)
library(clusterProfiler)
library(enrichplot)

#fSym <- read.table("func_annotation/gene_info_2.txt", header=TRUE, sep="\t", quote="")
#colnames(fSym) <- c("GID","SYMBOL","GENENAME")
#head(fSym)

#fGO <- read.table("func_annotation/go_annotation_2.tsv", header=TRUE, sep="\t", quote="")
#colnames(fGO) <- c("GID","GO","EVIDENCE")
#head(fGO)

#fChr <- read.table("func_annotation/chr_data.txt", header=TRUE, sep="\t", quote="")
#colnames(fChr) <- c("GID","CHROMOSOME")
#head(fChr)

#**WARNING!**: This step is essential for the construction of the OrgDB (annotation database) that will be used to run clusterProfiler. We will not execute it because this code has already been run along with the command "R CMD build org.Vcorymbosum.eg.db/" in the terminal. The resulting file is located at "ordb/org.Vcorymbosum.eg.db_0.0.3.tar.gz".

#makeOrgPackage(
#  gene_info=fSym, chromosome=fChr, go=fGO,
#  version = "0.0.3",
#  author = "Edson7",
#  maintainer = "Edson6 <mariodeandradee@gmail.com>",
#  outputDir = ".",
#  genus="Vaccinium",
#  species="corymbosum",
#  tax_id = "69266",
#  goTable="go"
#)

## then you can call install.packages based on the return value.
## run the command on terminal: R CMD build org.Vcorymbosum.eg.db/
#install.packages("ordb/org.Vcorymbosum.eg.db_0.0.3.tar.gz", repos=NULL)

library(org.Vcorymbosum.eg.db)


gene_list <- readLines("../vcorimb/func_annotation/gene_list_UP.txt")


gene_list

# Realizar a análise de enriquecimento GO com IDs GO válidos
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Vcorymbosum.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # "BP" para Biological Process, "CC" para Cellular Component, "MF" para Molecular Function
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05
)

mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

dotplot(ego, showCategory=15) + ggtitle("dotplot for ORA")


#Performing the test using the downregulates genes data.

gene_list <- readLines("../vcorimb/func_annotation/gene_list_DOWN.txt")

head(gene_list)

# Realizar a análise de enriquecimento GO com IDs GO válidos
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Vcorymbosum.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # "BP" para Biological Process, "CC" para Cellular Component, "MF" para Molecular Function
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05
)

print(ego)


mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

dotplot(ego, showCategory=15) + ggtitle("dotplot for ORA")

