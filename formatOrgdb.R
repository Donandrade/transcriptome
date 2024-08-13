library(GO.db)
library(AnnotationDbi)
library(AnnotationForge)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
# Suponha que você tenha os seguintes arquivos

fSym <- read.table("func_annotation/gene_info_2.txt", header=TRUE, sep="\t", quote="")

colnames(fSym) <- c("GID","SYMBOL","GENENAME")

fGO <- read.table("func_annotation/go_annotation_2.tsv", header=TRUE, sep="\t", quote="")
colnames(fGO) <- c("GID","GO","EVIDENCE")

fChr <- read.table("func_annotation/chr_data.txt", header=TRUE, sep="\t", quote="")
colnames(fChr) <- c("GID","CHROMOSOME")


makeOrgPackage(
  gene_info=fSym, chromosome=fChr, go=fGO,
  version = "0.0.3",
  author = "Edson7",
  maintainer = "Edson6 <mariodeandradee@gmail.com>",
  outputDir = ".",
  genus="Vaccinium",
  species="corymbosum",
  tax_id = "69266",
  goTable="go"
)


## then you can call install.packages based on the return value. run the command on terminal: R CMD build org.Vcorymbosum.eg.db/
install.packages("org.Vcorymbosum.eg.db_0.0.3.tar.gz", repos=NULL)

library(org.Vcorymbosum.eg.db)

gene_list <- readLines("../vcorimb/func_annotation/gene_list_UP.txt")

head(gene_list)

# Realizar a análise de enriquecimento GO com IDs GO válidos
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Vcorymbosum.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # "BP" para Biological Process, "CC" para Cellular Component, "MF" para Molecular Function
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)



mutate(ego, qscore = -log(p.adjust, base=20)) %>% 
  barplot(x="qscore")

dotplot(ego, showCategory=30) + ggtitle("dotplot for ORA")

#################DOWN


gene_list <- readLines("../vcorimb/func_annotation/gene_list_DOWN.txt")

head(gene_list)

# Realizar a análise de enriquecimento GO com IDs GO válidos
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Vcorymbosum.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # "BP" para Biological Process, "CC" para Cellular Component, "MF" para Molecular Function
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

print(ego)

mutate(ego, qscore = -log(p.adjust, base=20)) %>% 
  barplot(x="qscore")

dotplot(ego, showCategory=30) + ggtitle("dotplot for ORA")

