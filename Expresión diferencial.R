HISAT2 # instalar 
hisat2 -p 4 -x /home/alumno26/dc_workshop/Genome/genome \-1 /home/alumno26/dc_workshop/data/trimmed_fastq/JC1A_R1.trim.fastq.gz \-2 /home/alumno26/dc_workshop/data/trimmed_fastq/JC1A_R2.trim.fastq.gz \-S /home/alumno26/dc_workshop/mapping_results/JC1A.sam
featureCounts -p -t exon -g gene_id \
-a /home/alumno5/dc_workshop/Genome/Arabidopsis_thaliana.cop.TAIR10.57.gtf \
-o JC1A.txt /home/alumno5/dc_workshop/mapping_results/JC1A_sorted.bam #Cuantificar a nivel de genes el nivel de expresión

/home/alumno26/dc_workshop

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
install.packages(c("ggplot2", "pheatmap"))

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

setwd("/home/alumno26/Curso_RNAseq/")
countData <- read.table("matriz_arabidopsis_2023.txt", header=TRUE, row.names=1, sep="\t")
column_order <- c("control_1", "control_2", "control_3", "treatment_1", "treatment_2","treatment_3")
countData <- countData[, column_order]

condition <- factor(rep(c("Control", "Treatment"), each=3))
colData <- data.frame(row.names = colnames(countData), condition)

dim(countData)

plot(log2(countData$control_1), log2(countData$treatment_1))
plot(log2(countData$control_1), log2(countData$treatment_2))

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
dds <- dds[rowSums(counts(dds)) > 10,]
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
summary(res)
head(res)  

rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup = "condition") +
  geom_text(aes(label = name), vjust = 0.5) +
  theme_minimal()

EnhancedVolcano(res, lab = rownames(res), x = "log2FoldChange", y = "pvalue")
plotMA(res, alpha = 0.05, main = "Diferencias de Expresión: Tratamiento vs Control")

names(res)
plotMA(res, alpha = 0.05, main = "Diferencias de Expresión: Tratamiento vs Control",ylim=c(-10,10))

top_genes=row.names(res)(1-20)


library(limma)
library(edgeR)

setwd("/home/alumno26/Curso_RNAseq")
list.files()

outpath = "/home/alumno26/Curso_RNAseq/"
dir.create(outpath, showWarnings=FALSE)

counts = read.table("matriz_arabidopsis_2023.txt", header=TRUE, row.names = 1, sep="\t", comment.char="") 
head(counts)
dim(counts)
counts = counts[rowSums(cpm(counts) >= 1) >=3,]
head(counts)
dim(counts)

colnames(counts)

plot(log2(counts[,c("control", "tratamiento")]), col="gray")

grp = sub("..$", "", colnames(counts)) 
grp

dge = DGEList(counts=counts, group=grp)
dge
dgeNorm = calcNormFactors(dge)
dgeNorm

plotMDS(dge)

dgeNorm = estimateCommonDisp(dgeNorm)
dgeNorm = calcNormFactors(dge)
dgeNorm$common.dispersion

diff_exp = exactTest(dgeNorm, dispersion = dgeNorm$common.dispersion, pair = c("control", "treatment" ))
diff_exp

dim(diff_exp)
topTags(diff_exp)
dim(topTags(diff_exp))
deTab = topTags(diff_exp, n=Inf)$table
deTab[c(15,30),]
row.names(deTab)[deTab$logFC > 5] 

deTab["YNL284C-A",]
deGenes = rownames(deTab)[deTab$FDR < 0.05 & abs(deTab$logFC) > 1]
down=row.names(deTab)[deTab$logFC< -1]
over=row.names(deTab)[deTab$logFC> 1]
print(paste("total de diferenciales:", length(deGenes)))
print(paste("número de genes inducidos:", length(over)))
print(paste("número de genes reprimidos:", length(down)))

write.table(deTab, file=paste(outpath, "diff_gene_wt.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(down, file=paste(outpath, "down.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(over, file=paste(outpath, "up.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")




