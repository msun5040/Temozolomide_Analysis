library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(gprofiler2)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pathview)
library(cluster)

# ======================== DATA LOADING ===============================
#cell.count <- read.table('/Users/michaelsun/LawlerLab/H3G34R_data/G34R_counts_t.txt', header=TRUE, row.names = 1)
cell.count <- read.csv('/Users/michaelsun/LawlerLab/Temozolomide_study/Ntafoulis_string_filtered.csv', header=TRUE, row.names=1)
metadata <- read.csv("/Users/michaelsun/LawlerLab/Temozolomide_study/Ntafoulis_metadata.csv", header=TRUE, row.names=1)

# ======================= PREPROCESSING ===============================
# remove rows with all zeros
# cleaned_counts <- cell.count[apply(cell.count, 1, var) != 0,]
threshold <- 50
cleaned_counts <- cell.count[rowMeans(cell.count) > threshold, ]

# cleaning out trash genes
patterns_to_remove <- c("^MT", "^RP","^H[0-9]")
combined_pattern <- paste(patterns_to_remove, collapse="|")
cleaned_counts <- cleaned_counts[!grepl(combined_pattern, rownames(cleaned_counts)),]

dds <- DESeq(dds)
norm_counts <- counts(dds, normalized=TRUE)
pca_res <- prcomp(t(log2(norm_counts + 1)))

# data frame for plotting PCA
pca_data <- data.frame(Sample = colnames(cleaned_counts), PC1 = pca_res$x[,1], PC2 = pca_res$x[,2])
p <- ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) + 
  geom_point() + 
  geom_text(nudge_x = 0.02, nudge_y = 0.02, check_overlap = TRUE) + 
  theme_minimal()
ggsave("preSeq_PCA.png", plot = p, width = 10, height = 8, dpi = 300)

# ======================= DESeq RUN ===================================
# first round of deseq
dds <- DESeqDataSetFromMatrix(countData = cleaned_counts,
                              colData = metadata,
                              design = ~1)
dds <- DESeq(dds)

# ======================= RESULTS =====================================
norm_counts <- counts(dds, normalized=TRUE)
keep <- rowSums(norm_counts >= 10) >= round(ncol(cleaned_counts) * 0.1)
dds <- dds[keep,]

# variance stabilizing transform
vsd <- vst(dds, blind = FALSE)

# plot PCA - regular
png(filename = "regular_pca.png")
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
# regular plotting
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  geom_text_repel(aes(label = name), max.overlaps=Inf)
dev.off()

# PCA with clustering
png(filename = "clustering_pca.png")
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)

# regular plotting
k <- 2
clusters <- kmeans(pca_data[, c("PC1", "PC2")], centers = k)
pca_data$cluster <- as.factor(clusters$cluster)

ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster, label = name)) +
  geom_point() +
  geom_text_repel(aes(label = name), max.overlaps=Inf)
dev.off()
clustered_samples <- split(pca_data$name, pca_data$cluster)

# silhouette analysis for clustering
silhouette_scores <- silhouette(clusters$cluster, dist(t(log2(norm_counts + 1))))

png(filename = "silhouette_plot.png", width = 800, height = 600)
plot(silhouette_scores)
dev.off()

avg_silhouette_width <- mean(silhouette_scores[, 3])
cat("Average silhouette width:", avg_silhouette_width, "\n")

clustering_metadata <- read.csv("/Users/michaelsun/LawlerLab/Temozolomide_study/Ntafoulis_cluster_roster.csv")
clustering_metadata$cluster <- NA

# creating metadata for clustering deseq
for(group_name in names(clustered_samples)) {
  samples_in_group <- clustered_samples[[group_name]]
  
  for(sample in samples_in_group) {
    row_idx <- which(clustering_metadata$sample_name == sample)
    clustering_metadata$cluster[row_idx] <- paste0("Cluster", group_name)
  }
}

write.csv(clustering_metadata, "/Users/michaelsun/LawlerLab/Temozolomide_study/Ntafoulis_cluster_roster.csv", row.names = FALSE)

# second round of deseq
dds_clustering <- DESeqDataSetFromMatrix(countData = cleaned_counts,
                                         colData = clustering_metadata,
                                         design = ~cluster)
dds_clustering <- DESeq(dds_clustering)

# Get results table
res12 <- results(dds_clustering, contrast=c("cluster","Cluster1","Cluster2"))
res21 <- results(dds_clustering, contrast=c("cluster","Cluster2","Cluster1"))

# ============= WRITING =============================
sig_genes_12 <- subset(res12, padj < 0.05)
sig_genes_21 <- subset(res21, padj < 0.05)


sig_genes_12 <- sig_genes_12[order(sig_genes_12$log2FoldChange, decreasing = TRUE), ]
sig_genes_21 <- sig_genes_21[order(sig_genes_21$log2FoldChange, decreasing = TRUE), ]

write.csv(sig_genes_12, "/Users/michaelsun/LawlerLab/Temozolomide_study/sig_genes_12.csv")
write.csv(sig_genes_21, "/Users/michaelsun/LawlerLab/Temozolomide_study/sig_genes_21.csv")

upreg_12 <- subset(sig_genes_12, log2FoldChange>2)
upreg_21 <- subset(sig_genes_21, log2FoldChange>1.5)

# write results to csv
write.csv(upreg_12, "/Users/michaelsun/LawlerLab/Temozolomide_study/cluster1_upreg.csv")
write.csv(upreg_21, "/Users/michaelsun/LawlerLab/Temozolomide_study/cluster2_upreg.csv")

# ============= VISUALIZING =============================
gene_list <- rownames(upreg_12)
#gene_list <- rownames(upreg_21)
entrez_ids <- mapIds(org.Hs.eg.db, keys=gene_list, column="ENTREZID", keytype="SYMBOL", multiVals="first")

enrich_res <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)

barplot(enrich_res)
dotplot(enrich_res)

pathway_id <- as.character(enrich_res$ID[1])  # Get the ID of the first significant pathway
pathview(gene.data = entrez_ids, pathway.id = pathway_id, species = "hsa")

# MA Plot for LFC visualization
ma_data <- as.data.frame(res12)
ma_data$A <- (log2(ma_data$baseMean) + log2(ma_data$baseMean + 1)) / 2 
ma_data$M <- ma_data$log2FoldChange

g <- ggplot(ma_data, aes(x=A, y=M)) +
  geom_point(alpha=0.4, color="blue") +
  theme_minimal() +
  labs(x="Average Log Expression (A)", y="Log2 Fold Change (M)", title="MA Plot") +
  geom_hline(yintercept=0, linetype="dashed", color="red")

# Save the plot
ggsave("MA_plot.png", plot = g, width = 10, height = 8, dpi = 300)

