### DESeq2 Analysis for Final Project ###

# Load Libraries
library("DESeq2")
library(tidyr)
library(dplyr)
library(fgsea)
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library(enrichR)
library(ggrepel)
### DESeq2 Analysis ###

### DID NOT PERFORM PREFILTERING (aside from removing rows where one sample had count = 0) ###

# Read in filtered Counts Data
data <- read.csv("/projectnb/bf528/students/nmf35/bf528-individual-project-NathanielFisher1/results/verse_concat_filtered.csv")

# make coldata dataframe based on conditions
coldata <- data.frame(Condition = strtrim(colnames(data)[-1], 2))
rownames(coldata) <- colnames(data)[-1]
coldata$Condition <- factor(coldata$Condition)

# format counts matrix
cts <- as.matrix(data[,-1])
rownames(cts) <- data[,1]

# make DESeq object with cts and coldata with design formula on Condition 
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Condition)


# Now add gene names and column data to dds
names <- read.table("/projectnb/bf528/students/nmf35/bf528-individual-project-NathanielFisher1/results/id2gene.txt")
colnames(names) <- c("id", "gene")
names_ordered <- tibble(id = rownames(cts))
add_gene_names <- tibble(names) 
final <- inner_join(names_ordered, add_gene_names, by = "id")
mcols(dds) <- DataFrame(final)


#calculate differential expression and look at results
dds <- DESeq(dds)
res <- results(dds)
final_results <- data.frame(res, gene = mcols(dds)$gene)
head(final_results)
dim(final_results)
summary(res)
plotMA(res, ylim=c(-5,5))

# first normalize data to remove dependent mean-variance relationship 
# transform data with variance stabilizing tranformation
vsd <- vst(dds, blind=FALSE)

# transform data with regularized logarithm
rld <- rlog(dds, blind=FALSE)

# sample to sample distance plot
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
# Result (one KO appears more similar to the CTL than the other two KOs)

# PCA
colData(dds)
plotPCA(vsd, intgroup=c("Condition"))


# one of the KOs is clearly different

# setting FDR (decided to use standard 0.01)

plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)


# get only significant results
final_results_ordered <- final_results[order(final_results$pvalue),]
final_results_Sig <- subset(final_results_ordered, padj < 0.1)

# write all results to csv
write.csv(final_results_ordered, 
          file="condition_treated_results.csv")


# make histogram of sig results
ggplot(final_results_Sig, aes(x = log2FoldChange)) + 
  geom_histogram(binwidth = .05, color = "black", fill = "red") +
  labs(title = "Log2FoldChange of Significantly DE Genes", y = "# Genes") 
  
# make a volcano plot of all results
volc_tibble <- tibble(final_results_ordered)

volc_df <- volc_tibble %>% mutate(DE = case_when(
  padj < 0.1 & log2FoldChange >= 1 ~ "UP",
  padj < 0.1 & log2FoldChange <= -1 ~ "DOWN",
  TRUE ~ "NS"
), log10p = -log10(padj))

color_scheme <- c("DOWN" = "red", "NS" = "lightgrey", "UP" = "green")


ggplot(volc_df, aes(x = log2FoldChange, y = log10p, color = DE)) + geom_point(size = .7) + scale_y_continuous(limits = c(-10,112))+
  scale_x_continuous(limits = c(-9,5))+
  labs(title = "Volcano Plot Showing Significantly Up and Down regulated Genes", y = "-log10(p-value)") +
geom_label_repel(data = subset(volc_df, log2FoldChange < -3.135 | log2FoldChange > 3.135),aes(label = gene),
                     show.legend = FALSE, 
                 size = 3,
                 box.padding = 1,
                     segment.color = "grey50",
                 max.overlaps = Inf,
                 colour = 'black'
                     )+
  scale_color_manual(values = color_scheme)




### GSEA ###
# Starting with sign gene set
# Reformat data and look at most and least DE genes
ordered_results <- final_results_Sig[order(-final_results_Sig$log2FoldChange),]
head(ordered_results)
rnk_list <- setNames(ordered_results$log2FoldChange, ordered_results$gene)
head(rnk_list)
tail(rnk_list)
length(rnk_list)
hallmark_pathways_fgsea <- fgsea::gmtPathways('/projectnb/bf528/students/nmf35/bf528-individual-project-NathanielFisher1/h.all.v7.5.1.symbols.gmt')
fgsea_results <- fgsea(hallmark_pathways_fgsea, rnk_list, minSize = 15, maxSize=500)


fgsea_results %>%
  mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot() +
  geom_bar(aes(x=pathway, y=NES, fill = padj < .25), stat='identity') +
  scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
  theme_minimal() +
  ggtitle('fgsea results for Hallmark MSigDB gene sets 
          for signidicantly DE Genes') +
  ylab('Normalized Enrichment Score (NES)') +
  xlab('') +
  coord_flip()

# now doing it for all genes
ordered_results <- final_results_ordered[order(-final_results_ordered$log2FoldChange),]
head(ordered_results)
rnk_list <- setNames(ordered_results$log2FoldChange, ordered_results$gene)
head(rnk_list)
tail(rnk_list)
length(rnk_list)
fgsea_results <- fgsea(hallmark_pathways_fgsea, rnk_list, minSize = 15, maxSize=500)


fgsea_results %>%
  mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
  ggplot() +
  geom_bar(aes(x=pathway, y=NES, fill = padj < .25), stat='identity') +
  scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
  theme_minimal() +
  ggtitle('fgsea results for Hallmark MSigDB gene sets 
  for all DE Genes') +
  ylab('Normalized Enrichment Score (NES)') +
  xlab('') +
  coord_flip()


# Gene functional enrichment

#used these dbs as a comprehensive RNAseq resource to start with some exploration, but then just used MSidDB HallMark
dbs <- c("RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",	 
"RNAseq_Automatic_GEO_Signatures_Human_Down", 
"RNAseq_Automatic_GEO_Signatures_Human_Up")

enriched <- enrichr(final_results_Sig$gene, "MSigDB_Hallmark_2020")
bp <- enriched[["MSigDB_Hallmark_2020"]]
#make plot
plotEnrich(enriched[[1]], showTerms = 40, numChar = 60, y = "Count", orderBy = "P.value", title = "Enrichment Analysis by EnrichR ")

# try separating them out for down and upregulated genes
UP <- final_results_Sig[final_results_Sig['log2FoldChange'] > 0,]
DOWN <- final_results_Sig[final_results_Sig['log2FoldChange'] < 0,]

enriched <- enrichr(UP$gene, "MSigDB_Hallmark_2020")
plotEnrich(enriched[[1]], showTerms = 40, numChar = 60, y = "Count", orderBy = "P.value", title = "Enrichment Analysis by EnrichR (Upregulated)")

enriched <- enrichr(DOWN$gene, "MSigDB_Hallmark_2020")
plotEnrich(enriched[[1]], showTerms = 40, numChar = 60, y = "Count", orderBy = "P.value", title = "Enrichment Analysis by EnrichR (Downregulated)")
