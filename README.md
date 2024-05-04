# ATACseq Project Proposal
## Methods

Three replicates for RNAseq of a WT and KO cell line will be downloaded from public accession database.  Read QC will be performed with Fastqc (v0.11.7). STAR (v2.7.11b) will be used to align reads and Samtools (v1.19.2) flagstats will assess alignment quality.  

Reads will be counted using Verse (v0.1.5) and filtered for genes that have counts > 0 for at least two of three replicates.  DESeq2 (v3.19) will be used for differential expression analysis in R and fgsea (v3.19) will be used for Gene Set enrichement analysis while ernichR (v3.2) will be used for gene enrichment analysis.

## Questions

1. Are there any concerning aspects of the quality control of your sequencing reads?
2. Are there any concerning aspects of the quality control related to alignment?
3. Based on quality control, will any samples be exluded from further analysis?
4. How many genes are significant at your chosen statistical threshold?
5. How similar are the results from GSEA and gene set enrichment? Are there any notable differences?
6. Are differences expected? If so, why?
7. What do the results imply about potentional biological functions of the factor of interest?

## Deliverables

1. Produce either a sample-to-sample distance plot from the counts matrix or a PCA biplot
2. A CSV containing all of the results from your DE analysis
3. A histogram showing the distribution of log2FoldChanges of your DE genes
4. A volcano plot that distinguishes between significant and non-significant DE genes as well as labels up- and downregulated genes
5. Label the top ten most significant genes with their associated gene name / gene symbol.
6. Perform a GSEA (FGSEA) analysis on all the genes discovered in the experiment
7. Create a single table / figure that reports the most interesting results
8. Use your list of DE genes at your chosen statistical threshold and perform a basic enrichment analysis using a tool of your choice
9. Create a single table / figure that reports the most interesting results
