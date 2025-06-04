# RNAseq Demo Project
### Methods

Three replicates for RNAseq of a WT and KO cell line will be downloaded from public accession database.  Read QC will be performed with Fastqc (v0.11.7). STAR (v2.7.11b) will be used to align reads and Samtools (v1.19.2) flagstats will assess alignment quality.  

Reads will be counted using Verse (v0.1.5) and filtered for genes that have counts > 0 for at least two of three replicates.  DESeq2 (v3.19) will be used for differential expression analysis in R and fgsea (v3.19) will be used for Gene Set enrichement analysis while ernichR (v3.2) will be used for gene enrichment analysis.

### Subsetting Samples 

I needed to start by subsetting samples to test my script to make testing more efficient.
I opted to used seqtk to subsample 10000 reads from each of the .fq.gz files:

```bash
seqtk sample -s100 [sample_name]_R1.fastq.gz 10000 > sub_[sample_name]_R1.fastq  
seqtk sample -s100 [sample_name]_R2.fastq.gz 10000 > sub_[sample_name]_R2.fastq 
```

### Quality of Reads and Alignment

Quality of Sequencing Reads (fastqc) 

![](/images/read_qc.png)
The read quality looks normal in terms of GC content and total reads (80-120 M).  While the duplication levels would appear to be concerning, this is normal for RNAseq data since there tend to be some sequences that are highly expressed relative to most leading to high duplicates in fastQC.  Read quality (per sequence and mean) also are all high ~35 so I would say that the read quality is good. 



Quality of Alignment (flagstats)

![](/images/alignment_qc.png)
The alignment does not show glaring issues.  For all samples >90% of reads are uniquely mapped and there are very few unmapped reads.

Based on quality, I will not exclude any samples from the analysis.

### PCA and Clustering

![](/images/pca.png)
Clearly there is one KO sample that is different from the other two in terms of overall expression profile.  This also is clear in the sample distance plot.

![](/images/sample_distance.png)

### Setting FDR and Determining Signifiance Genes

DESeq by default calculated the Benjamini-Hochberg adjusted p-values for results so I will use this method.  The default alpha (false discovery rate) is 0.1.  The method ranks all p-values and then calculated the adjusted p-value as  (p-value)/(rank)*(number of tests).  I tested a few different rates and saw slight variation in number of significantly DE genes, but will choose to use the default with the knowledge that this assumes 10% of the significnatly DE may be false positives.

With an FDR threshold of 0.1 there are 686 significantly upregulated genes and 822 significantly downregulated genes.

### Histogram of DE Genes
![](/images/log2_fold_change_histogram.png)

### Volcano Plot of DE Genes
![](/images/differentially_expressed_genes_volcano.png)

### GSEA Results
![](/images/fgsea_sig.png)
![](/images/fgsea_all.png)

### enrichR Results

![](/images/enrich_UP.png)
![Alt text](/images/enrich_DOWN.png)

### Discussion of Results

I used the MSigDB Hallmark 2020 database for both the enrichR enrichment analysis and the GSEA.  I used log2(foldchange) as a ranking metric for GSEA.  Given that the enrichment analysis used the exact same genes as the GSEA with only significantly DE genes whereas the GSEA with all genes included more data, I am not surprised that the GSEA with only DE genes is more similar.

Generally the results for all three analyses are highly similar with some minor differences.  We can see that P53 Pathways, KRAS signaling, Inflammatory response, coagulation, and Xenobiotic metabolism are among the top enriched.  Among the least enriched that the three analyses share are Pancreas Beta Cells, UV Response, KRAS Signaling, and Estrogen response.  I think these results support the claim that the enrichment analysis, while perhaps not as sensitive as GSEA, still provides consistent results for which pathways are enriched based on RNAseq results.  I would still defer to the GSEA with all genes since it uses the full dataset, but it seems highly similar to the results of the other two analyses.  Also, since more genes are searched using the full gene set, there is more statistical power and more gene sets are significantly enriched.

The analysis indicates that KO cells relative to controls have upregulation of the P53 tumor suppression pathway and Epithelial Mesenchymal Transition which relates to wound healing and downregulation of pancreas beta cells and KRAS signaling, which is involved in cellular communication via signal transduction.  These results generally are concordant with the identity of the KO, which is a KO in TYK2, a kinase involved in intracellular signaling for insulin regulation.  Chandra et al. 2022 also saw significantly reduced enrichment of pancreas beta cell pathways, although they used a different gene set for their enrichment analysis.  Furthermore, they observed strong effects of the TYK2 KO on KRAS signaling, which is also borne out in my results.  