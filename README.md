# ATACseq Project Proposal
## Methods

Two replicates for atacseq will be downloaded from public accession database.  Read QC will be performed with Fastqc (v0.11.7) and adapter trimming will be performed with trimmomatic (0.39).  BowTie2 (v2.5.1) (args: -X 2000) will be used for alignment.  Mitochondrial reads will be removed using samtools (v1.12) view and alignment QC will be performed with samtools flagstats.  Tagmentation bias of transposase preferentially cleaving certain motifs will be normalized using deeptools (v3.5.1) (args: alignmentSieve --ATACshift).  Then ATACSeqQC (v4.3 via Bioconductor v3.18) will be used to filter out nuceosome bound reads and duplicate reads.

Peak calling will be performed using MACS3 (v3.0.1) and reproducible peaks will be defined as those existing in both replicates based on analysis using bedtools intersect (v2.31.0).  Analysis of peaks will include searching for motif enrichment, peak annotation, dsitribution of peaks in genomic regions, and potentially integration of peak data with RNAseq data.

## Questions

1. Are there any concerning aspects of the quality control of the sequencing reads?
2. Are there any concerning aspects of the quality control related to alignment?
3. Based on quality control, will any samples be excluded from further analysis?
4. How many peaks are present in each of the replicates?
5. How many peaks are present in the set of reproducible peaks? What strategy was used to determine “reproducible” peaks?
6. What can chromatin accessibility let us infer biologically?

## Deliverables

1. Produce a fragment length distribution plot for each of the samples
2. Produce a table of how many alignments for each sample before and after filtering alignments falling on the mitochondrial chromosome
3. Create a signal coverage plot centered on the TSS (plotProfile) for the nucleosome-free regions (NFR) and the nucleosome-bound regions (NBR)
4. A table containing the number of peaks called in each replicate, and the number of reproducible peaks
5. A single BED file containing the reproducible peaks you determined from the experiment.
6. Perform motif finding on reproducible peaks
7. Create a single table / figure with the most interesting results
8. Perform a gene enrichment analysis on the annotated peaks using a well-validated gene enrichment tool
9. Create a single table / figure with the most interesting results
10. Produce a figure that displays the proportions of regions that appear to have accessible chromatin called as a peak (Promoter, Intergenic, Intron, Exon, TTS, etc.)
