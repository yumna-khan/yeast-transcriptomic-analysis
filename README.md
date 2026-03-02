# yeast-transcriptomic-analysis
## General Overview


## Table of Contents
- [Introduction](#introduction)
- [Methods](#methods)
  - [1. Data Description](#1-data-description)
  - [2. Quality Control](#2-quality-control)
  - [3. Quantification](#3-quantification)
  - [4. Differential Expression](#4-differential-expression)
  - [5. Visualization & Functional Annotation](#5-visualization-&-functional-annotation)
- [Results](#results)
- [Discussion](#discussion)
- [References](#references)

## Introduction
Flor yeast _Saccharomyces cerevisiae_ is widely used in fermentation industries, particularly in the production of sherry-type wines. During biological ageing, specific flor strains transition from planktonic (free-living, suspended) growth to form a floating biofilm known as velum on the surface of fermented liquid. Velum develops under harsh conditions, including high acetaldehyde concentrations, oxidative stress from metabolizing non-fermentable carbon sources, ethanol exposure, and water stress. These environmental pressures require adaptive responses that alter cell morphology, surface hydrophobicity, and metabolic activity [^1]. Prolonged stress may also lead to genomic variation, such as single nucleotide polymorphisms (SNPs), insertions and deletions (InDels), and chromosomal rearrangements, which can influence regulatory pathways and ultimately affect gene expression patterns [^2]. Characterizing transcriptional changes (transcriptional reprogramming) across different stages of velum development is therefore essential for understanding how flor yeast adapts to winemaking environments. This knowledge may help clarify the functions of previously uncharacterized genes, improve sherry wine production, and identify targets for strain optimization [^2].

Thus, comparative transcriptome analysis using bulk RNA sequencing (RNA-seq) provides a genome-wide method for quantifying gene expression changes across biological conditions. Initial data quality assessment can be performed using tools such as FastQC and MultiQC to detect sequencing errors, GC bias, or contamination [^3][^4]. Because excessive trimming may introduce bias in pseudo-alignment methods, trimming is generally recommended only when necessary [^5].

For transcript quantification, pseudo-alignment tools such as Salmon are commonly used due to their speed and computational efficiency. Salmon is used for quantification because its quasi-mapping approach skips mismatched bases, removing the need for trimming. It is faster and more accurate than comparable tools such as kallisto and eXpress, and handles both alignment and quantification in a single step. Salmon also corrects for technical errors introduced during library preparation, leading to more reliable estimates of gene activity in yeast metabolic studies [^6]. However, it is important to note that Salmon can show decreased accuracy for very short or lowly expressed transcripts, where limited sequence data makes reliable quantification more difficult [^8].

For statistical comparison across biofilm stages, DESeq2 is widely applied due to its ability to model count data using a negative binomial distribution and stabilize dispersion estimates through empirical Bayes shrinkage. It performs well with small sample sizes and controls the false discovery rate (FDR), making it suitable for multi-group transcriptomic comparisons [^8]. Although it does not easily accommodate complex mixed-effect models, it remains a robust and validated choice for standard differential expression analysis.

To interpret differentially expressed genes, functional enrichment tools such as clusterProfiler can be applied. Given the well-annotated genome of _S. cerevisiae_, approaches such as Over-Representation Analysis (ORA) allow for the identification of enriched biological processes involved in ethanol metabolism, oxidative stress response, and cell wall remodeling. ClusterProfiler improves interpretability by reducing redundancy among enriched Gene Ontology terms [^9]. However, its results are dependent on the quality and completeness of existing annotations, and very small or sparse gene lists may reduce statistical power and biological interpretability.

Overall, the objective of this study is to perform differential expression analysis and functional annotation across three stages of velum development, each with three biological replicates, in order to characterize the transcriptional programs underlying yeast biofilm formation during wine ageing.



## Methods
### 1. Data Description
The dataset was obtained from the NCBI BioProject PRJNA592304 and includes three biological replicates for each developmental stage: early biofilm, thin biofilm, and mature biofilm. Transcript quantification with Salmon required a reference transcriptome (cDNA), which was retrieved from Ensembl along with the corresponding GTF annotation file. Initial quantification performed using only protein-coding cDNA sequences resulted in a suboptimal mapping rate of approximately 41%, whereas expected mapping rates are typically around 80–90% [^12]. In yeast, the transcriptome annotation includes both coding (cDNA) and non-coding RNA (ncRNA) sequences [^13]. Therefore, to improve mapping sensitivity and overall alignment performance, the reference was expanded to include ncRNA sequences in addition to cDNA.

### 2. Quality Control
Quality control of the raw sequencing reads was performed in Ubuntu using FastQC (v0.12.1), which generates HTML reports summarizing key metrics such as per-base sequence quality, GC content, adapter contamination, and sequence duplication levels. FastQC is widely used in RNA-seq workflows to identify potential technical issues prior to downstream analyses [^3]. 

```
fastqc *.fastq.gz
```

To facilitate simultaneous inspection of all samples, MultiQC (v1.25.1) was used to aggregate individual FastQC reports into a single consolidated summary file, enabling efficient cross-sample comparison.

```
multiqc .
```

### 3. Quantification
Salmon (v1.10.3) requires two main steps: indexing the reference transcriptome and quantifying each sample. 

First, the _Saccharomyces cerevisiae_ cDNA and ncRNA reference files were combined to ensure quantification against the complete transcriptome, as yeast transcriptomes include both coding and non-coding RNA [^13].
```
zcat Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz > yeast_transcriptome_combined.fa
```

Next, the combined transcriptome was indexed. The `-t` parameter specifies the transcript FASTA file, and `-i` defines the name of the index directory. While Salmon’s default k-mer size is 31 (optimized for reads ≥75 bp), this was reduced to 21 because the reads in this dataset are 50 bp long [^7]. Using a smaller k-mer improves sensitivity for shorter reads.
```
salmon index -t yeast_transcriptome_combined.fa -i combined_index -k 21
```

Before quantifying all samples, a preliminary test run was performed on a single file to confirm successful indexing and assess the mapping rate. The `-l A` parameter allows Salmon to automatically detect library type (stranded or unstranded), `-r` specifies single-end reads, and `--validateMappings` enables selective alignment, which improves accuracy.
```
salmon quant -i combined_index \
             -l A \
             -r SRR10551657.fastq.gz \
             --validateMappings \
             -o test_boost_quant

# Check mapping rate
grep "Mapping rate" test_boost_quant/logs/salmon_quant.log
```
Once the preliminary test was complete, all nine samples were quantified through a loop script. 

Finally, before proceeding to differential expression analysis, a sanity check was performed to confirm successful quantification. For yeast RNA-seq datasets, mapping rates are typically expected to fall between 80–90% [^12].
```
grep "Mapping rate" quants/*_quant/logs/salmon_quant.log
```
 

### 4. Differential Expression


### 5. Visualization & Functional Annotation



## Results



**Figure 1:** 



**Figure 2:** 



**Figure 3:** 



**Figure 4:** 


## Discussion
### 1. 


### 6. Conclusion


## References
[^1]: Bakker, J. (2013). Flor Yeasts. ScienceDirect. https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/flor-yeasts
[^2]: Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome profile of yeast strain used for biological wine aging revealed dynamic changes of gene expression in course of flor development. Frontiers in Microbiology, 11, 538. https://doi.org/10.3389/fmicb.2020.00538
[^3]: Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data. (n.d.). https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[^4]: Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
[^5]: Williams, C. R., Baccarella, A., Parrish, J. Z., & Kim, C. C. (2016b). Trimming of sequence reads alters RNA-Seq gene expression estimates. BMC Bioinformatics, 17(1), 103. https://doi.org/10.1186/s12859-016-0956-2
[^6]: Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197
[^7]: Wu, D. C., Yao, J., Ho, K. S., Lambowitz, A. M., & Wilke, C. O. (2018). Limitations of alignment-free tools in total RNA-seq quantification. BMC Genomics, 19(1), 510. https://doi.org/10.1186/s12864-018-4869-5
[^8]: Wang, T., Li, B., Nelson, C. E., & Nabavi, S. (2019). Comparative analysis of differential gene expression analysis tools for single-cell RNA sequencing data. BMC Bioinformatics, 20(1), 40. https://doi.org/10.1186/s12859-019-2599-6
[^9]: De Oliveira, F. H. S., Gomes, F. A., & Feltes, B. C. (2026). Benchmarking multiple gene ontology enrichment tools reveals high biological significance, ranking, and stringency heterogeneity among datasets. Frontiers in Bioinformatics, 6. https://doi.org/10.3389/fbinf.2026.1755664
[^10]: Mistry, M. P. a. M. (2017, August 21). Quantitation of transcript abundance using Salmon. Introduction to RNA-seq using high performance computing (Orchestra) - ARCHIVED. https://hbctraining.github.io/Intro-to-rnaseq-hpc-orchestra/lessons/09_salmon.html
[^11]: I Love, M., Anders, S., & Huber, W. (2025, December 1). Analyzing RNA-seq data with DESEQ2. Bioconductor. https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
[^12]: Lu, Z., & Lin, Z. (2019). Pervasive and dynamic transcription initiation in Saccharomyces cerevisiae. Genome Research, 29(7), 1198–1210. https://doi.org/10.1101/gr.245456.118
[^13]: Jiang, Z., Zhou, X., Li, R., Michal, J. J., Zhang, S., Dodson, M. V., Zhang, Z., & Harland, R. M. (2015). Whole transcriptome analysis with sequencing: methods, challenges and potential solutions. Cellular and Molecular Life Sciences, 72(18), 3425–3439. https://doi.org/10.1007/s00018-015-1934-y
[^14]: Lancelle, L. J., Potru, P. S., Spittau, B., & Wiemann, S. (2025). DgeaHeatmap: an R package for transcriptomic analysis and heatmap generation. Bioinformatics Advances. https://pmc.ncbi.nlm.nih.gov/articles/PMC12401572/
[^15]: Lynch, D. B., Logue, M. E., Butler, G., & Wolfe, K. H. (2010). Chromosomal G + C Content Evolution in yeasts: Systematic interspecies differences, and GC-Poor troughs at Centromeres. Genome Biology and Evolution, 2, 572–583. https://doi.org/10.1093/gbe/evq042
[^16]: Alexandre, H. (2013, October 15). Flor yeasts of Saccharomyces cerevisiae—Their ecology, genetics and metabolism. ScienceDirect. https://www.sciencedirect.com/science/article/abs/pii/S0168160513004078
[^17]: Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020b). Transcriptome profile of yeast strain used for biological wine aging revealed dynamic changes of gene expression in course of flor development. Frontiers in Microbiology, 11, 538. https://doi.org/10.3389/fmicb.2020.00538
[^18]: David-Vaizant, V., & Alexandre, H. (2018). Flor yeast diversity and dynamics in biologically aged wines. Frontiers in Microbiology, 9, 2235. https://doi.org/10.3389/fmicb.2018.02235
[^19]: Zara, S., Gross, M. K., Zara, G., Budroni, M., & Bakalinsky, A. T. (2010). Ethanol-Independent Biofilm Formation by a Flor Wine Yeast Strain ofSaccharomyces cerevisiae. Applied and Environmental Microbiology, 76(12), 4089–4091. https://doi.org/10.1128/aem.00111-10
[^20]: Váchová, L., Čáp, M., & Palková, Z. (2012). Yeast Colonies: A model for studies of aging, environmental adaptation, and Longevity. Oxidative Medicine and Cellular Longevity, 2012, 1–8. https://doi.org/10.1155/2012/601836
[^21]: Caron-Godon, C. A., Collington, E., Wolf, J. L., Coletta, G., & Glerum, D. M. (2024). More than Just Bread and Wine: Using Yeast to Understand Inherited Cytochrome Oxidase Deficiencies in Humans. International Journal of Molecular Sciences, 25(7), 3814. https://doi.org/10.3390/ijms25073814
[^22]: Baumann, L., Doughty, T., Siewers, V., Nielsen, J., Boles, E., & Oreb, M. (2021). Transcriptomic response of Saccharomyces cerevisiae to octanoic acid production. FEMS Yeast Research, 21(2). https://doi.org/10.1093/femsyr/foab011
[^23]: Pisithkul, T., Schroeder, J. W., Trujillo, E. A., Yeesin, P., Stevenson, D. M., Chaiamarit, T., Coon, J. J., Wang, J. D., & Amador-Noguez, D. (2019). Metabolic Remodeling during Biofilm Development of Bacillus subtilis. mBio, 10(3). https://doi.org/10.1128/mbio.00623-19
[^24]: Yi, H., Ahn, Y., Song, G. C., Ghim, S., Lee, S., Lee, G., & Ryu, C. (2016). Impact of a Bacterial Volatile 2,3-Butanediol on Bacillus subtilis Rhizosphere Robustness. Frontiers in Microbiology, 7, 993. https://doi.org/10.3389/fmicb.2016.00993
[^25]: GonzáLez, E., FernáNdez, M. R., Marco, D., Calam, E., Sumoy, L., ParéS, X., Dequin, S., & Biosca, J. A. (2009). Role of Saccharomyces cerevisiae Oxidoreductases Bdh1p and Ara1p in the Metabolism of Acetoin and 2,3-Butanediol. Applied and Environmental Microbiology, 76(3), 670–679. https://doi.org/10.1128/aem.01521-09
