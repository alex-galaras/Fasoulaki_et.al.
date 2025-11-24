# README.md
This repository contains the scripts for performing RNA-seq & ATAC-seq analysis for the paper Fasoulaki *et.al.* 2025, currently found in the link: 
https://www.biorxiv.org/content/10.1101/2025.11.21.689716v1.full.pdf

# **RNA-seq analysis**

## Alignment

## **MetaseqR2 Analysis Pipeline for the Delidakis Quant-Seq Dataset**


The analysis is performed using **metaseqR2**, integrating multiple statistical methods and generating QC reports, differential expression results, and UCSC Genome Browser tracks.

## **1. Overview**

This analysis uses **Quant-Seq 3’ mRNA sequencing** data from Drosophila melanogaster (*dm6*) to evaluate transcriptome‐wide expression differences across four biological contrasts, including:

* trxRNAi-Suc vs cnt-Suc
* trxRNAi-DSS vs cnt-DSS
* cnt-Suc vs cnt-DSS
* trxRNAi-Suc vs trxRNAi-DSS

The workflow performs:

* Read count normalization
* Differential expression using multiple statistical engines
* Meta-analysis using Pandora
* Comprehensive QC plotting
* Export of statistics, normalized counts, fold-changes, and annotations
* Generation of UCSC browser track hubs

---

## **2. Requirements**

### **Software**

* **R ≥ 4.0**
* **metaseqR2**
* The following R packages (loaded automatically by metaseqR2):

  * DESeq
  * DESeq2
  * edgeR
  * NBPSeq
  * NOISeq
  * others (QC and plotting dependencies)

### **External Files**

* **targets_delidakis.txt**
  A tab-delimited file describing samples, groups, and count files.

* **annotation.sqlite**
  SQLite annotation DB for dm6 (local reference database).

---

## **3. Usage**

### **Run the analysis**

From the project directory:

```bash
Rscript run_metaseqR2_delidakis.R
```

This will:

* Perform normalization and differential expression
* Integrate statistics using Pandora
* Generate QC plots
* Export all results to the specified output directory
* Produce UCSC Genome Browser tracks
* Create summary reports

### **Input Target File**

The script expects the target file:

```
targets_delidakis.txt
```
