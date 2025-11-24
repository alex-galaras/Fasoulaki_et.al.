# README.md
This repository contains the scripts for performing RNA-seq & ATAC-seq analysis for the paper Fasoulaki *et.al.* 2025, currently found in the link: 
https://www.biorxiv.org/content/10.1101/2025.11.21.689716v1.full.pdf

# **RNA-seq analysis**

## **Read Alignment Pipeline (HISAT2 + Bowtie2) for dm6**

This repository contains a Bash script used to perform read alignment against the *Drosophila melanogaster* **dm6** reference genome.
The alignment is executed in **two stages**:

1. **Primary alignment** using **HISAT2**
2. **Realignment of unmapped reads** using **Bowtie2**, followed by merging and coordinate sorting

This workflow is optimized for Quant-Seq / single-end libraries and produces BAM files ready for downstream quantification (e.g., with metaseqR2).

---

### **1. Overview of the Pipeline**

For each FASTQ file, the script:

1. Aligns reads with **HISAT2**
2. Collects unmapped reads and realigns them with **Bowtie2**
3. Converts all outputs into BAM format
4. Normalizes Bowtie2 headers to match HISAT2 (required for merging)
5. Merges HISAT2 and Bowtie2 BAM files
6. Performs coordinate sorting and indexing with **samtools**
7. Cleans intermediate files

Final output per sample:

```
sample_name/
    sample_name.bam
    sample_name.bam.bai
```

---

### **2. Requirements**

#### **Software**

You must have the following tools installed:

| Tool     | Version (example) | Used for                                        |
| -------- | ----------------- | ----------------------------------------------- |
| HISAT2   | 2.1.0             | Primary alignment                               |
| Bowtie2  | 2.x               | Secondary alignment of unmapped reads           |
| Samtools | ≥ 1.9             | BAM conversion, sorting, merging, indexing      |
| Bedtools | Optional          | (Included as variable, not used in this script) |

#### **Reference files**

You must provide:

* HISAT2 index files for **dm6**
* Transcriptome splice site file (optional)
* Bowtie2 genome index for **dm6**

---

### **3. Directory Structure**

Expected layout:

```
HOME_PATH/
    fastq_files/
        sample1.fastq
        sample2.fastq
        ...
        hisat_files/        # Created automatically
```

**Output directory:**
`fastq_files/hisat_files/<sample_name>/`

---

### **4. Script Description**

The script contains the following main components:

#### **A. Directory and tool definitions**

Sets the paths to:

* Input FASTQ files
* Output directory
* HISAT2 / Bowtie2 / Samtools / Bedtools executables
* Genome indices

### **B. Loop through FASTQ files**

For each file:

1. Derives sample name
2. Creates a sample-specific output folder

### **C. HISAT2 alignment**

Produces:

* `hisat2.sam`
* `unmapped.fastq` (unmapped reads)

### **D. Bowtie2 alignment of HISAT2-unmapped reads**

Bowtie2 output is:

* Converted directly to BAM (`unmapped_remap_tmp.bad`)

### **E. Header correction**

Bowtie2 and HISAT2 headers differ; script enforces HISAT2 headers for downstream compatibility.

### **F. Merging, sorting, and indexing**

Creates final:

```
sample_name.bam
sample_name.bam.bai
```

### **G. Cleanup**

Intermediate SAM, BAM, and FASTQ files are removed.


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
# **ATAC-seq analysis**
