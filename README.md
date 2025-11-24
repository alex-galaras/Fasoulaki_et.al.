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
  SQLite annotation DB for dm6 (local reference database). It can be downloaded from [here](https://drive.google.com/drive/folders/15lOY9PBggCcaoohO_0rQTvExXenqah55) as suggested in the metaseqR2 package

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

This repository contains the complete workflow used to process, align, QC, and analyze ATAC-seq data for the Delidakis project, including trimming, alignment, duplicate removal, peak calling, differential accessibility analysis, peak annotation, and visualization.

The workflow is implemented in an R Markdown file and includes both **bash** and **R** code, enabling full reproducibility.

---

# **Pipeline Overview**

The ATAC-seq analysis is performed in the following major steps:

1. **Read trimming** (Trim Galore! + FastQC)
2. **Alignment** (Bowtie2 + Samtools)
3. **Duplicate removal + mitochondrial filtering**
4. **Track generation** (bigWig normalized tracks using CPM)
5. **Peak calling** (MACS3)
6. **Differential accessibility** (DiffBind + DESeq2)
7. **Peak annotation** (TSS and gene-body intersections using bedtools)

All paths and file structures correspond to the configuration used on the Hatzis Lab servers.

---

## **1. Input Requirements**

### **Raw Data**

* Paired or single-end gzipped FASTQ files:

  ```
  /media/samba/hatzis_lab/delidakis_atac/fastq_files/*.gz
  ```

### **Reference Files**

* Bowtie2 index for **dm6**
* dm6 genome FASTA
* dm6 Ensembl GTF (BDGP6.54.115)
* Blacklist regions: `DBA_BLACKLIST_DM6` (from DiffBind)

### **Software**

| Tool         | Version              | Used for                   |
| ------------ | -------------------- | -------------------------- |
| Trim Galore! | ≥ 0.6                | Adapter trimming           |
| FastQC       | ≥ 0.11               | QC                         |
| Bowtie2      | ≥ 2.4                | ATAC-seq alignment         |
| Samtools     | ≥ 1.9                | BAM processing             |
| Bedtools     | ≥ 2.25               | Peak annotation            |
| MACS3        | ≥ 3.0                | Peak calling               |
| DiffBind     | ≥ 3.0                | Differential accessibility |
| R ≥ 4.0      | Statistical analysis |                            |

---

# **2. Workflow Summary**

## **Step I — Trimming**

Trim Galore! removes adapters and performs quality trimming:

```bash
for i in /media/samba/hatzis_lab/delidakis_atac/fastq_files/*.gz; do
    /media/raid/resources/ngstools/trimgalore/trim_galore --max_n 5 --fastqc $i
done
```

---

## **Step II — Alignment**

Reads are aligned to the **dm6** genome using Bowtie2 in very-sensitive local mode.
Alignment report includes total reads, mapped reads, and high-quality reads.

Outputs:

* Sorted BAM files
* Alignment summary table (`bowtie2_report.txt`)

---

## **Step III — Duplicate Removal & Filtering**

Reads mapped to:

* chrM
* unplaced scaffolds
* haplotypes
* random chromosomes

are removed.

Duplicates are removed using **Picard MarkDuplicates**.

Outputs:

* *_cleaned.bam
* *_deduplicated.bam
* BAM indexes
* duplication metrics

---

## **Step IV — Normalized Track Generation**

Each cleaned BAM file is converted to bigWig with CPM normalization:

```bash
bamCoverage --bam sample.bam --normalizeUsing CPM --extendReads 200 -o sample.bw
```

Track lines for the UCSC Genome Browser are auto-generated:

```
track type=bigWig name=Sample bigDataUrl=http://...
```

---

## **Step V — Peak Calling (MACS3)**

MACS3 is used with:

* 200 bp extension
* 100 bp shift
* q = 0.01
* summit calling

Outputs:

* NarrowPeak files
* bedGraph & bigWig signal
* Summit files

---

## **Step VI — Differential Accessibility (DESeq2)**

A fully automated pipeline performs all pairwise comparisons:

* suc_trx vs. suc_control
* DSS_trx vs. DSS_control
* DSS_control vs. suc_control
* DSS_trx vs. suc_trx

For each comparison, the pipeline exports:

* All peaks (`_all_deseq2.txt`)
* FDR-filtered peaks (`_fdr_deseq2.txt`)
* P-value-filtered peaks (`_pval_deseq2.txt`)
* Volcano plot images (PNG)

---

## **Step VII — Peak Annotation (±5 kb TSS + Gene Body)**

Peaks are annotated based on:

* TSS within ±5 kb of peak center
* Overlaps with gene bodies

Using bedtools:

* A 5-kb sliding window is generated per peak
* Overlaps with TSS and genes are computed
* Merged annotation tables are produced


# ** Contact**

For questions or issues regarding the workflow:

📧 **[galaras@fleming.gr](mailto:galaras@fleming.gr)**
