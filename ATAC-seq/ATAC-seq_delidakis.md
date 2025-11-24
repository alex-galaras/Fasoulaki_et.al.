**Step I: Trimming**

    for i in /media/samba/hatzis_lab/delidakis_atac/fastq_files/*.gz; do /media/raid/resources/ngstools/trimgalore/trim_galore --max_n 5 --fastqc $i; done

**Step II: Alignment**

    # Paths and vars
    BOWTIE2_PATH=/media/raid/resources/ngstools/bowtie2
    FASTQ_PATH=/media/samba/hatzis_lab/delidakis_atac/fastq_trimmed
    DIR=/media/samba/hatzis_lab/delidakis_atac/bam_files
    BOWTIE2_INDEX=/media/raid/resources/genomes/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome
    SAMTOOLS=/media/raid/resources/ngstools/samtools
    BEDTOOLS=/media/raid/resources/ngstools/bedtools/bin/
    THREADS=8

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "name" "total" "mapped" "qc_filtered" "unique" "clean" > $DIR/bowtie2_report.txt
    mkdir -p $DIR
    for FILE in $FASTQ_PATH/*.fq.gz
    do
        BASE=`basename $FILE | sed s/\.fq.gz//`
        echo "========== Processing $BASE..."   
        echo "     ===== Mapping..."
        $BOWTIE2_PATH/bowtie2 -p $THREADS --local --very-sensitive-local --mm -x $BOWTIE2_INDEX -U $FASTQ_PATH/$BASE".fq.gz" | $SAMTOOLS/samtools view -bS -o $DIR/$BASE".uns" -
        echo "     ===== Sorting..."
        $SAMTOOLS/samtools sort -o $DIR/$BASE".bam" $DIR/$BASE".uns"
        echo "     ===== Processing the output..."
        MAPPED=$DIR/$BASE".bam"
        printf "%s\t" $BASE >> $DIR/bowtie2_report.txt
        echo "     ===== Counting total reads..."
        printf "%d\t" `zcat $FASTQ_PATH/$BASE".fastq.gz" | wc -l | awk '{print $1/4}'` >> $DIR/bowtie2_report.txt
        echo "     ===== Counting mapped reads..."
        printf "%d\t" `$SAMTOOLS/samtools view -c -F 4 $MAPPED` >> $DIR/bowtie2_report.txt
        echo "     ===== Counting mapped and high-quality reads..."
        printf "%d\t" `$SAMTOOLS/samtools view -c -Fsamtoo $MAPPED` >> $DIR/bowtie2_report.txt
        rm $DIR/$BASE".tmp.bed" $DIR/$BASE".uns"
        printf "%d\n" `zcat $DIR/$BASE".bed.gz" | wc -l` >> $DIR/bowtie2_report.txt
    done

    #Index bam files
    samtools index -@ 25 *.bam

**Step III: Remove duplicates**

    DIR=/media/samba/hatzis_lab/delidakis_atac/bam_files
    CLEAN=/media/samba/hatzis_lab/delidakis_atac/bam_clean_files
    PICARD=/media/raid/resources/ngstools/picard/picard.jar
    SAMTOOLS=`which samtools`
    GENOME=/media/raid/resources/genomes/igenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa

    mkdir -p $CLEAN

    for FILE in $DIR/*.bam
    do
    BASE=$(basename "$FILE" ".bam")
    $SAMTOOLS view -h "$FILE" | \
    grep -vP "chrG|chrM|chrU|rand|hap|map|cox|loc|chrEBV" | \
    samtools view -@ 16 -q 20 -bT "$GENOME" > "$CLEAN/${BASE}_cleaned.bam" 
    done

    $SAMTOOLS index -@ 16 -M "$CLEAN"/*.bam

    for FILE in $CLEAN/*.bam
    do
    BASE=$(basename "$FILE" "_cleaned.bam")
    java -jar "$PICARD" MarkDuplicates \
        I="$FILE" \
        O="$CLEAN/${BASE}_deduplicated.bam" \
        M="${BASE}_metrics.txt" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true
    done

**Generate normalized tracks**

    for bam in /media/samba/hatzis_lab/delidakis_atac/bam_clean_files/*_trimmed_deduplicated.bam
    do
      bamCoverage --bam $bam --normalizeUsing CPM --extendReads 200 -o ${bam/_trimmed_deduplicated.bam/.bw}
    done

Create track links

    # Output file
    OUTFILE="/media/samba/hatzis_lab/delidakis_atac/track.txt"

    # Overwrite existing file
    > "$OUTFILE"

    # Loop through all .bw files
    for FILE in /media/samba/hatzis_lab/delidakis_atac/bam_clean_files/*.bw; do
        SAMPLE=$(basename "$FILE" .bw)
        echo -e "\ntrack type=bigWig name=${SAMPLE} bigDataUrl=http://epigenomics.fleming.gr/~alexandros/delidakis_ATAC/${SAMPLE}.bw color=0,0,0" >> "$OUTFILE"
    done

**Peak calling**

    BAM_DIR=/media/samba/hatzis_lab/delidakis_atac/bam_clean_files
    OUT_DIR=/media/samba/hatzis_lab/delidakis_atac/peaks
    GENOME_SIZE=dm  # effective genome size for dm6 (Drosophila)

    mkdir -p "$OUT_DIR"

    for bam in "$BAM_DIR"/*_deduplicated.bam
    do
      BASE=$(basename "$bam" "_deduplicated.bam")
      echo "==> Calling peaks for $BASE ..."

      macs3 callpeak \
        -t "$bam" \
        -f BAM \
        -g "$GENOME_SIZE" \
        -n "$BASE" \
        -B \
        -q 0.01 \
        --nomodel \
        --shift -100 --extsize 200 \
        --call-summits \
        --outdir "$OUT_DIR"
    done

**Differential accessibility analysis**

    # Load required libraries
    library("DiffBind")  # For ATAC-seq differential binding analysis
    library("ggplot2")   # For plotting


    # ==============================================================
    # DiffBind Automated DESeq2 Analysis for All Pairwise Comparisons
    # with Safe Handling for Empty Contrasts
    # ==============================================================

    # DiffBind Automated DESeq2 Analysis for All Pairwise Comparisons

    # Load required libraries
    library("DiffBind")
    library("ggplot2")

    # === Global Parameters ===
    template <- "/media/samba/hatzis_lab/delidakis_atac/diffbind/templates/template_broad.csv"
    base_dir <- "/media/samba/hatzis_lab/delidakis_atac/diffbind"
    output_root = "/media/samba/hatzis_lab/delidakis_atac/diffbind/results_broad/"

    dir.create(output_root, showWarnings = FALSE)

    # Define all experimental conditions
    conditions <- c("DSS_trx", "DSS_control", "suc_trx", "suc_control")

    # Set working directory
    setwd(output_root)

    # === Import and preprocess sample data once ===
    samples.ATAC <- dba(sampleSheet = template)
    samples.ATAC <- dba.blacklist(samples.ATAC, blacklist = DBA_BLACKLIST_DM6, greylist = FALSE)

    # Save peaks info
    write.table(dba.show(samples.ATAC), "peaks_information.txt",
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

    # Correlation heatmap for peaks
    tiff("correlation_heatmap_peaks.tiff")
    plot(samples.ATAC)
    dev.off()

    # Count reads (shared across comparisons)
    ATAC.unique_counts <- dba.count(samples.ATAC,
                                    minOverlap = 2,
                                    filter = 10,
                                    bParallel = TRUE,
                                    fragmentSize = 0,
                                    summits=FALSE,
                                    bRemoveDuplicates = FALSE)

    # Save counts info
    write.table(dba.show(ATAC.unique_counts), "counts_information.txt",
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

    # Correlation + PCA plots
    tiff("correlation_heatmap_counting.tiff")
    plot(ATAC.unique_counts)
    dev.off()

    tiff("pca_counting.tiff")
    dba.plotPCA(ATAC.unique_counts, c(DBA_TISSUE, DBA_CONDITION), DBA_ID)
    dev.off()

    #Perform differential accessibility analysis with DESEQ2
    # ------------------------------
    # Differential Accessibility Function with Volcano Plot
    # ------------------------------

    diff.accessibility <- function(condB, condA, 
                                   fc_thresh = 0.5, 
                                   fdr_thresh = 0.05,
                                   output_root = "/media/samba/hatzis_lab/delidakis_atac/diffbind/results_broad/",
                                   dba_object = ATAC.unique_counts) {
      
      # Comparison setup
      comparison_name <- paste0(condB, "_vs_", condA)
      comparison_dir <- file.path(output_root, comparison_name)
      
      # Create output directory if not exists
      if (!dir.exists(comparison_dir)) {
        dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
      }
      
      message("🔹 Running comparison: ", comparison_name)
      
      # ------------------------------
      # Define contrast and analyze
      # ------------------------------
      
      comp <- dba.contrast(ATAC.unique_counts, contrast = c("Condition", condB, condA))
      comp <- dba.analyze(comp, method = DBA_DESEQ2)
      
      # ------------------------------
      # Extract and save results
      # ------------------------------
      
      # All results
      dba.all <- dba.report(comp, method = DBA_DESEQ2, th = 1, bUsePval = FALSE, fold = 0, bCount = TRUE)
      write.table(as.data.frame(dba.all),
                  file = file.path(comparison_dir, paste0(comparison_name, "_all_deseq2.txt")),
                  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      
      # FDR-filtered
      dba.fdr <- dba.report(comp, method = DBA_DESEQ2, th = fdr_thresh, bUsePval = FALSE, fold = 0, bCount = TRUE)
      dba.fdr <- subset(dba.fdr, abs(Fold) > fc_thresh)
      write.table(as.data.frame(dba.fdr),
                  file = file.path(comparison_dir, paste0(comparison_name, "_fdr_deseq2.txt")),
                  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      
      # P-value-filtered
      dba.pval <- dba.report(comp, method = DBA_DESEQ2, th = 0.05, bUsePval = TRUE, fold = 0, bCount = TRUE)
      dba.pval <- subset(dba.pval, abs(Fold) > fc_thresh)
      write.table(as.data.frame(dba.pval),
                  file = file.path(comparison_dir, paste0(comparison_name, "_pval_deseq2.txt")),
                  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      
      # ------------------------------
      # Volcano Plot (FDR-based)
      # ------------------------------
      
      dba.all$Significant <- ifelse(dba.all$FDR <= fdr_thresh & dba.all$Fold > fc_thresh, "Upregulated",
                                    ifelse(dba.all$FDR <= fdr_thresh & dba.all$Fold < -fc_thresh, "Downregulated",
                                           "Unregulated"))
      
      g1 <- ggplot(as.data.frame(dba.all), aes(x = Fold, y = -log10(FDR), color = Significant)) +
        geom_point(alpha = 0.8, size = 1.5) +
        
        ggtitle(paste("Volcano Plot -", comparison_name)) +
        xlab("log2 Fold Change") +
        ylab("-log10(FDR)") +
        
        scale_color_manual(values = c("Upregulated" = "#E41A1C", 
                                      "Downregulated" = "#377EB8", 
                                      "Unregulated" = "black")) +
        
        geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "grey40") +
        geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed", color = "grey40") +
        
        theme_classic(base_family = "Arial") +
        theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 16, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face = "bold"),
          legend.position = "top"
        )
      
      ggsave(filename = file.path(comparison_dir, paste0(comparison_name, "_deseq2_volcano_fdr.png")),
             plot = g1, width = 10, height = 8, dpi = 300)
      
      message("✅ Finished: ", comparison_name)
      
    }

    # ------------------------------
    # Run Comparisons
    # ------------------------------

    diff.accessibility("suc_trx", "suc_control")
    diff.accessibility("DSS_trx", "DSS_control")
    diff.accessibility("DSS_control", "suc_control")
    diff.accessibility("DSS_trx", "suc_trx")

\#Peak annotation based on everything that resides 5kb around the peak
center

\##Prepare dm6 genome

    wget https://ftp.ensembl.org/pub/release-115/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.115.gtf.gz
    DIR=/media/samba/hatzis_lab/delidakis_atac/peak_annotation/bedtools_approach

    awk 'BEGIN{OFS="\t"}{
        if ($1 !~ /^chr/) $1 = "chr"$1;
        print
    }' /media/raid/users/agalaras/genomes/Drosophila_melanogaster.BDGP6.54.115.gtf | tail +5 > /media/raid/users/agalaras/genomes/dm6_genes.BDGP6.54.115.gtf

    #Gene gtf
    GTF=/media/raid/users/agalaras/genomes/dm6_genes.BDGP6.54.115.gtf
    grep -w "gene" $GTF | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$4-1,$5,$10","$12}' | sort -k1,1 -k2,2n > ${DIR}/dm6_genes.bed

    #TSS gtf
    awk 'BEGIN{OFS="\t"}
         $3 == "gene" {
             if ($7 == "+") {
                 tss = $4 - 1;   # BED is 0-based
             } else if ($7 == "-") {
                 tss = $5 - 1;
             } else {
                 next;
             }
             print $1, tss, tss+1, $10","$12
         }' /media/raid/users/agalaras/genomes/dm6_genes.BDGP6.54.115.gtf | sort -k1,1 -k2,2n > ${DIR}/tss.bed

\##Process peak file

    DIR=/media/samba/hatzis_lab/delidakis_atac/peak_annotation/bedtools_approach
    TSS=${DIR}/tss.bed
    GENEBODY=${DIR}/dm6_genes.bed

    for FILE in /media/samba/hatzis_lab/delidakis_atac/diffbind/results_broad/*/*_fdr_deseq2.txt
    do
    SAMPLE=$(basename $FILE "_fdr_deseq2.txt")
    #TSS
    tail +2 $FILE | \
    awk 'BEGIN{OFS="\t"}{
        orig_id = $1 ":" $2 "-" $3      # chr:start-end
        center = int(($2+$3)/2)
        start  = center - 5000
        end    = center + 5000
        print $1, start, end, orig_id, $9, $10
    }' | sort -k1,1 -k2,2n  > ${DIR}/${SAMPLE}_5kb_window.bed

    #Intersection to find TSSs inside 5kb from the peak center
    bedtools intersect -a ${DIR}/${SAMPLE}_5kb_window.bed -b $TSS -wa -wb | cut -f4-10 | sort -k1,1 > ${DIR}/${SAMPLE}_overlap_tss.bed

    #Gene body
    tail +2 $FILE | \
    awk 'BEGIN{OFS="\t"} {
    orig_id = $1 ":" $2 "-" $3      # chr:start-end
    print $1, $2, $3, orig_id
    }' | sort -k1,1 -k2,2n  > ${DIR}/${SAMPLE}_peaks.bed

    #Intersection to find overlapping genes with the peak
    bedtools intersect -a ${DIR}/${SAMPLE}_peaks.bed -b $GENEBODY -wa -wb | cut -f4-10 | sort -k1,1 > ${DIR}/${SAMPLE}_overlap_genebody.bed

    join -t $'\t' -1 1 -2 1 ${DIR}/${SAMPLE}_overlap_tss.bed ${DIR}/${SAMPLE}_overlap_genebody.bed > ${SAMPLE}_merged.bed
    rm ${DIR}/${SAMPLE}_overlap_tss.bed ${DIR}/${SAMPLE}_overlap_genebody.bed ${DIR}/${SAMPLE}_peaks.bed ${DIR}/${SAMPLE}_5kb_window.bed
    done
