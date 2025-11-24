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
    echo "     ===== Converting BAM file $MAPPED to BED file..."
    $BEDTOOLS/bedtools bamtobed -i $MAPPED > $DIR/$BASE".tmp.bed"
    echo "     ===== Cleaning the BED file..."
    printf "%d\t" `sort -k1,1 -k2g,2 -u $DIR/$BASE".tmp.bed" | wc -l` >> $DIR/bowtie2_report.txt
    echo "     ===== Creating the final BED file..."
    grep -vP 'chrG|chrM|chrU|rand|hap|map|cox|loc|chrEBV' $DIR/$BASE".tmp.bed" | sort -k1,1 -k2g,2 -u | pigz -p $THREADS > $DIR/$BASE".bed.gz"
    rm $DIR/$BASE".tmp.bed" $DIR/$BASE".uns"
    printf "%d\n" `zcat $DIR/$BASE".bed.gz" | wc -l` >> $DIR/bowtie2_report.txt
done

#Index bam files
samtools index -@ 25 *.bam
