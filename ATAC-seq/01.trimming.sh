for i in /media/samba/hatzis_lab/delidakis_atac/fastq_files/*.gz; do /media/raid/resources/ngstools/trimgalore/trim_galore --max_n 5 --fastqc $i; done
