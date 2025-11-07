#Run FastQC on all FASTQ Files in the folder
for i in *.gz; do fastqc $i; done

#Optional MultiQC
multiqc .
