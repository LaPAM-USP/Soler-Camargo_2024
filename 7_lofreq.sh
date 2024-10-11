#!/bin/bash

# Indelqual
for bamfile in *-sorted-named-dupl.bam; do 
    lofreq indelqual --dindel -f NC0009623.fasta -o "${bamfile%.}_lofreq.bam" "$bamfile"
done

# Aln qual
for bamfile in *_lofreq.bam; do 
    output_aln_qual="${bamfile%.*}_aln_qual.bam"
    lofreq alnqual -b -r "$bamfile" NC0009623.fasta > "$output_aln_qual"
done

# Variant calling
for aln_qual_file in *_aln_qual.bam; do 
    output_vcf="${aln_qual_file%.*}_LofreqCall.vcf"
    lofreq call -f NC0009623.fasta -o "$output_vcf" --call-indels "$aln_qual_file"
done
