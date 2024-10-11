#!/bin/bash

# Diretório com os arquivos BAM
bam_dir="run/map"
gvcf_dir="/run/map/gvcf"

# Cria o diretório gvcf se ele não existir
mkdir -p "$gvcf_dir"

# Loop sobre os arquivos BAM
for bam_file in "${bam_dir}"/*-sorted-named-dupl.bam; do
    # Extrair o nome do arquivo sem a extensão
    filename=$(basename "$bam_file")
    sample_name="${filename%-sorted-named-dupl.bam}"
    
    # Nome do arquivo de saída para o gVCF
    output_gvcf="${gvcf_dir}/${sample_name}.g.vcf.gz"
    
    # Comando para chamar as variantes com GATK HaplotypeCaller em modo GVCF, incluindo MNPs
    /java -jar /gatk-package-4.4.0.0-local.jar HaplotypeCaller \
    -R NC0009623.fasta \
    -I "$bam_file" \
    -O "$output_gvcf" \
    -ERC GVCF \
    --max-mnp-distance 1 \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A QualByDepth \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest \
    -A DepthPerAlleleBySample \
    -A Coverage
    
    echo "gVCF gerado para $sample_name"
done

