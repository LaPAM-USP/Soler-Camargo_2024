import subprocess
import glob
import os
import pandas as pd

def filter_vcf_with_bcftools(patient_vcf, path_filtered):

    # Nome do arquivo de saída
    output_vcf = os.path.join(path_filtered, os.path.basename(patient_vcf).replace('.vcf.gz', '_filtered.vcf.gz'))

    # Comando bcftools para filtrar o VCF só onde tem anotação
    bcftools_command = f"bcftools view -i 'INFO/ANN!=\".\"' {patient_vcf} -Oz -o {output_vcf}"

    # Execute o comando bcftools usando subprocess
    subprocess.run(bcftools_command, shell=True, check=True)

    # Indexar o VCF filtrado usando bcftools
    index_command = f"bcftools index {output_vcf}"
    subprocess.run(index_command, shell=True, check=True)

    return output_vcf

def genProcessVcf(path_gatk, path_filtered):
     
    # Encontrar todos os arquivos VCF no caminho especificado
    patient_vcfs = glob.glob(os.path.join(path_gatk, "*.vcf.gz"))
    
    # Processar cada VCF sequencialmente
    for patient_vcf in patient_vcfs:

        filter_vcf_with_bcftools(patient_vcf, path_filtered)
