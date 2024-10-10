#!/usr/bin/env python3

import pysam
import pandas as pd
import glob

# Caminho para os arquivos VCF
data = []  # Lista para armazenar os dados finais

for file in glob.glob('*.gz'):
    # Abrir o arquivo VCF
    vcf = pysam.VariantFile(file)
    sample_name = file.split('-sorted')[0]  # Ajuste conforme o padrão de nome do arquivo

    found_position = False  # Flag para verificar se a posição foi encontrada
    af = 0  # Frequência alélica inicial

    # Iterar sobre cada registro no arquivo VCF
    for record in vcf.fetch(region="NC_000962.3"):
        if record.pos == 1472307:  # Verifica apenas essa posição de interesse
            found_position = True
            if 'AF' in record.info:
                af = record.info['AF'][0]  # Assume que AF está presente e pega o primeiro valor
            break

    data.append([sample_name, af])  # Armazenar nome da amostra e frequência alélica

# Criar DataFrame e salvar em Excel
df = pd.DataFrame(data, columns=['Sample Name', 'NTM <20% 1472307'])
df.to_excel('allelic_frequency.xlsx', index=False)

