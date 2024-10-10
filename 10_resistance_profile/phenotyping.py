import pandas as pd
import pysam
import os


def recover_annotation(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    annotations = []

    for record in vcf:
        position = record.pos
        ref = record.ref
        alt_alleles = [alt for alt in record.alts if alt != '<NON_REF>']  # Remove NON_REF
        zygosity = record.samples[0]['GT']  # Genótipo do primeiro indivíduo

        if 'ANN' in record.info:
            for annotation in record.info['ANN']:
                ann_fields = annotation.split('|')
                gene_name = ann_fields[3] if ann_fields[3] else "NA"
                annotation_genic = ann_fields[9] if ann_fields[9] else "NA"
                annotation_protein = ann_fields[10] if ann_fields[10] else "NA"

                nt_change = f"{gene_name}_{annotation_genic}" if gene_name != "NA" and annotation_genic != "NA" else "NA"
                aa_change = f"{gene_name}_{annotation_protein}" if gene_name != "NA" and annotation_protein != "NA" else "NA"

                # Para cada alelo alternativo, criamos uma linha separada
                for alt in alt_alleles:
                    annotations.append([position, ref, alt, nt_change, aa_change, zygosity])
        else:
            for alt in alt_alleles:
                annotations.append([position, ref, alt, "NA", "NA", zygosity])

    # Criar DataFrame com os dados coletados
    df_annotations = pd.DataFrame(annotations, columns=['position', 'ref', 'alt', 'nt_change', 'aa_change', 'zygosity'])

    return df_annotations

def normalize_values(df):
    # Remover espaços em branco e transformar em maiúsculas para evitar erros de formatação
    df['ref'] = df['ref'].str.strip().str.upper()
    df['alt'] = df['alt'].str.strip().str.upper()
    # Garantir que a coluna 'position' seja do tipo inteiro
    df['position'] = df['position'].astype(int)
    return df

def matching(df_annotations, coord_file_path, master_file_path):
    # Carregar o genomic coordinates
    coord_file = pd.read_csv(coord_file_path)
    
    # Carregar o master_file
    master_file = pd.read_csv(master_file_path)
    
    # Normalizar os valores no df_annotations e coord_file
    df_annotations = normalize_values(df_annotations)
    
    # Garantir que a coluna 'position' no coord_file seja do tipo inteiro
    coord_file['position'] = coord_file['position'].astype(int)
    coord_file['reference_nucleotide'] = coord_file['reference_nucleotide'].str.strip().str.upper()
    coord_file['alternative_nucleotide'] = coord_file['alternative_nucleotide'].str.strip().str.upper()
    
    # Realizar a anotação preliminar com o coord_file usando merge
    df_annotations = df_annotations.merge(
        coord_file[['position', 'reference_nucleotide', 'alternative_nucleotide', 'variant']],
        left_on=['position', 'ref', 'alt'],
        right_on=['position', 'reference_nucleotide', 'alternative_nucleotide'],
        how='left'
    )
    
    # Renomear a coluna 'variant' da coord_file para 'master_change'
    df_annotations.rename(columns={'variant': 'master_change'}, inplace=True)

    # Agora prosseguir com o match no master_file usando três critérios:
    # 1. master_change vs variant
    # 2. nt_change vs variant
    # 3. aa_change vs variant
    
    # Primeiro match: Comparar master_change com a coluna variant
    master_match = pd.merge(df_annotations, master_file, left_on="master_change", right_on="variant", how="left")
    
    # Segundo match: Comparar nt_change com a coluna variant
    nt_match = pd.merge(df_annotations, master_file, left_on="nt_change", right_on="variant", how="left")
    
    # Terceiro match: Comparar aa_change com a coluna variant
    aa_match = pd.merge(df_annotations, master_file, left_on="aa_change", right_on="variant", how="left")
    
    # Combinar os três resultados
    combined_df = pd.concat([master_match, nt_match, aa_match])
    
    # Remover duplicatas (caso uma variante tenha dado match em mais de um critério)
    combined_df.drop_duplicates(inplace=True)
    
    # Manter apenas as linhas onde houve interseção (match)
    matched_df = combined_df.dropna(subset=["drug", "variant", "tier", "effect", "FINAL CONFIDENCE GRADING"])
    
    # Selecionar as colunas no formato desejado
    resistance_df = matched_df[['position', 'ref', 'alt', 'nt_change', 'aa_change', 'zygosity', 'drug', 'variant', 'tier', 'effect', 'FINAL CONFIDENCE GRADING', 'master_change']]
    
    return resistance_df


def get_filter_status_from_excel(filter_df, chrom, pos, ref, alt):
    """
    Verifica se uma combinação de cromossomo, posição, ref e alt está presente no dataframe do filter.xlsx.
    Se encontrada, retorna o valor da coluna FILTER. Caso contrário, retorna 'Not Found'.
    """
    matching_row = filter_df[(filter_df['#CHROM'] == chrom) & (filter_df['POS'] == pos) &
                             (filter_df['REF'] == ref) & (filter_df['ALT'] == alt)]
    
    if not matching_row.empty:
        return matching_row.iloc[0]['FILTER']
    return 'Not Found'

def filter_passed_criteria(record):
    """
    Aplica os filtros QD, QUAL e DP. Calcula QD manualmente como QUAL / DP.
    """
    filter_passed = "PASS"

    # Obter DP e QUAL, calculando QD manualmente
    DP = record.info.get("DP", 1)  # Use 1 para evitar divisão por zero se DP não estiver presente
    QUAL = record.qual if record.qual is not None else 0
    QD = QUAL / DP if DP > 0 else 0

    # Aplicar os filtros para QD, QUAL e DP
    if QD < 2.0:
        filter_passed = "QD2"
    elif QUAL < 30.0:
        filter_passed = "QUAL30"
    elif DP < 10:  # Filtro mínimo de profundidade de cobertura (ajuste conforme necessário)
        filter_passed = "DP10"

    return filter_passed


def filtering(resistance_df, filter_excel_path, vcf_file_path):
    # Carregar o arquivo filter.xlsx
    filter_df = pd.read_excel(filter_excel_path)

    # Abrir o VCF correto com base no caminho fornecido
    vcf = pysam.VariantFile(vcf_file_path)

    filter_status = []
    for index, row in resistance_df.iterrows():
        chrom = "NC_000962.3"  # Ajustar conforme necessário
        pos = row['position']
        ref = row['ref']
        alt = row['alt']

        # Verificar se a variante está presente no filter.xlsx
        excel_filter = get_filter_status_from_excel(filter_df, chrom, pos, ref, alt)

        if excel_filter == "Not Found":
            # Se a variante não foi encontrada no arquivo Excel, aplicar os critérios de filtro
            filter_passed = "PASS"
            
            # Encontrar a variante no arquivo VCF para extrair as informações de qualidade
            for record in vcf.fetch(chrom, pos - 1, pos):  # Ajuste para 0-based
                if record.pos == pos and record.ref == ref and alt in record.alts:
                    # Aplicar os critérios de filtragem baseados em QD, QUAL, e DP
                    filter_passed = filter_passed_criteria(record)
                    filter_status.append(filter_passed)
                    break
        else:
            # Se foi encontrada, usar o filtro do arquivo Excel
            filter_status.append(excel_filter)

    # Adicionar a nova coluna de filtro ao dataframe
    resistance_df['Filter_Status'] = filter_status
    
    return resistance_df

def genPhenotyping(vcf_file):

    master_file_path = 'db/tbdr_catalogue_master_file.csv'
    coord_file_path = 'db/tbdr_genomic_coordinates.csv'
    filter_excel_path = 'filter/filter.xlsx'

    # Recuperar as anotações do VCF
    df_annotations = recover_annotation(vcf_file)
    
    # Substituir apenas a extensão .vcf.gz por _ANN.xlsx para o arquivo de anotações
    ann_output_file = vcf_file.replace('.vcf.gz', '_ANN.xlsx')
    
    # Criar a pasta 'results' se ela ainda não existir
    results_dir = os.path.join(os.path.dirname(vcf_file), 'results')
    os.makedirs(results_dir, exist_ok=True)
    
    # Criar o caminho completo para salvar o arquivo _ANN.xlsx na pasta 'results'
    ann_output_path = os.path.join(results_dir, os.path.basename(ann_output_file))
    
    # Salvar o DataFrame df_annotations em um arquivo Excel (.xlsx)
    df_annotations.to_excel(ann_output_path, index=False, sheet_name='annotations')
    
    print(f"Anotações salvas como {ann_output_path}")
    
    # Executar a função de matching por anotação
    resistance_df = matching(df_annotations, coord_file_path, master_file_path)
    
    # Aplicar a função de filtro, passando o caminho do vcf_file como argumento
    final_df = filtering(resistance_df, filter_excel_path, vcf_file)

    # Substituir apenas a extensão .vcf.gz por .xlsx para o arquivo de resultados finais
    output_file = vcf_file.replace('.vcf.gz', '_target.xlsx')
    
    # Criar o caminho completo para salvar o arquivo final na pasta 'results'
    output_path = os.path.join(results_dir, os.path.basename(output_file))
    
    # Salvar o resultado em um novo arquivo Excel (.xlsx)
    final_df.to_excel(output_path, index=False, sheet_name='resistance')
    
    print(f"Arquivo final salvo como {output_path}")

    return output_path
