#!/usr/bin/env python3

import pandas as pd
import time
import gc  # Garbage collector

# Clock start
#start_time = time.time()

def _genNrPositions(genomic_coordinates):

    # MNPs (where len(ref) == len(alt) and len(ref) > 1) for nr_positions and range_positions files
    mnps = genomic_coordinates[
        (genomic_coordinates['reference_nucleotide'].str.len() == genomic_coordinates['alternative_nucleotide'].str.len()) &
        (genomic_coordinates['reference_nucleotide'].str.len() > 1)
    ].copy()
    
    # Expand MNPs into individual positions
    mnps_positions = mnps.apply(
        lambda row: pd.Series(range(row['position'], row['position'] + len(row['reference_nucleotide']))),
        axis=1
    ).stack().reset_index(drop=True)

    # SNPs and indels for nr_positions and range_positions files
    snps_indels = genomic_coordinates[
        (genomic_coordinates['reference_nucleotide'].str.len() != genomic_coordinates['alternative_nucleotide'].str.len()) | 
        ((genomic_coordinates['reference_nucleotide'].str.len() == genomic_coordinates['alternative_nucleotide'].str.len()) & 
        (genomic_coordinates['reference_nucleotide'].str.len() == 1))
    ].copy()
    snps_indels_positions = snps_indels['position']

    # Combine MNPs and SNPs/indels into a single list of positions for nr_positions file
    positions_of_interest = pd.concat([mnps_positions, snps_indels_positions], ignore_index=True).drop_duplicates()

    # Convert positions to integers to avoid floating point numbers
    positions_of_interest = positions_of_interest.astype(int)
    
    # Sort by 'position' before saving
    positions_of_interest = positions_of_interest.sort_values()

    # Criar um bed file com intervalos de posições para não ficarem muitas linhas
    bed_intervals = []
    start = positions_of_interest.iloc[0]
    end = start

    for pos in positions_of_interest.iloc[1:]:
        if pos == end + 1:
            end = pos
        else:
            if start == end:
                end = start + 1  # Ajustar para evitar start e end iguais
            bed_intervals.append(['NC_000962.3', start, end])
            start = pos
            end = start
    
    # Adicionar o último intervalo
    if start == end:
        end = start + 1  # Ajustar para evitar start e end iguais
    bed_intervals.append(['NC_000962.3', start, end])
    
    # Converter para DataFrame
    bed_df = pd.DataFrame(bed_intervals, columns=['chrom', 'start', 'end'])
    
    # Salvar no formato BED
    bed_df.to_csv('db/tbdr_nr_positions.bed', sep='\t', header=False, index=False)

    # Clean up memory
    del genomic_coordinates, mnps, snps_indels, positions_of_interest, bed_df
    gc.collect()

    return

def genDbFormat():
    
    # Open TB drugs excel database
    genomic_coordinates = pd.read_excel('db/WHO-UCN-TB-2023.7-eng.xlsx', sheet_name='Genomic_coordinates')

    _genNrPositions(genomic_coordinates)

    # Save the DataFrames to new CSV files
    genomic_coordinates.to_csv('db/tbdr_genomic_coordinates.csv', index=False)

    # Clean up memory
    del genomic_coordinates
    gc.collect()

    master_file = pd.read_excel('db/WHO-UCN-TB-2023.7-eng.xlsx', sheet_name='Catalogue_master_file', header=2)
    master_file = master_file[['drug', 'variant', 'tier', 'effect', 'FINAL CONFIDENCE GRADING']]
    
    # Save the DataFrames to new CSV files
    master_file.to_csv('db/tbdr_catalogue_master_file.csv', index=False)

    # Clean up memory
    del master_file
    gc.collect()

    return


#end_time = time.time()
#print(f"Execution time dbFormat: {(end_time - start_time) / 60:.2f} minutes")
