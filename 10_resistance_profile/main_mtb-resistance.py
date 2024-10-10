#!/usr/bin/python3 -u

import time
import os
import glob

from vcf_processing import genProcessVcf
from dbFormat import genDbFormat
from phenotyping import genPhenotyping

# Clock start
start_time = time.time()


#Database path
path_db = 'db'

#'tbdr_nr_positions.bed',

#Check if csv files from db exist
expected_files = [
    'tbdr_catalogue_master_file.csv',
    'tbdr_genomic_coordinates.csv'
]

missing_files = [file for file in expected_files if not os.path.exists(os.path.join(path_db, file))]

if missing_files:
    print("Preparing database...")
    genDbFormat()
else:
    print('Database already formated.')

# Run - VCF files from gatk
path_snpeff = '/home/nailasoler/Lapam/IAL/4_0_bcftools_norm/snpeff_norm205'

#Pre-processing annotated VCFs
path_filtered = os.path.join(path_snpeff, "filtered")
os.makedirs(path_filtered, exist_ok=True)
genProcessVcf(path_snpeff, path_filtered)

#Matching DB patient
for vcf_file in glob.glob(f'{path_filtered}/*.gz'):

    print(vcf_file, flush=True)
    genPhenotyping(vcf_file)


end_time = time.time()
print(f"Execution time full pipeline: {(end_time - start_time) / 60:.2f} minutes")
