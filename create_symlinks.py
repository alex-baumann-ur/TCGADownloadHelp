import os
import pandas as pd

# sample_sheet = pd.read_table('sample_sheets/2024-03-25_sample_sheet_lusc_paths_dna.tsv') # for DNA
sample_sheet = pd.read_table('sample_sheets/2024-04-08_sample_sheet_lusc_paths_rna.tsv') # for RNA

path_sample_data = '/media/md0/abaumann/01_OLCIR/genomic_instab_test_lusc/01_sample_data/'
path_raw_data = '/media/md0/abaumann/01_OLCIR/genomic_instab_test_lusc/00_raw_data/'

case_names = list(sample_sheet['Case ID'].unique())
# methods = ['BRASS','CaVEMan','CNV_segment','Pindel'] # for DNA
methods = ['Splicing', 'STAR_counts']

# create folders for methods
for i in methods:
    os.makedirs(path_sample_data+i, exist_ok=True)

# create symlinks for samples to categorize them in 
for raw_path, sample_path in zip(list(sample_sheet['Path_raw']), list(sample_sheet['Path_sample'])):
    try:
        os.symlink(raw_path, sample_path)
    except FileExistsError:
        print(sample_path + ' already exists.')




