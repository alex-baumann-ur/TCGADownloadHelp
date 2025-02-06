### Run the pipeline with the following command in the directory with the Snakefile in a Snakemake conda environment: snakemake

import pandas as pd
import os
import yaml
import time

start_time = time.strftime('%Y-%m-%d_%H:%M:%S', time.localtime())

# Read in file paths and specifications from config file
with open('data/config.yaml', 'r') as streamfile:
    conf_yaml = yaml.load(streamfile, Loader=yaml.FullLoader)

analysis_path = conf_yaml['analysis_path']

# Create log file folders
os.makedirs(analysis_path+'log_files', exist_ok=True)
os.makedirs(analysis_path+'log_files/correct_files', exist_ok=True)

rule all:
    input:
        analysis_path + f'log_files/correct_files/{start_time}_config_file_correct.txt',
        analysis_path + f'log_files/correct_files/{start_time}_combine_manifest_sample_sheet_correct.txt',
        analysis_path + f'log_files/correct_files/{start_time}_download_TCGA_data_correct.txt',
        analysis_path + f'log_files/correct_files/{start_time}_rename_files_correct.txt'


rule CheckConfigFile:
    input:
        conf = 'data/config.yaml'
    output:
        out_correct_file = analysis_path + f'log_files/correct_files/{start_time}_config_file_correct.txt',
    log:
        log_file = analysis_path + f'log_files/{start_time}_check_config_file_log.txt'
    shell:
        'python3 scripts_TCGA_pipeline/01_check_config_file.py {input.conf} {output.out_correct_file} {log.log_file}'


rule CombineManifestSampleSheet:
    input:
        conf = 'data/config.yaml',
        prev_files = analysis_path + f'log_files/correct_files/{start_time}_config_file_correct.txt'
    output:
        out_correct_file = analysis_path + f'log_files/correct_files/{start_time}_combine_manifest_sample_sheet_correct.txt',
    log:
        log_file = analysis_path + f'log_files/{start_time}_combine_manifest_sample_sheet_log.txt'
    shell:
        'PYTHONPATH=scripts_TCGA_pipeline python3 scripts_TCGA_pipeline/02_combine_manifest_sample_sheet.py {input.conf} {output.out_correct_file} {log.log_file}'


rule DownloadTCGAData:
    input:
        conf = 'data/config.yaml',
        prev_files = analysis_path + f'log_files/correct_files/{start_time}_combine_manifest_sample_sheet_correct.txt',
    output:
        out_correct_file = analysis_path + f'log_files/correct_files/{start_time}_download_TCGA_data_correct.txt',
    log:
        log_file = analysis_path + f'log_files/{start_time}_download_TCGA_data_log.txt'
    shell:
        'python3 scripts_TCGA_pipeline/03_download_TCGA_data.py {input.conf} '\
        f'{analysis_path}log_files/{start_time}_gdc_client_log '\
        '{output.out_correct_file} {log.log_file}'


rule RenameFiles:
    input:
        conf = 'data/config.yaml',
        prev_files = analysis_path + f'log_files/correct_files/{start_time}_download_TCGA_data_correct.txt'
    output:
        out_correct_file = analysis_path + f'log_files/correct_files/{start_time}_rename_files_correct.txt'
    log:
        log_file = analysis_path + f'log_files/{start_time}_rename_files_log.txt'
    shell:
        'python3 scripts_TCGA_pipeline/04_rename_files.py {input.conf} {output.out_correct_file} {log.log_file}'
