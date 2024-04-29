# on Curie: run pipeline with: snakemake --cores 20 --use-conda --conda-prefix /media/md0/abaumann/99_snakemake_tmp
# locally: run pipeline with: snakemake --cores 6 --use-conda --conda-prefix /home/alex/Dokumente/98_snakemake_tmp

import pandas as pd
import os
import yaml

with open('data/config.yaml', 'r') as streamfile:
    conf_yaml = yaml.load(streamfile, Loader=yaml.FullLoader)

an_path = conf_yaml['an_path']
sample_sheet = pd.read_table(conf_yaml['sample_sheet'])

SAMPLES = list(sample_sheet['Case ID'].unique())
METHODS = conf_yaml['methods']

os.makedirs(an_path+'02_results_raw', exist_ok=True)
for m in METHODS:
    os.makedirs(an_path+'02_results_raw/'+m, exist_ok=True)
os.makedirs(an_path+'03_results_combined', exist_ok=True)

rule all:
    input:
        #expand(an_path+'02_results_raw/{method}/{sample}_{method}_out.txt', method=METHODS, sample=SAMPLES)
        expand(an_path+'02_results_raw/BRASS/{sample}_{method}_out.txt', method=METHODS, sample=SAMPLES)


rule BRASS:
    input:
        in_file = an_path+'01_sample_data/BRASS/{sample}.BRASS.raw_structural_variation.vcf.gz',
        bedpe = an_path+'01_sample_data/BRASS/{sample}.BRASS.rerun_structural_variation.bedpe.gz',
    output:
        out_txt = an_path+'02_results_raw/BRASS/{sample}_BRASS_out.txt',
    threads: 2
    conda: 'envs/vcf_env.yaml'
    script:
        'scripts/BRASS_analysis.py'


rule CaVEMan:
    input:
        in_file = an_path+'01_sample_data/CaVEMan/{sample}.CaVEMan.raw_somatic_mutation.vcf.gz',
    output:
        out_txt = an_path+'02_results_raw/CaVEMan/{sample}_CaVEMan_out.txt',
    threads: 2
    script:
        'scripts/CaVEMan_analysis.py'


rule CNV_segment:
    input:
        in_file = an_path+'01_sample_data/CNV_segment/{sample}.ASCAT.copy_number_variation.seg.txt',
    output:
        out_txt = an_path+'02_results_raw/CNV_segment/{sample}_CNV_segment_out.txt',
    threads: 2
    script:
        'scripts/CNV_segment_analysis.py'


rule Pindel:
    input:
        in_file = an_path+'01_sample_data/Pindel/{sample}.sanger_raw_pindel.raw_somatic_mutation.vcf.gz',
    output:
        out_txt = an_path+'02_results_raw/Pindel/{sample}_Pindel_out.txt',
    threads: 2
    script:
        'scripts/Pindel_analysis.py'


rule Splicing:
    input:
        in_file = an_path+'01_sample_data/Splicing/{sample}.star_splice_junctions.tsv.gz',
    output:
        out_txt = an_path+'02_results_raw/Splicing/{sample}_Splicing_out.txt',
    threads: 2
    script:
        'scripts/Splicing_analysis.py'


rule STAR_counts:
    input:
        in_file = an_path+'01_sample_data/STAR_counts/{sample}.augmented_star_gene_counts.tsv',
    output:
        out_txt = an_path+'02_results_raw/STAR_counts/{sample}_STAR_counts_out.txt',
    threads: 2
    script:
        'scripts/STAR_counts_analysis.py'
