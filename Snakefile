import pandas as pd
import os
import yaml

with open('data/config.yaml', 'r') as streamfile:
    conf_yaml = yaml.load(streamfile, Loader=yaml.FullLoader)

analysis_path = conf_yaml['analysis_path']
snakemake_caseids = conf_yaml['snakemake_caseids']

METHODS_LIST = conf_yaml['snakemake_methods']

def Prepare_Samples_Lists():
    if snakemake_caseids == True:
        manifest_list = conf_yaml['manifest_for_download']
        manifests = [pd.read_table(analysis_path+'sample_sheets/manifests_merged_sample_sheet/'+i.rsplit('_ready.',1)[0]+'_merged_sample_sheet.txt') 
                    for i in manifest_list]
        sample_sheet = pd.concat(manifests)
    else:
        sample_sheet = pd.read_table(snakemake_caseids)
    
    # check if files exist
    sample_sheet['is_file_caseid'] = sample_sheet['Path_sample'].apply(os.path.isfile)

    SAMPLES_LIST = []
    for method in METHODS_LIST:
        sample_sheet_met = sample_sheet[(sample_sheet['Folder']==method)&(sample_sheet['is_file_caseid']==True)].copy()
        SAMPLES_LIST.append(list(sample_sheet_met['Case ID'].unique()))
    
    return SAMPLES_LIST

SAMPLES_LIST = Prepare_Samples_Lists()

os.makedirs(analysis_path+'02_results_raw', exist_ok=True)
for m in METHODS_LIST:
    os.makedirs(analysis_path+'02_results_raw/'+m, exist_ok=True)
os.makedirs(analysis_path+'03_results_combined', exist_ok=True)

rule all:
    input:
        expand(analysis_path+'02_results_raw/'+METHODS_LIST[0]+'/{sample}_'+METHODS_LIST[0]+'_out.txt', sample=SAMPLES_LIST[0]),
        expand(analysis_path+'02_results_raw/'+METHODS_LIST[1]+'/{sample}_'+METHODS_LIST[1]+'_out.txt', sample=SAMPLES_LIST[1]),
        #expand(analysis_path+'02_results_raw/'+METHODS_LIST[2]+'/{sample}_'+METHODS_LIST[2]+'_out.txt', sample=SAMPLES_LIST[2]),
        # expand(analysis_path+'02_results_raw/BRASS/{sample}_{method}_out.txt', method=METHODS, sample=SAMPLES),


rule BRASS:
    input:
        in_file = analysis_path+'01_sample_data/BRASS/{sample}.wgs.BRASS.raw_structural_variation.vcf.gz',
        bedpe = analysis_path+'01_sample_data/BRASS/{sample}.wgs.BRASS.rerun_structural_variation.bedpe.gz',
    output:
        out_txt = analysis_path+'02_results_raw/BRASS/{sample}_BRASS_out.txt',
    threads: 2
    conda: 'envs/vcf_env.yaml'
    script:
        'scripts_snakemake/BRASS_analysis.py'


rule CaVEMan:
    input:
        in_file = analysis_path+'01_sample_data/CaVEMan/{sample}.wgs.CaVEMan.raw_somatic_mutation.vcf.gz',
    output:
        out_txt = analysis_path+'02_results_raw/CaVEMan/{sample}_CaVEMan_out.txt',
    threads: 2
    script:
        'scripts_snakemake/CaVEMan_analysis.py'


rule CNV_segment:
    input:
        in_file = analysis_path+'01_sample_data/CNV_segment/{sample}.wgs.ASCAT.copy_number_variation.seg.txt',
    output:
        out_txt = analysis_path+'02_results_raw/CNV_segment/{sample}_CNV_segment_out.txt',
    threads: 2
    script:
        'scripts_snakemake/CNV_segment_analysis.py'


rule Pindel:
    input:
        in_file = analysis_path+'01_sample_data/Pindel/{sample}.wgs.sanger_raw_pindel.raw_somatic_mutation.vcf.gz',
    output:
        out_txt = analysis_path+'02_results_raw/Pindel/{sample}_Pindel_out.txt',
    threads: 2
    script:
        'scripts_snakemake/Pindel_analysis.py'


rule Splicing:
    input:
        in_file = analysis_path+'01_sample_data/Splicing/{sample}.rna_seq.star_splice_junctions.tsv.gz',
    output:
        out_txt = analysis_path+'02_results_raw/Splicing/{sample}_Splicing_out.txt',
    threads: 2
    script:
        'scripts_snakemake/Splicing_analysis.py'


rule STAR_counts:
    input:
        in_file = analysis_path+'01_sample_data/STAR_counts/{sample}.rna_seq.augmented_star_gene_counts.tsv',
    output:
        out_txt = analysis_path+'02_results_raw/STAR_counts/{sample}_STAR_counts_out.txt',
    threads: 2
    script:
        'scripts_snakemake/STAR_counts_analysis.py'
