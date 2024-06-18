import yaml
import os
import subprocess
import shlex
import pandas as pd
from datetime import datetime

##### Check validity of important configuration file entries #####
print('##### Check validity of important configuration file entries #####')

with open('data/config.yaml', 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

# Check analysis path
analysis_path = config_file['analysis_path']
# if analysis_path is valid
if not os.path.exists(analysis_path):
    print('Please make sure the analysis path exists and includes the manifest files. Your current analysis path is the following:')
    print(analysis_path)
else:
    # if analysis_path has "/" at the end
    if analysis_path[-1] != '/':
        analysis_path = analysis_path+'/'
        change_config_cmd = f"sed -i '/^analysis_path: .*/s/$/\//' data/config.yaml"
        subprocess.run(shlex.split(change_config_cmd))
    # if sample sheet and manifests folders exist
    if not os.path.exists(analysis_path+'sample_sheets/manifests') or not os.path.exists(analysis_path+'sample_sheets/sample_sheets_prior'):
        print('Please make sure the paths sample_sheets/manifests and sample_sheets/sample_sheets_prior exist in your analysis folder.')
    else:
        # if files in manifest and sample sheet folders exist
        manifests_prior = config_file['manifests_prior']
        sample_sheets_prior = config_file['sample_sheets_prior']

        for category_list, category_name in zip([manifests_prior, sample_sheets_prior], ['manifests', 'sample_sheets_prior']):
            for check_file in category_list:
                if not os.path.isfile(analysis_path+f'sample_sheets/{category_name}/'+check_file):
                    print(f'{check_file} ist not existent. Please check again.')
                else:
                    continue
        
        # if sample sheet for filtering is existent
        sample_sheet_filtering = config_file['sample_sheet_filtering']

        if sample_sheet_filtering != False:
            if not os.path.isfile(sample_sheet_filtering):
                print(f'{sample_sheet_filtering} is not existent. Please check again or input: False')
        
        # if manifest for download is existent
        manifest_for_download = config_file['manifest_for_download']

        if manifest_for_download != False:
            for check_file in manifest_for_download:
                if not os.path.isfile(analysis_path+f'sample_sheets/manifests/'+check_file):
                    print(f'{check_file} is not existent. Please check again or input: False')

# Check if TCGA user token file is existent
tcga_user_token_file = config_file['tcga_user_token_file']
if tcga_user_token_file != False:
    if not os.path.isfile(tcga_user_token_file):
        print(f'User token file not existent, please check again or if you do not want to use one, input: False. Your current analysis path is the following:')
        print(tcga_user_token_file)



##### Combine manifest with sample sheet, filter for relevant files to download #####
print('##### Combine manifest with sample sheet, filter for relevant files to download #####')

# load config file
with open('data/config.yaml', 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

# analysis path
analysis_path = config_file['analysis_path']

# original (prior) manifest files
manifests = config_file['manifests_prior']
manifests_dfs = [pd.read_table(analysis_path+'sample_sheets/manifests/'+i) for i in manifests]

# create folder for downloadable manifest files and merged with sample sheet
os.makedirs(analysis_path+'sample_sheets/manifests_for_download', exist_ok=True)
os.makedirs(analysis_path+'sample_sheets/manifests_merged_sample_sheet', exist_ok=True)

# original (prior) sample sheets, merge them
sample_sheets = [pd.read_table(analysis_path+'sample_sheets/sample_sheets_prior/'+i) for i in config_file['sample_sheets_prior']]
sample_sheets_merge = pd.concat(sample_sheets)

sample_sheets_merge['Case ID'] = sample_sheets_merge['Case ID'].str.split(', ', expand=True)[0]

# merge each manifest file with the sample sheet
merge_manifest_sample_sheet = []
for manifest in manifests_dfs:
    df = pd.merge(manifest, sample_sheets_merge, how='left', left_on=['id','filename'], right_on=['File ID','File Name'], indicator=True)

    
    if len(df[df['_merge']=='left_only'])>0:
        print('It seems like a manifest entry cannot be mapped to the sample sheet. Please check your manifest file and sample sheet again.')
        print('Only mapped samples are further processed.')
        df = df[df['_merge']=='both'].copy()
    merge_manifest_sample_sheet.append(df.drop(columns='_merge'))


# Check whether files in manifest have already been downloaded an exclude them from the manifest
def Check_Already_Downloaded(manifest_sample_sheet_df, method_dict):
    df = manifest_sample_sheet_df.copy().drop_duplicates()
    
    df['Folder'] = ''
    for met in method_dict.keys():
        df.loc[df['File Name'].str.contains(met), 'Folder'] = method_dict[met]
    
    # 36 characters file ID
    df['File Suffix'] = df['File Name'].str[36:]
    df['Path_raw'] = analysis_path+'00_raw_data/'+df['File ID']+'/'+df['File Name']
    df['Path_sample'] = analysis_path+'01_sample_data/'+df['Folder']+'/'+df['Case ID']+df['File Suffix']

    # check if file was already downloaded
    df['is_file_raw'] = df['Path_raw'].apply(os.path.isfile)
    df['is_file_caseid'] = df['Path_sample'].apply(os.path.isfile)
    
    return df


# if filtering of manifest file on case IDs of previous sample sheet is wanted
# include step of checking whether files are already downloaded
sample_sheet_filtering = config_file['sample_sheet_filtering']
method_dict = config_file['methods_dict']

if sample_sheet_filtering != False:
    sample_sheet_filter = pd.read_table(analysis_path+'sample_sheets/'+sample_sheet_filtering)
    case_ids = list(sample_sheet_filter['Case ID'].unique())
    list_manifests_for_download = []
    for manifest, file_name_manifest in zip(merge_manifest_sample_sheet, manifests):
        # only include primary tumor samples
        filtered_manifest = manifest[(manifest['Case ID'].isin(case_ids))&(manifest['Sample Type']=='Primary Tumor')].copy()
        
        # check for already downloaded files
        filtered_manifest = Check_Already_Downloaded(filtered_manifest, method_dict)

        # save filtered and merged manifest with old and new file names
        filtered_manifest.to_csv(analysis_path+'sample_sheets/manifests_merged_sample_sheet/'+file_name_manifest.rsplit('.', 1)[0]+'_merged_sample_sheet.txt', 
                                 sep='\t', index=False)

        # only include not already downloaded entries
        filtered_checked_manifest = filtered_manifest[(filtered_manifest['is_file_raw']==False)&(filtered_manifest['is_file_caseid']==False)]
        filtered_manifest_name = analysis_path+'sample_sheets/manifests_for_download/'+file_name_manifest.rsplit('.', 1)[0]+'_ready.txt'
        list_manifests_for_download.append(filtered_manifest_name)
        # save filtered manifest for download
        filtered_checked_manifest[['id','filename','md5','size','state']].to_csv(filtered_manifest_name, sep='\t', index=False)
            
else:
    list_manifests_for_download = []
    for manifest, file_name_manifest in zip(merge_manifest_sample_sheet, manifests):
        # check for already downloaded files
        manifest = Check_Already_Downloaded(manifest, method_dict)

        # save merged manifest with old and new file names
        manifest.to_csv(analysis_path+'sample_sheets/manifests_merged_sample_sheet/'+file_name_manifest.rsplit('.', 1)[0]+'_merged_sample_sheet.txt', sep='\t', index=False)

        # only include not already downloaded entries
        checked_manifest = manifest[(manifest['is_file_raw']==False)&(manifest['is_file_caseid']==False)]
        manifest_name = analysis_path+'sample_sheets/manifests_for_download/'+file_name_manifest.rsplit('.', 1)[0]+'_ready.txt'
        list_manifests_for_download.append(manifest_name)
        # save manifest for download
        checked_manifest[['id','filename','md5','size','state']].to_csv(manifest_name, sep='\t', index=False)
        

# manifests for download
# read in configuration file again because it was changed the prior step
with open('data/config.yaml', 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

# filtered manifest(s) in there, but only names in 'sample_sheets/manifests/
if config_file['manifest_for_download'] == False:
    manifests_pipeline_files = [m.split('/')[-1] for m in list_manifests_for_download]
    manifests_pipeline = ', '.join(manifests_pipeline_files) # only file names
    
    manifest_file_change_cmd = f"sed -i 's/^manifest_for_download: .*/manifest_for_download: [{manifests_pipeline}]/' data/config.yaml"
    subprocess.run(shlex.split(manifest_file_change_cmd))



##### Download TCGA data via a manifest document and the GDC-client tool #####
print('##### Download TCGA data via a manifest document and the GDC-client tool #####')

# Read in manifest download files (either from manual input in config.yaml or from previous pipeline steps)
def Create_Manifest_Download_List():
    with open('data/config.yaml', 'r') as streamfile:
        config_file = yaml.load(streamfile, Loader=yaml.FullLoader)
    
    manifest_download_list = config_file['manifest_for_download']
    analysis_path = config_file['analysis_path']

    if manifest_download_list == False:
        manifest_for_download = False
        # merge_manifest_sample_sheet = False
        print('Please execute the previous part of the pipeline or input your files manually in the data/config.yaml file')
    else:
        manifest_for_download = [analysis_path+'sample_sheets/manifests_for_download/'+i for i in manifest_download_list]
        # merge_manifest_sample_sheet = []
    
    return manifest_for_download


# Download the gdc-client in a new conda environment and run the gdc-client, if accepted
def Download_gdc_client():
    with open('data/config.yaml', 'r') as streamfile:
        config_file = yaml.load(streamfile, Loader=yaml.FullLoader)
    
    conda_gdc = config_file['conda_gdc']
    name_conda_gdc_env = False

    if conda_gdc == False:
        print('Please execute the TCGA data download in your own environment or '+
            'set "conda_gdc" in the data/config.yaml file to "First_install" to create a conda environment with the gdc-client')
    elif conda_gdc == True:
        name_conda_gdc_env = 'gdc_client'
    elif conda_gdc == 'First_install':
        gdc_client_conda_cmd = shlex.split(f'conda create --name gdc_client --file envs/gdc_client.txt')
        subprocess.run(gdc_client_conda_cmd)
        
        conda_gdc_status_change_cmd = f"sed -i 's/^conda_gdc: .*/conda_gdc: True/' data/config.yaml"
        subprocess.run(shlex.split(conda_gdc_status_change_cmd))

        name_conda_gdc_env = 'gdc_client'
    else:
        name_conda_gdc_env = conda_gdc
    
    return name_conda_gdc_env


# Prepare commands to download TCGA data, dependent on available user token file
def TCGA_Data_Download():
    with open('data/config.yaml', 'r') as streamfile:
        config_file = yaml.load(streamfile, Loader=yaml.FullLoader)
    
    tcga_user_token_file = config_file['tcga_user_token_file']
    analysis_path = config_file['analysis_path']
    raw_data_path = analysis_path + '00_raw_data'

    os.makedirs(raw_data_path, exist_ok=True)

    manifest_for_download = Create_Manifest_Download_List()

    name_conda_gdc_env = Download_gdc_client()

    if manifest_for_download == False:
        print('No TCGA data are download due to no manifest files.')
    elif name_conda_gdc_env == False:
        print('No TCGA data are downloaded due to the specifications in the gdc client conda environment.')
    else:
        for manifest_file in manifest_for_download:
            if tcga_user_token_file == False:
                print(f'Download TCGA data with TCGA manifest {manifest_file.split("/")[-1]} without TCGA user token')
                command_download_tcga_data = f'conda run -n {name_conda_gdc_env} gdc-client download -m {manifest_file}'
            else:
                print(f'Download TCGA data with TCGA manifest {manifest_file.split("/")[-1]} with TCGA user token file {tcga_user_token_file.split("/")[-1]}')
                command_download_tcga_data = f'conda run -n {name_conda_gdc_env} gdc-client download -m {manifest_file} -t {tcga_user_token_file}'
            process = subprocess.Popen(command_download_tcga_data, cwd=raw_data_path, shell=True)
            process.wait()

TCGA_Data_Download()


##### Rename file names with symlinks #####
print('##### Rename file names with symlinks #####')

with open('data/config.yaml', 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

analysis_path = config_file['analysis_path']
os.makedirs(analysis_path+'01_sample_data', exist_ok=True)

### work with already existing merged files in folder 'manifests_merged_sample_sheet'
downloaded_manifests = config_file['manifest_for_download']
downloaded_manifests = pd.concat([pd.read_table(analysis_path+'sample_sheets/manifests_merged_sample_sheet/'+i.rsplit('_ready.',1)[0]+'_merged_sample_sheet.txt') 
                                  for i in downloaded_manifests])

# create folders for methods
methods = list(downloaded_manifests['Folder'].unique())
for m in methods:
    os.makedirs(analysis_path+'01_sample_data/'+m, exist_ok=True)

# create symlinks for samples to sort them into method folders
list_rename_successful = []
for raw_path, sample_path in zip(list(downloaded_manifests['Path_raw']), list(downloaded_manifests['Path_sample'])):
    try:
        os.symlink(raw_path, sample_path)
        list_rename_successful.append('Success')
    except FileExistsError:
        list_rename_successful.append('Already exists')
    except FileNotFoundError:
        list_rename_successful.append('Not found')

downloaded_manifests['rename_success'] = list_rename_successful

downloaded_manifests[['Path_raw','Path_sample','rename_success']].to_csv((analysis_path+'01_sample_data/'+datetime.now().strftime('%Y-%m-%d_%H-%M-%S')+'_rename_files.txt'), 
                                                                         sep='\t', index=False)



##### Analyze files with Snakemake pipeline #####
print('##### Analyze files with Snakemake pipeline #####')

with open('data/config.yaml', 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

conda_snakemake = config_file['conda_snakemake']
snakemake_threads = config_file['snakemake_threads']

# run Snakemake pipeline
command_snakemake = f'conda run -n {conda_snakemake} snakemake --cores {snakemake_threads} --use-conda -k'
process_snakemake = subprocess.Popen(command_snakemake, shell=True)
process_snakemake.wait()

print('##### Pipeline was run successfully #####')
