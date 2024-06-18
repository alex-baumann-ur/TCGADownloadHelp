import yaml
import os
import subprocess
import shlex

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