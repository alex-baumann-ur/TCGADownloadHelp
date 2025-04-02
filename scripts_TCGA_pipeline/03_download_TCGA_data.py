import yaml
import os
import subprocess
import shlex
import sys
import time

log_messages = []
error_count = 0

##### Download TCGA data via a manifest document and the GDC-client tool #####
print('##### Download TCGA data via a manifest document and the GDC-client tool #####')

msg0 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'##### Download TCGA data via a manifest document and the GDC-client tool #####')
print(msg0)
log_messages.append(msg0)

# Read in manifest download files (either from manual input in config.yaml or from previous pipeline steps)
def Create_Manifest_Download_List():
    with open(sys.argv[1], 'r') as streamfile:
        config_file = yaml.load(streamfile, Loader=yaml.FullLoader)
    
    manifest_download_list = config_file['manifest_for_download']
    analysis_path = config_file['analysis_path']

    if manifest_download_list == False:
        manifest_for_download = False
        error_count += 1
        msg1 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'Please execute the previous part of the pipeline (again) or input your manifest files for download manually in the data/config.yaml file.')
        print(msg1)
        log_messages.append(msg1)
    else:
        manifest_for_download = [analysis_path+'sample_sheets/manifests_for_download/'+i for i in manifest_download_list]
        msg2 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'Your manifest files for download are being processed further.')
        print(msg2)
        log_messages.append(msg2)
    
    return manifest_for_download


# Download the gdc-client in a new conda environment and run the gdc-client, if accepted
def Download_gdc_client():
    with open(sys.argv[1], 'r') as streamfile:
        config_file = yaml.load(streamfile, Loader=yaml.FullLoader)
    
    conda_gdc = config_file['conda_gdc']
    name_conda_gdc_env = False

    if conda_gdc == False:
        error_count += 1
        msg3 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'Please execute the TCGA data download in your own environment or '+
                'set "conda_gdc" in the data/config.yaml file to "First_install" to create a conda environment with the gdc-client')
        print(msg3)
        log_messages.append(msg3)
    elif conda_gdc == True:
        name_conda_gdc_env = 'gdc_client'
        msg4 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'The previously installed conda environment "gdc-client" will be used.')
        print(msg4)
        log_messages.append(msg4)
    elif conda_gdc == 'First_install':

        gdc_client_conda_cmd = shlex.split(f'conda create --name gdc_client --file envs/gdc_client.txt')
        subprocess.run(gdc_client_conda_cmd)
        
        conda_gdc_status_change_cmd = f"sed -i 's/^conda_gdc: .*/conda_gdc: True/' data/config.yaml"
        subprocess.run(shlex.split(conda_gdc_status_change_cmd))

        msg5 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'The conda environment "gdc-client" was created.')
        print(msg5)
        log_messages.append(msg5)

        name_conda_gdc_env = 'gdc_client'
    else:
        name_conda_gdc_env = conda_gdc
    
    return name_conda_gdc_env


# Prepare commands to download TCGA data, dependent on available user token file
def TCGA_Data_Download():
    with open(sys.argv[1], 'r') as streamfile:
        config_file = yaml.load(streamfile, Loader=yaml.FullLoader)
    
    tcga_user_token_file = config_file['tcga_user_token_file']
    analysis_path = config_file['analysis_path']
    raw_data_path = analysis_path + '00_raw_data'

    os.makedirs(raw_data_path, exist_ok=True)

    manifest_for_download = Create_Manifest_Download_List()

    name_conda_gdc_env = Download_gdc_client()

    if manifest_for_download == False:
        error_count += 1
        msg6 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'No TCGA data are download due to no manifest files.')
        print(msg6)
        log_messages.append(msg6)
    elif name_conda_gdc_env == False:
        error_count += 1
        msg7 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'No TCGA data are downloaded due to the specifications in the gdc client conda environment.')
        print(msg7)
        log_messages.append(msg7)
    else:
        msg8 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'TCGA data will be downloaded.')
        print(msg8)
        log_messages.append(msg8)
        for i,manifest_file in enumerate(manifest_for_download):
            name_log_file = sys.argv[2]+str(i)+'.txt'
            if tcga_user_token_file == False:
                msg9 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                        f'Download TCGA data with TCGA manifest {manifest_file.split("/")[-1]} without TCGA user token. '+
                        'Please check the log file')
                print(msg9)
                log_messages.append(msg9)
                command_download_tcga_data = f'conda run -n {name_conda_gdc_env} gdc-client download -m {manifest_file} --log-file {name_log_file}'
            else:
                msg10 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                        f'Download TCGA data with TCGA manifest {manifest_file.split("/")[-1]} with TCGA user token file {tcga_user_token_file.split("/")[-1]}')
                print(msg10)
                log_messages.append(msg10)
                command_download_tcga_data = f'conda run -n {name_conda_gdc_env} gdc-client download -m {manifest_file} -t {tcga_user_token_file} --log-file {name_log_file}'
            process = subprocess.Popen(command_download_tcga_data, cwd=raw_data_path, shell=True)
            process.wait()
        
        msg11 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                f'The TCGA data download step was finished. Please check the log files ({sys.argv[2]}<i>.txt) for a successful download process.')
        print(msg11)
        log_messages.append(msg11)

        with open(sys.argv[3], 'w') as o3:
            o3.write(msg11)

    with open(sys.argv[4], 'w') as o2:
        for msg in log_messages:
            o2.write(msg+'\n')

TCGA_Data_Download()