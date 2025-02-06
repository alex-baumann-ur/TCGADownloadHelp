import yaml
import os
import subprocess
import shlex
import pandas as pd
import sys
import time
from exceptions_renaming import file_rename_exceptions

log_messages = []

##### Combine manifest with sample sheet, filter for relevant files to download #####
msg0 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'##### Combine manifest with sample sheet, filter for relevant files to download #####')
print(msg0)
log_messages.append(msg0)

# load config file
with open(sys.argv[1], 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

# TCGA study name
tcga_study = config_file['tcga_study']

# sample type(s) to download
sample_type = config_file['sample_type']

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
msg1 = None
for manifest in manifests_dfs:
    df = pd.merge(manifest, sample_sheets_merge, how='left', left_on=['id','filename'], right_on=['File ID','File Name'], indicator=True)
    
    if len(df[df['_merge']=='left_only'])>0:
        msg1 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'It seems like a manifest entry cannot be mapped to the sample sheet. Please check your manifest file and sample sheet again.'+
                'Only mapped samples are further processed.')
        print(msg1)
        log_messages.append(msg1)
        df = df[df['_merge']=='both'].copy()
    else:
        msg2 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'Manifest file has been mapped to the sample sheet.')
        print(msg2)
        log_messages.append(msg2)
    
    merge_manifest_sample_sheet.append(df.drop(columns='_merge'))
    
    
# Check whether files in manifest have already been downloaded and exclude them from the manifest and prepare renaming
def Check_Already_Downloaded(manifest_sample_sheet_df, method_dict, tcga_study=tcga_study):
    df = manifest_sample_sheet_df.copy().drop_duplicates()

    df['Folder'] = ''
    for met in method_dict.keys():
        df.loc[df['File Name'].str.lower().str.contains(met.lower()), 'Folder'] = method_dict[met]
    
    # exceptions for renaming the file (read in from 'exceptions_renaming.py')
    df = file_rename_exceptions(df, tcga_study)

    # build new names for the files
    df['Path_raw'] = analysis_path+'00_raw_data/'+df['File ID']+'/'+df['File Name']
    df['Path_sample'] = analysis_path+'01_sample_data/'+df['Folder']+'/'+df['Case ID']+df['File Suffix']

    # check if file was already downloaded
    df['is_file_raw'] = df['Path_raw'].apply(os.path.isfile)
    df['is_file_caseid'] = df['Path_sample'].apply(os.path.isfile)

    msg3 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'Checked whether files have already been downloaded.')
    print(msg3)
    log_messages.append(msg3)
    
    return df


# if filtering of manifest file on case IDs of previous sample sheet is wanted
# include step of checking whether files are already downloaded
sample_sheet_filtering = config_file['sample_sheet_filtering']
method_dict = config_file['methods_dict']

if sample_sheet_filtering != False:
    sample_sheet_filter = pd.read_table(analysis_path+'sample_sheets/sample_sheets_prior/'+sample_sheet_filtering)
    case_ids = list(sample_sheet_filter['Case ID'].unique())
    list_manifests_for_download = []
    for manifest, file_name_manifest in zip(merge_manifest_sample_sheet, manifests):
        # only include primary tumor samples
        filtered_manifest = manifest[(manifest['Case ID'].isin(case_ids))&(manifest['Sample Type'].isin(sample_type))].copy()

        msg4 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'Manifest has been filtered for provided Case IDs.')
        print(msg4)
        log_messages.append(msg4)
        
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

        msg5 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'Manifest has been filtered for already downloaded files and new manifest file prepared for download.')
        print(msg5)
        log_messages.append(msg5)
            
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

        msg6 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'Manifest has been filtered for already downloaded files and new manifest file prepared for download.')
        print(msg6)
        log_messages.append(msg6)
        

# manifests for download
# read in configuration file again because it was changed the prior step
with open(sys.argv[1], 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

# filtered manifest(s) in there, but only names in 'sample_sheets/manifests/
if config_file['manifest_for_download'] == False:
    manifests_pipeline_files = [m.split('/')[-1] for m in list_manifests_for_download]
    manifests_pipeline = ', '.join(manifests_pipeline_files) # only file names
    
    manifest_file_change_cmd = f"sed -i 's/^manifest_for_download: .*/manifest_for_download: [{manifests_pipeline}]/' data/config.yaml"
    subprocess.run(shlex.split(manifest_file_change_cmd))

    msg7 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'Config file has been updated with new manifest files to download.')
    print(msg7)
    log_messages.append(msg7)

with open(sys.argv[3], 'w') as o2:
    for msg in log_messages:
        o2.write(msg+'\n')

if msg1 == None:
    with open(sys.argv[2], 'w') as o3:
        o3.write(time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'Manifest files have been merged with sample sheet and been prepared for download.')
