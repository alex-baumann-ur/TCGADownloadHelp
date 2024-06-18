import yaml
import os
import subprocess
import shlex
import pandas as pd

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