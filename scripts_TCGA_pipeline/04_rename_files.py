import os
import pandas as pd
import yaml
from datetime import datetime

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