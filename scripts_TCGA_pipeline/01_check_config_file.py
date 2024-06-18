import yaml
import os
import subprocess
import shlex

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