import yaml
import os
import subprocess
import shlex
import time
import sys

with open(sys.argv[1], 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

analysis_path = config_file['analysis_path']

log_messages = []

##### Check validity of important configuration file entries #####
msg0 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'##### Check validity of important configuration file entries #####')
print(msg0)
log_messages.append(msg0)

with open('data/config.yaml', 'r') as o:
    conf_lines = [l.rstrip() for l in o.readlines()]

with open('data/config.yaml', 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

# Check analysis path
analysis_path = config_file['analysis_path']

# Count errors
error_count = 0

# if analysis_path is valid
if not os.path.exists(analysis_path):
    error_count += 1
    msg1 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
            'Please make sure the analysis path exists and includes the manifest files. Your current analysis path is the following:' + analysis_path)
    print(msg1)
    log_messages.append(msg1)
else:
    # if analysis_path has "/" at the end
    if analysis_path[-1] != '/':
        analysis_path = analysis_path+'/'
        conf_lines2 = [l+'/\n' if l.startswith('analysis_path') else l+'\n' for l in conf_lines]
        with open('data/config.yaml', 'w') as w:
            w.writelines(conf_lines2)

    # if sample sheet and manifests folders exist
    if not os.path.exists(analysis_path+'sample_sheets/manifests') or not os.path.exists(analysis_path+'sample_sheets/sample_sheets_prior'):
        error_count += 1
        msg2 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                'Please make sure the paths sample_sheets/manifests and sample_sheets/sample_sheets_prior exist in your analysis folder.')
        print(msg2)
        log_messages.append(msg2)

    else:
        # if files in manifest and sample sheet folders exist
        manifests_prior = config_file['manifests_prior']
        sample_sheets_prior = config_file['sample_sheets_prior']

        for category_list, category_name in zip([manifests_prior, sample_sheets_prior], ['manifests', 'sample_sheets_prior']):
            for check_file in category_list:
                msg3 = None
                if not os.path.isfile(analysis_path+f'sample_sheets/{category_name}/'+check_file):
                    error_count += 1
                    msg3 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+f'{check_file} ist not existent. Please check again.')
                    print(msg3)
                    log_messages.append(msg3)
                else:
                    continue
        
        # if sample sheet for filtering is existent
        sample_sheet_filtering = config_file['sample_sheet_filtering']

        if sample_sheet_filtering != False:
            if not os.path.isfile(analysis_path+f'sample_sheets/sample_sheets_prior/'+sample_sheet_filtering):
                error_count += 1
                msg4 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                        f'{sample_sheet_filtering} is not existent. Please check again or input: False')
                print(msg4)
                log_messages.append(msg4)
        
        # if manifest for download is existent
        manifest_for_download = config_file['manual_manifest_download']

        if manifest_for_download != False:
            for check_file in manifest_for_download:
                msg5 = None
                if not os.path.isfile(analysis_path+f'sample_sheets/manifests_for_download/'+check_file):
                    error_count += 1
                    msg5 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+
                            f'{check_file} is not existent. Please check again or input: False')
                    print(msg5)
                    log_messages.append(msg5)

# Check if TCGA user token file is existent
tcga_user_token_file = config_file['tcga_user_token_file']
if tcga_user_token_file != False:
    if not os.path.isfile(tcga_user_token_file):
        error_count += 1
        msg6 = ('User token file not existent, please check again or if you do not want to use one, input: False.' +
                f' Your current analysis path is the following: {tcga_user_token_file}')
        print(msg6)
        log_messages.append(msg6)

# Error check
if error_count == 0:
    msg7 = time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'Your configuration file seems to be correct.'
    print(msg7)
    log_messages.append(msg7)
    with open(sys.argv[3], 'w') as o:
        for msg in log_messages:
            o.write(msg+'\n')
    with open(sys.argv[2], 'w') as o2:
        o2.write((time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+'Your config file seems correct.\n'))
else:
    msg7 = (time.strftime('%Y-%m-%d %H:%M:%S: ', time.localtime())+f'Your configuration file seems to have {error_count} errors. '+
            'If you are using the Snakemake pipeline, the pipeline will be interrupted in a few seconds.')
    print(msg7)
    log_messages.append(msg7)
    with open(sys.argv[3], 'w') as o:
        for msg in log_messages:
            o.write(msg+'\n')
