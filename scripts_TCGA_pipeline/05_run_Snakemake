import subprocess
import yaml

with open('data/config.yaml', 'r') as streamfile:
    config_file = yaml.load(streamfile, Loader=yaml.FullLoader)

conda_snakemake = config_file['conda_snakemake']
snakemake_threads = config_file['snakemake_threads']

# run Snakemake pipeline
command_snakemake = f'conda run -n {conda_snakemake} snakemake --cores {snakemake_threads} --use-conda -k'
process_snakemake = subprocess.Popen(command_snakemake, shell=True)
process_snakemake.wait()
