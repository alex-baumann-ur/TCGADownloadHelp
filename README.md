# Handling TCGA data

## 1. Overall goals
- Work with TCGA data (e.g. of a specific study or primary site)
- Downoad data (or further data after a while)
- Download data matching specific cases from previous analyses

## 2. Your input
### 2.1. Prerequisites
- To use this TCGADownloadHelper with the Jupyter Notebook or the Snakemake pipeline, you have to [create a conda environment](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) with the required packages, which can be found in the environment yaml file [TCGADownloadHelper_env.yaml](envs/TCGADownloadHelper_env.yaml). Use this command to create a matching conda environment:

```
conda env create --name TCGAHelper -f envs/TCGADownloadHelper_env.yaml
```

### 2.2 Quick example
- To test out the pipeline and its functionality very quickly, there is a "test_example" folder including a test manifest file and sample sheet for 4 files of the TCGA database.
- After you downloaded the whole GitHub repository, create a local folder on your computer where the TCGA data should be downloaded.
- Copy the "sample_sheets" subfolder of the "test_example" folder to your analysis folder.
- Add the path of your analysis folder in the "config.yaml" file under "analysis_path:". If you are using the Jupyter Notebook, you have to add the analysis path in the first line as well.
- Activate your TCGAHelper conda environment and run the pipeline (further explained in 2.8).
- Important: From time to time, the GDC server is overloaded or slow, which might lead to aborted downloads. If that happens, please execute the whole pipeline again.

### 2.3. Filter and select TCGA data
- Go to: [https://portal.gdc.cancer.gov/analysis_page](https://portal.gdc.cancer.gov/analysis_page)
- Click on "Repository" ([image](figures/tcga_gdc_repository_files.png))

- Filter by files ([image](figures/tcga_gdc_files_filter.png)) or by cases ([image](figures/tcga_gdc_cases_filter.png)).

- Select your files of interest or put all (relevant) data into the cart ([image](figures/tcga_gdc_select_file_cart.png)).

- Click on the Cart symbol at the top to view the contents of the cart ([image](figures/tcga_gdc_cart_overview.png)).

- Click on "Download Associated Data" ([image](figures/tcga_gdc_cart_ass_data.png)). Download the Sample Shee (and optional files like Clinical data or Metadata for additional information).

- Click on "Download Cart" ([image](figures/tcga_gdc_cart_manifest.png)). Download the Manifest for a data download with the gdc-client tool.

- Follow these steps every time for your new analyses, also when you have new aspects or file types to consider later on.

### 2.4. File locations for manifests and sample sheets
- In your local analysis folder, a "sample_sheets" folder with subdirectories "clinical_data", "manifests", and "sample_sheets_prior" is created by the pipeline. Please save the previously downloaded data (especially Manifest file, Sample Sheet) in the according folders.

&emsp;&emsp;&ensp;<analysis_path>

&emsp;&emsp;&ensp;└── sample_sheets
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;├── clinical_data
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;├── manifests
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;└── sample_sheets_prior    

### 2.5. Optional: Download TCGA access token for restricted access files
- For restricted access files, login at TCGA (with NIH account) and download an access token, save as a secured file.

### 2.6. Adapt the configuration file
- The configuration file ["data/config.yaml"](data/config.yaml) encompasses all necessary information for a run of the pipeline to download TCGA data from a given manifest file and sample sheet.
- Information on how to fill out the configuration file are prepared within the [configuration file](data/config.yaml).

### 2.7. Adapt the Snakemake sample analysis pipeline
- This step is optional, as this Snakemake sample analysis pipeline is only a template and is not ready to use for the analysis of your TCGA samples yet.
- If you have decided on what to analyze, you can define the rules in the [Snakefile_sample_analysis](Snakefile_sample_analysis).
- Each rule requires a Python script with the analysis methods or an adapted shell command in the rule.
- The Python scripts for the analysis are located in the folder [scripts_snakemake](scripts_snakemake).

### 2.8. Start the pipeline
- If you have done all previous steps, you can start the TCGADownloadHelp pipeline after activating your TCGAHelper conda environment with the required packages.
```
conda activate TCGAHelper
```
1. If you want to use the Jupyter Notebook, open Jupyter in a browser by running:
```
jupyter notebook
```
- Go to the ["TCGA_steps_code.ipynb"](TCGA_steps_code.ipynb) file, put in your analysis_path and follow the steps for adapting the config.file.
- After that, run the script cell by cell.

2. If you prefer using Snakemake, you have to decide on the amount of cores to use. You can run the pipeline in the TCGAHelper environment using the default parameters with:
```
snakemake --cores <cores>
```

## 3. Pipeline steps explained
### 3.1. Check validity of configuration file entries
- Checks whether all files listed in the configuration file are existent.

### 3.1. Combine manifest with sample sheet, filter for relevant files
- This is how a manifest file looks like:

<img src="figures/tcga_manifest_file_example.png" style="width:1000px; position: relative; left: 40px">

- This is how a sample sheet looks like:

<img src="figures/tcga_sample_sheet_example.png" style="width:1000px; position: relative; left: 40px">

- This script merges the manifest(s) and the sample sheet.
- If previous selection of case IDs is wanted, the script filters for specific case IDs of a previous analysis. To use that option, you have to adapt the [config.yaml](data/config.yaml) file and input a sample sheet under "sample_sheet_filtering" with the sample IDs you want to use.
- Creates adapted filtered manifest file for gdc-client download.

### 3.2. Download TCGA data via a manifest document and the GDC-client tool
- The script downloads the TCGA data from manifest(s) specified in previous steps and/or the configuration file via the gdc-client tool, which is listed in the conda environment yaml file. The following command is used within the pipeline:
```
gdc-client download -m manifest.txt (-t user-token.txt)
```
- The files from the manifest are downloaded into the following folder: analysis_path + '00_raw_data'
- Optional: for restricted access files, you need an access token, which can be added in the [config.yaml](data/config.yaml) file under "tcga_user_token_file".


### 3.3. Rename the downloaded files as case_id.file_suffix
- The filename consists of a suffix of 36 different characters as a unique id and is saved in a separate folder with another unique id of 36 characters. This is an example:

&emsp;&emsp;&ensp;00_raw_data
<br>
&emsp;&emsp;&ensp;├── 45c6656d-a8f3-4968-aa40-142f3a340dde
<br>
&emsp;&emsp;&ensp;│&emsp;&ensp;└── ec765dd2-6541-4cdc-a26c-6d116398dc87.rna_seq.augmented_star_gene_counts.tsv
<br>
&emsp;&emsp;&ensp;└── 5d703777-b3db-4ed5-952b-203a5641767e
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;└── 1138323c-5a60-4535-bfb3-31703a106798.rna_seq.augmented_star_gene_counts.tsv

- The manifest and sample sheet are merged to map the case id to the unique file id.
this script changes the suffix to the case id.
- The downloaded files are renamed and sorted in new folders for each analysis. This is an example:

&emsp;&emsp;&ensp;01_sample_data
<br>
&emsp;&emsp;&ensp;└── STAR_counts
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;├── TCGA-18-3410.rna_seq.augmented_star_gene_counts.tsv
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;└── TCGA-33-4589.rna_seq.augmented_star_gene_counts.tsv

### 3.4. Analyze files
- A Snakemake pipeline can be used to analyze all downloaded data at once (if wanted).
- The Snakemake pipeline is a template and is not ready to use for your analysis.
- Please adapt the rules in the [Snakemake file](Snakefile) to use the Snakemake pipeline.
- After that, you have to set the "snakemake_sample_analysis" parameter of the [config.yaml](data/config.yaml) file to True.
- You can run the sample analysis pipeline in the Jupyter Notebook by just running the last cell.
- If you want to use Snakemake on the command line, you can run the following command:
```
snakemake -s Snakefile_sample_analysis --cores <cores>
```
