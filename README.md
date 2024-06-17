# Handling TCGA data

# --- Still under construction ---

## 1. Overall goals
- Work with TCGA data (e.g. of a specific study or primary site)
- Downoad data (or further data after a while)
- Download data matching specific cases from previous analyses

## 2. Your input
### 2.1. Prerequisites
- To use this script as a one-touch pipeline, you have to [create a conda environment](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) with the [Jupyter Notebook](https://anaconda.org/anaconda/jupyter) and the [pandas](https://anaconda.org/anaconda/pandas) package
- If you want to do the further analysis step with the [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline, you will have to create a conda environment with the Snakemake and pandas package. You can also use the [given Snakemake environment file](envs/snakemake_env.txt) with this command to create a conda environment:

```
conda create --name Snakemake --file envs/snakemake_env.txt
```

### 2.2. Filter and select TCGA data
- Go to: [https://portal.gdc.cancer.gov/analysis_page](https://portal.gdc.cancer.gov/analysis_page)
- Click on "Repository" ([image](figures/tcga_gdc_repository_files.png))

- Filter by files ([image](figures/tcga_gdc_files_filter.png)) or by cases.

- Select your files of interest or put all (relevant) data into the cart ([image](figures/tcga_gdc_select_file_cart.png)).

- Click on the Cart symbol at the top to view the contents of the cart ([image](figures/tcga_gdc_cart_overview.png)).

- Click on "Download Associated Data" ([image](figures/tcga_gdc_cart_ass_data.png)). Download the following files:
    - Sample Sheet
    - Metadata
    - Optional: Clinical: TSV data for additional information

- Click on "Download Cart" ([image](figures/tcga_gdc_cart_manifest.png)). Download the Manifest for a data download with the gdc-client tool.

- Follow these steps every time for your new analyses, also when you have new aspects or file types to consider later on.

### 2.3. File locations for manifests and sample sheets
- Create a local analysis folder and in this folder, a "sample_sheets" folder with sub-folders "clinical_data", "manifests", and "sample_sheets_prior"

&emsp;&emsp;&ensp;<analysis_path>

&emsp;&emsp;&ensp;└── sample_sheets
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;├── clinical_data
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;├── manifests
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;└── sample_sheets_prior    

### 2.4. Optional: Download TCGA access token for restricted access files
- For restricted access files, login at TCGA (with NIH account) and download an access token, save as a secured file.

### 2.5. Adapt the configuration file
- The configuration file ["data/config.yaml"](data/config.yaml) encompasses all necessary information for a run of the pipeline to download TCGA data from a given manifest file and sample sheet.
- Information on how to fill out the configuration file are prepared within the [configuration file](data/config.yaml).

### 2.6. Adapt the Snakemake pipeline
- This step is optional, as this Snakemake pipeline is only a template and is not ready to use for the analysis of your TCGA samples yet.
- If you have decided what to analyze, you can define the rules in the Snakemake pipeline.

### 2.7. Start the pipeline
- If you have done all previous steps, you can start the pipeline.
- If you feel comfortable with using the command line only, you can run the script ["TCGA_pipeline_python.py"](TCGA_pipeline_python.py) in an environment with working Python, [pandas](https://anaconda.org/anaconda/pandas), and [pyyaml](https://anaconda.org/conda-forge/pyyaml) packages.
- Otherwise, you can activate your conda environment with the installed Jupyter Notebook package.
```
conda activate <name_of_environment>
```
- Then, open the Jupyter Notebook ["TCGA_steps_code.ipynb"](TCGA_steps_code.ipynb) and either run the script cell by cell or everything at once.

## 3. Pipeline steps explained
### 3.1. Check validity of configuration file entries
- Checks whether all files listed in the configuration file are existent.

### 3.1. Combine manifest with sample sheet, filter for relevant files
- This is how a manifest file looks like:

<img src="figures/tcga_manifest_file_example.png" style="width:1000px; position: relative; left: 40px">

- This is how a sample sheet looks like:

<img src="figures/tcga_sample_sheet_example.png" style="width:1000px; position: relative; left: 40px">

- This script merges the manifest(s) and the sample sheet.
- If previous selection of case IDs is wanted, the script filters for specific case IDs of a previous analysis.
- Creates adapted filtered manifest file for gdc-client download.

### 3.2. Download TCGA data via a manifest document and the GDC-client tool
- The script creates a new conda environment called "gdc_client" and downloads the gdc-client tool. If you have already installed the gdc-client in a separete conda environment, you can specify that in the configuration file.
- The script downloads the TCGA data from manifest(s) specified in previous steps and/or the configuration file via the gdc-client.
- The files from the manifest are downloaded into the following folder: analysis_path + '00_raw_data'
- Download TCGA data via the gdc-client tool.
- Optional: for restricted access files, you need an access token.

```
gdc-client download -m manifest.txt (-t user-token.txt)
```

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
- 


