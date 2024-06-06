# Handling TCGA data
## 1. Overall goals
- Work with TCGA data (e.g. of a specific study or primary site)
- Downoad data (or further data after a while)
- Download data matching specific cases from previous analyses

## 2. Your input
### 2.1. Filter and select TCGA data
- Go to: [https://portal.gdc.cancer.gov/analysis_page](https://portal.gdc.cancer.gov/analysis_page)
- Click on "Repository"

&emsp;&emsp;&ensp;<img src="figures/tcga_gdc_repository_files.png" style="width:1000px; position: relative; left: 40px">

- Filter by files:
    - Experimental strategy
    - Data type
    - Data format
    - Access
    - ...

&emsp;&emsp;&ensp;<img src="figures/tcga_gdc_files_filter.png" style="width:300px; position: relative; left: 40px">

- Filter by cases:
    - Primary site
    - Disease Type
    - Project
    - Gender
    - ...

&emsp;&emsp;&ensp;<img src="figures/tcga_gdc_cases_filter.png" style="width:1000px; position: relative; left: 40px">

- Select your files of interest or put all relevant data into the cart

&emsp;&emsp;&ensp;<img src="figures/tcga_gdc_select_file_cart.png" style="width:800px; position: relative; left: 40px">

- Click on the Cart symbol at the top to view the contents of the cart

&emsp;&emsp;&ensp;<img src="figures/tcga_gdc_cart_overview.png" style="width:1000px; position: relative; left: 40px">

- Click on "Download Associated Data". Download the following files:
    - Sample Sheet
    - Metadata
    - Optional: Clinical: TSV data for additional information

&emsp;&emsp;&ensp;<img src="figures/tcga_gdc_cart_ass_data.png" style="width:180px; position: relative; left: 40px">

- Click on "Download Cart". Download the Manifest for a data download with the gdc-client tool

&emsp;&emsp;&ensp;<img src="figures/tcga_gdc_cart_manifest.png" style="width:120px; position: relative; left: 40px">

- Follow these steps every time for your new analyses, also when you have new aspects or file types to consider later on


### 2.2. File locations for manifests and sample sheets
- Create a local analysis folder and in this folder, a "sample_sheets" folder with sub-folders "clinical_data", "manifests", and "sample_sheets_prior"

&emsp;&emsp;&ensp;<analysis_path>

&emsp;&emsp;&ensp;└── sample_sheets
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;├── clinical_data
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;├── manifests
<br>
&emsp;&emsp;&emsp;&emsp;&ensp;└── sample_sheets_prior    


### 2.3. Optional: Download TCGA access token for restricted access files
- For restricted access files:
    - Login at TCGA (with NIH account) for restricted access files
    - Download access token, save as secured file


### 2.3. Adapt the configuration file
- The configuration file ["data/config.yaml"](data/config.yaml) encompasses 


## 3. Pipeline steps explained
### 3.1. Combine manifest with sample sheet, filter for relevant files
- This is how a manifest file looks like:

&emsp;&emsp;&ensp;<img src="figures/tcga_manifest_file_example.png" style="width:1000px; position: relative; left: 40px">

- This is how a sample sheet looks like:

&emsp;&emsp;&ensp;<img src="figures/tcga_sample_sheet_example.png" style="width:1000px; position: relative; left: 40px">

- Merge manifest and sample sheet
- If previous selection of case IDs wanted, filter for specific case IDs of previous analysis
- create adapted filtered manifest file for gdc-client download

### 3.2. Download TCGA data via a manifest document and the GDC-client tool
- The script creates a new conda environment called "gdc_client" and downloads the gdc-client tool. If you have already installed the gdc-client in a separete conda environment, you can specify that in the configuration file.
- The script downloads the TCGA data from manifest specified in previous steps and/or the configuration file via the gdc-client
- The files from the manifest are downloaded into the following folder: analysis_path + '00_raw_data'

- Optional: for restricted access files:
    - Login at NIH for restricted access files
    - Download access token, save as secured file
- Download gdc-client tool

```
gdc-client download -m manifest.txt -t user-token.txt
```

### 3.3. Rename the downloaded files as case_id.file_suffix
- in manifest only id, filename with 36 different characters
- take merged manifest and sample sheet
- rename downloaded files and put them in new folders for each analysis

### 3.4. Analyze files
- with Snakemake pipeline


