# StructMuMMI 
StructMuMMI is a framework organized in two blocks: the preprocessing and the models training.

## Requirements

- [conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)

## Installation
### 1. Clone the repository

```bash
git clone git@github.com:bca-95/structmummi.git
cd structmummi
```

### 2. Install dependency tools, references, conda environments

```
cd structmummi
bash setup.sh
```

### 3.  Run the scripts from anywhere

Run this from the folder that contains `structmummi/`:

```bash
echo "export PATH=\"$(pwd)/structmummi:\$PATH\"" >> ~/.bashrc
echo "export PATH=\"$(pwd)/structmummi/preprocessing:\$PATH\"" >> ~/.bashrc
echo "export PATH=\"$(pwd)/structmummi/training:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

*(optional)* Set the temporary file.
The pipeline writes temporary files to "$TMPDIR" (default: "/tmp").
On clusters, "/tmp" is often small: point "TMPDIR" to a larger disk, either for one run:
(Important to clean it regularly)
```
Add "export TMPDIR=/path/to/largerdiskfor/tmp" to your "~/.bashrc"
Then source, 
source ~/.bashrc
```

### Preprocessing
The preprocessing steps can be executed (1) fully at once or (2) step by steps for more controllable process, using the function```bash --start_from```.\
If executed step by step, files created by previous steps has to be already created as only the options below are modifiable, not files name.\
(1) 
```bash
pre_structmummi.sh --input_control /path/to/control --input_case /path/to/case --ref_hg_target hg38 <--perform_imputation no (default:yes)>
```
(2)
```bash
pre_structmummi.sh --input_control /path/to/control --input_case /path/to/case --ref_hg_target hg38 --start_from 2
```

| Option | Values | Default | Description |
|---|---|---|---|

Required:\
  |'--input_control'|'name_control'|'required'|  PLINK basename of the control dataset (without .bed/.bim/.fam)\
  |'--input_case'   |'name_case'   |'required'|  PLINK basename of the case dataset\
  |'--ref_hg_target'|'genome build ID'|'hg19' or 'hg38'|' Target genome build\
 
Options:\
  |'--perform_imputation'|'bool'|'yes' or 'no'| run phasing + imputation\
  **It is recommended to perform imputation in presence of non-genotype variants.
  Warning: Imputation takes days (~4 days for samples size of 4K)**\
  |'--eur_marginX'           |'fraction'|'0.3'| Margin around EUR samples on PC1, fraction of EUR range\
  |' --eur_marginY'          |'fraction'|'0.3'| Margin around EUR samples on PC2, fraction of EUR range \
  |'--cpu_threads_phasing'   |'number'  |'30' |Threads per chromosome for Eagle \
  |'--cpu_threads_imputation'|'number'  |'16' | Threads per chromosome for Minimac4 \
  |'--hc_threshold '         |'fraction'|'0.1'| Hard-call threshold after imputation \
  |'--maf'                   |'fraction'|'0.01'| Minor allele frequency threshold \
  |'--cleanup'               |'bool'    |'yes'| yes|no: delete intermediate files at the end \
  **<ins>cleanup</ins> delete unnecessary intermediate files, excluding phasing, imputation, images, the genotypes matrix,...**\
  |'--start_from'            |'number'  |'1'   | It executes steps from the indicated number to the end. Step to start from, 1-6.  \
  |1: Take the inputs and harmonize the datasets to 1000G reference panel\ 
  |2: perform stratification and extract EUR samples\
  |3: perform imputation \
  |4: perform quality control \
  |5: Prepare the genotype matrix and delete the intermediate files\
  |'--help'               | Description of the option \

### Models training
Giving a case-control genotype matrix, an ensemble of models is trained, SNPs contribution to the disease prediction are evaluated, and prioritized SNPs (a consensus across the ensemble of models) are listed and mapped to an interaction network.
For more controllable process, this training block can be performed as a whole or step by step, using the function ```bash --start_from```.

```bash
train_structmummi.sh <--start_from 1 (default: 1)>
```

| Option | Values | Default | Description |
Options:\
|'--start_from'            |'number'  |'1'   | It executes steps from the indicated number to the end. Step to start from, 1-2 \
| 1: Perform an ensemble of "default: 10" models \
| 2: Do the SNPs prioritization\
|'--n_repeats'            |'number'  |'10'   | Number of ensemble of models\
|'--interation_evals'    |'number'  |'200'  | Number of interation during tuning\


**An extensive option are describe and available by executing**

```bash
train_structmummi.sh --help <or -h>>
```

### Repository layout

```
structmummi/
├─- preprocessing/
│   ├── run_merge.sh                 # main script
│   ├── run_liftover.sh              # called by run_merge.sh
│   └── harmonization_to_1000g.sh    # allele harmonization to 1000G
├── external_tools/
│   ├── liftover/                    # UCSC chain files (included)
│   ├── plink/plink2                 # not included, see Installation
│   └── bcftools-1.23.1/             # not included, see Installation
├── gcta_env.yml
└── liftover_env.yml
```

## Data availability

The WTCCC genotype data are under controlled access and are **not** distributed with this repository.
They can be requested from the [European Genome-phenome Archive (EGA)](https://ega-archive.org/).
