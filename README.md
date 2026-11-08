# StructMuMMI 
StructMuMMI is a framework organized in two blocks: the preprocessing and the models training.


### Preprocessing
Use the following script: [WILL BE ADDED LATER]
OR for controllable process:

Preprocessing pipeline for W0TCCC genotype data (PLINK format):
### Models training
[WILL BE UPDATED LATER]

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
echo "export PATH=\"$(pwd)/structmummi/preprocessing:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

## Usage
run_merge.sh
| Option | Description | Default |\
|---|---|---|\
| `--input_control` | PLINK basename of the control dataset | required |\
| `--input_case` | PLINK basename of the case dataset | required |\
| `--ref_hg_target` | Target genome build: `hg19` or `hg38` | required |\
| `--liftover` | Lift data from hg18 to the target build (`yes` only for now) | optional, default: `yes` |\

simplified_run_stratification.sh
| Option | Description | Default |\
|---|---|---|\
| `--eur_marginX` | European (EUR) samples are filtered by setting a delimitation at the x-axis around the reference sample 1000G.
| optional, default : 0.3, a box is extended 30% wider around EUR samples from 1000G|\
| `--eur_marginY` | European (EUR) samples are filtered by setting a delimitation at the y-axis around the reference sample 1000G.
| optional, default : 0.3, a box is extended 30% wider around EUR samples from 1000G|\

Example : 

Merge the control and case datasets.
```bash
run_merge.sh --input_control WTCCC_control --input_case WTCCC_PS_case --ref_hg_target hg38
```

Perform stratification.
```bash
simplified_run_stratification.sh --ref_hg_target hg38 <--eur_marginX 0.3 (optional)>  <--eur_marginY 0.3 (optional)>
```

Run imputation. 
Perform stratification.
```bash
run_imputation.sh --ref_hg_target hg38 <--cpu_threads_phasing 40 (optional)> <--cpu_threads_imputation 14 (optional)>
```

Inputs are PLINK basenames (without `.bed/.bim/.fam`).

## Data availability

The WTCCC genotype data are under controlled access and are **not** distributed with this repository.
They can be requested from the [European Genome-phenome Archive (EGA)](https://ega-archive.org/).
