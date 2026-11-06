# StructMuMMI 
StructMuMMI is a framework organized in two blocks: the preprocessing and the models training

Preprocessing pipeline for WTCCC genotype data (PLINK format):

1. Sex check and variant QC on case and control datasets
2. Allele/strand alignment and merging of case and control
3. Liftover from hg18 (NCBI36) to hg19 (GRCh37) or hg38 (GRCh38), with reference-allele fixing (`bcftools +fixref`)
4. (Optional, in progress) Allele harmonization to the 1000 Genomes reference panel

## Repository layout

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

### 2. Create the conda environments

```bash
conda env create -f gcta_env.yml       # creates "gcta_env" (plink, gcta)
conda env create -f liftover_env.yml   # creates "liftover_env" (liftOver)
```

### 3. Install the external tools

```bash
# plink2
mkdir -p external_tools/plink
# download the Linux build from https://www.cog-genomics.org/plink/2.0/
# unzip it and place the binary at external_tools/plink/plink2
chmod +x external_tools/plink/plink2

# bcftools 1.23.1
cd external_tools
wget https://github.com/samtools/bcftools/releases/download/1.23.1/bcftools-1.23.1.tar.bz2
tar -xjf bcftools-1.23.1.tar.bz2
cd bcftools-1.23.1 && ./configure && make && cd ../..
```

### 4. Download the reference data
The processed 1000 Genomes reference panels (hg19 and hg38) are available on Zenodo: `<DOI>`.

```bash
bash download_references.sh
```

They are placed in `references/`.

### 5.  Run the scripts from anywhere

Run this from the folder that contains `structmummi/`:

```bash
chmod +x structmummi/preprocessing/*.sh
echo "export PATH=\"$(pwd)/structmummi/preprocessing:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

## Usage

Run the pipeline from a **separate working folder**: all outputs are written to the current directory.

```bash
mkdir -p ~/my_analysis && cd ~/my_analysis

run_merge.sh \
    --input_control /path/to/control \
    --input_case    /path/to/case \
    --ref_hg_target hg19
```

(Without step 5, use `bash /path/to/structmummi/preprocessing/run_merge.sh ...`.)

Inputs are PLINK basenames (without `.bed/.bim/.fam`).

### Options

| Option | Description | Default |
|---|---|---|
| `--input_control` | PLINK basename of the control dataset | required |
| `--input_case` | PLINK basename of the case dataset | required |
| `--ref_hg_target` | Target genome build: `hg19` or `hg38` | required |
| `--liftover` | Lift data from hg18 to the target build (`yes` only for now) | `yes` |

### Environment variables

| Variable | Description | Default |
|---|---|---|
| `PLINK2` | Path to the plink2 binary | `external_tools/plink/plink2`, otherwise `plink2` in PATH |
| `THREADS` | Number of threads | all CPUs |

Example: `THREADS=8 run_merge.sh ...`

## Output

`allvariants.{bed,bim,fam}`: merged case/control dataset, lifted to the target build, with variant IDs `chr:pos:ref:alt` and phenotypes (1 = control, 2 = case).

## Data availability

The WTCCC genotype data are under controlled access and are **not** distributed with this repository.
They can be requested from the [European Genome-phenome Archive (EGA)](https://ega-archive.org/).
