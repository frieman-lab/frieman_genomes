# Test dataset

This is a very small, test dataset that should run on most machines. 

Associated files:

## test_config.yml

This file contains the configuration required for the pipeline. The pipeline needs paths to a sample list (below), where to save the output, and a .fasta-formatted target to align to.

```
all:
  root_dir: "."
  output_dir: "output"
  sample_list: "test_samples.tsv"

align:
  target_fasta: "../reference/MN985325_1.fasta"
```

## test_samples.yml

This file contains sample information in a tab-separated format. Paths to paired read files as well as information on the data collection strategy (e.g. long or short read sequencing) are included below.

```
sample_name	paired	method	r1	r2
paired1	True	short	data/paired1_r1.fastq.gz	data/paired1_r2.fastq.gz
paired2	True	long	data/paired2_r1.fastq.gz	data/paired2_r2.fastq.gz
```

To run starting from zero on a standard Linux/Ubuntu distro:

```
# setup conda

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ${HOME}/miniconda3
export PATH=${PATH}:${HOME}/miniconda3/bin

# install dependencies through conda

conda install snakemake mamba -c conda-forge -c bioconda

# download and run the pipeline on the test dataset

git clone https://github.com/louiejtaylor/frieman_genomes
cd test/
snakemake all_summarize --use-conda -p --snakefile ../Snakefile --configfile test_config.yml
```

