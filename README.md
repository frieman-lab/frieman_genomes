# frieman_genomes

Snakemake pipeline for reference-guided genome assembly. Requires `conda` and `snakemake` installed. The pipeline handles the rest of the dependencies at runtime.

Each pipeline run needs an updated `config.yml` and `samples.tsv` file. Templates are provided, as well as an [example with more detailed instructions](example/)

Example run on a small dataset:

```
git clone https://github.com/louiejtaylor/frieman_genomes
cd example/
snakemake all_summarize --use-conda -p --snakefile ../Snakefile --configfile example_config.yml --cores 1
```

Alternately, on a cluster (Slurm in the below example) the command could look like this:

```
snakemake all_summarize --use-conda -p --snakefile ../Snakefile --configfile example_config.yml --jobs 1 --cores 1 --latency-wait 30 --cluster "sbatch --mem 10G -c {threads} "
```

Used by [Frieman lab](https://www.medschool.umaryland.edu/profiles/frieman-matthew/) members for genome assembly of sequenced laboratory stocks.
