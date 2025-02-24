# -*- mode: Snakemake -*-

# Workflow for viral genomes from short or long read sequencing.

from pathlib import Path

# initialize

# sample sheet format:
#sample_name     paired  method  r1      r2
#              True/False short/long

def read_samples(sample_list):
  sample_dict = {}
  added_samples = []
  for l in open(sample_list).readlines()[1:]: # use pandas
    fields = l.split()
    if fields[0] in added_samples:
      raise RuntimeError("Not all samples have unique names:" + fields[0] + "appears twice.")
    else:
      added_samples.append(fields[0])
    # don't need fields[1]--unpaired if no r2 given
    sample_dict[fields[0]] = {'paired':fields[1] == "True", 'method':fields[2], 'r1':Path(fields[3]), 'r2':Path(fields[4])}
  return sample_dict

sample_dict = read_samples(Path(config["all"]["sample_list"]))
output_dir = Path(config["all"]["output_dir"])

# build indices
rule build_bt2_index:
  input: config["align"]["target_fasta"]
  output: output_dir/'db'/'bt2'/'target.1.bt2'
  params: target = str(output_dir/'db'/'bt2'/'target')
  threads: 10
  conda: "envs/bowtie2.yml"
  shell:
    """
    bowtie2-build {input} {params.target} --threads {threads}
    """

rule build_mm2_index:
  input: config["align"]["target_fasta"]
  output: output_dir/'db'/'mm2'/'target.mmi'
  params: opt = "map-ont" # option to make user-configurable
  threads: 10
  conda: "envs/minimap2.yml"
  shell:
    """
    minimap2 -x {params.opt} -t {threads} -d {output} {input}
    """

# align
rule align_bt2:
  input:
    index_generated = rules.build_bt2_index.output,
    r1 = lambda wildcards: sample_dict[wildcards.sample]['r1'],
    r2 = lambda wildcards: sample_dict[wildcards.sample]['r2']
  output: temp(output_dir/'align'/'bt2'/'{sample}.sam')
  params:
    bt2_index = str(output_dir/'db'/'bt2'/'target')
  threads: 10
  conda: "envs/bowtie2.yml"
  shell:
    """
    bowtie2 -x {params.bt2_index} -1 {input.r1} -2 {input.r2} -S {output} --threads {threads}
    """

rule align_mm2:
  input:
    index = rules.build_mm2_index.output,
    r1 = lambda wildcards: sample_dict[wildcards.sample]['r1'],
    r2 = lambda wildcards: sample_dict[wildcards.sample]['r2']
  output: temp(output_dir/'align'/'mm2'/'{sample}.sam')
  params:
    opt = "map-ont"
  threads: 10
  conda: "envs/minimap2.yml"
  shell:
    """
    minimap2 -a {input.index} -t {threads} {input.r1} {input.r2} > {output}
    """

methods_map = {'short': 'bt2', 'long': 'mm2', 'ont': 'mm2', 'illumina': 'bt2'}

# fork between illumina/ont here
rule sam_to_bam:
  input: bam = lambda wildcards: str(output_dir/'align'/methods_map[sample_dict[wildcards.sample]['method']]/wildcards.sample)+'.sam'
  output: temp(str(output_dir/'align'/'{sample}.bam'))
  threads: 10
  conda: "envs/process_alignments.yml"
  shell:
    """
    samtools view -u -@ {threads} -o {output} {input}
    """

rule all_align:
  input: expand(output_dir/'align'/'{sample}.bam', sample = sample_dict.keys())

# summarize
rule sort_index_bam:
  input:
    bam = str(output_dir/'align'/'{sample}.bam')
  output:
    sorted_bam = str(output_dir/'align'/'{sample}.bam.sorted'),
    index = str(output_dir/'align'/'{sample}.bam.sorted.bai')
  threads: 10
  conda: "envs/process_alignments.yml"
  shell:
    """
    samtools sort -@ {threads} -o {output.sorted_bam} {input.bam}
    samtools index -@ {threads} {output.sorted_bam} {output.index}
    """

rule generate_vcf:
  input:
    ref = config["align"]["target_fasta"],
    sorted_bam = str(output_dir/'align'/'{sample}.bam.sorted')
  output:
    vcf = str(output_dir/'consensus'/'{sample}'/'{sample}_calls.vcf'),
    compressed_vcf = str(output_dir/'consensus'/'{sample}'/'intermediates'/'{sample}_calls.vcf.gz')
  params:
    prefix = str(output_dir/'consensus'/'{sample}'/'intermediates'/'{sample}')
  threads: 10
  conda: "envs/process_alignments.yml"
  shell:
    """
    bcftools mpileup -Ou --max-depth 10000 --max-idepth 10000 --annotate FORMAT/AD -f {input.ref} {input.sorted_bam} | bcftools call -mv -Oz -o {params.prefix}_unfilled_calls.vcf.gz
    bcftools index {params.prefix}_unfilled_calls.vcf.gz --threads {threads}
    bcftools norm -f {input.ref} {params.prefix}_unfilled_calls.vcf.gz -Ob -o {params.prefix}_calls.norm.bcf --threads {threads}
    bcftools filter --IndelGap 5 {params.prefix}_calls.norm.bcf -Ob -o {params.prefix}_calls.norm.flt-indels.bcf --threads {threads}
    bcftools +fill-tags {params.prefix}_unfilled_calls.vcf.gz -Oz -o {output.compressed_vcf} -- -t VAF
    bcftools index {output.compressed_vcf} --threads {threads}
    gzip -cd {output.compressed_vcf} > {output.vcf}
    """

# credit to https://www.biostars.org/p/367626/#417214
rule generate_consensus:
  input: 
    calls = str(output_dir/'consensus'/'{sample}'/'intermediates'/'{sample}_calls.vcf.gz'),
    ref = config["align"]["target_fasta"]
  output: str(output_dir/'consensus'/'{sample}'/'{sample}_raw.fasta')
  params:
    sample = "{sample}"
  conda: "envs/process_alignments.yml"
  shell:
    """
    bcftools consensus {input.calls} -p {params.sample}_raw -f {input.ref} > {output}
    """

rule symlink_summaries:
  input:
    fasta = str(output_dir/'consensus'/'{sample}'/'{sample}_raw.fasta'),
    sorted_bam = str(output_dir/'align'/'{sample}.bam.sorted'),
    index = str(output_dir/'align'/'{sample}.bam.sorted.bai'),
    vcf = str(output_dir/'consensus'/'{sample}'/'{sample}_calls.vcf')
  output:
    fasta = str(output_dir/'summary'/'{sample}'/'{sample}_raw.fasta'),
    sorted_bam = str(output_dir/'summary'/'{sample}'/'{sample}_sorted.bam'),
    index = str(output_dir/'summary'/'{sample}'/'{sample}_sorted.bam.bai'),
    vcf = str(output_dir/'summary'/'{sample}'/'{sample}_calls.vcf')
  shell:
    """
    ln -sr {input.fasta} {output.fasta}
    ln -sr {input.sorted_bam} {output.sorted_bam}
    ln -sr {input.index} {output.index}
    ln -sr {input.vcf} {output.vcf}
    """

rule all_summarize:
  input:
    fasta = expand(output_dir/'summary'/'{sample}'/'{sample}_raw.fasta', sample = sample_dict.keys()),
    sorted_bam = expand(output_dir/'summary'/'{sample}'/'{sample}_sorted.bam', sample = sample_dict.keys()),
    index = expand(output_dir/'summary'/'{sample}'/'{sample}_sorted.bam.bai', sample = sample_dict.keys()),
    vcf = expand(output_dir/'summary'/'{sample}'/'{sample}_calls.vcf', sample = sample_dict.keys())

