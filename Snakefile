# -*- mode: Snakemake -*-

# Workflow for viral genomes from short read sequencing.

from pathlib import Path

# initialize

# sample sheet format:
#sample_name     paired  method  r1      r2

def read_samples(sample_list):
  sample_dict = {}
  added_samples = []
  for l in open(sample_list).readlines()[1:]: # use pandas
    fields = l.split()
    if fields[0] in added_samples:
      raise RuntimeError("Not all samples have unique names:" + fields[0] + "appears twice.")
    else:
      added_samples.append(fields[0])
    sample_dict[fields[0]] = {'paired':fields[1] == "True", 'method':fields[2], 'r1':Path(fields[3]), 'r2':Path(fields[4])}
  return sample_dict

sample_dict = read_samples(Path(config["all"]["sample_list"]))
output_dir = Path(config["all"]["root_dir"])/Path(config["all"]["output_dir"])

# align
rule build_bt2_index:
  input: config["align"]["target_fasta"]
  output: output_dir/'db'/'bt2'/'target.1.bt2'
  params: target = str(output_dir/'db'/'bt2'/'target')
  shell:
    """
    bowtie2-build {input} {params.target}
    """

rule align_bt2:
  input:
    index_generated = rules.build_bt2_index.output,
    r1 = lambda wildcards: sample_dict[wildcards.sample]['r1'],
    r2 = lambda wildcards: sample_dict[wildcards.sample]['r2']
  output: output_dir/'align'/'bt2'/'{sample}.bam'
  params:
    bt2_index = str(output_dir/'db'/'bt2'/'target'),
    sam = str(output_dir/'align'/'bt2'/'{sample}.sam')
  threads: 10
  shell:
    """
    bowtie2 -x {params.bt2_index} -1 {input.r1} -2 {input.r2} -S {params.sam} --threads {threads}
    samtools view -u -@ {threads} -o {output} {params.sam}
    rm {params.sam}
    """

rule unify_alignments_bt2:
  input: rules.align_bt2.output
  output: str(output_dir/'align'/'{sample}.bam')
  shell: 
    """
    ln -sr {input} {output}
    """

rule all_align:
  input: expand(output_dir/'align'/'{sample}.bam', sample = sample_dict.keys())

# summarize
rule sort_index_bam:
  input: str(output_dir/'align'/'{sample}.bam')
  output:
    sorted_bam = str(output_dir/'align'/'{sample}.bam.sorted'),
    index = str(output_dir/'align'/'{sample}.bam.sorted.bai')
  shell:
    """
    samtools sort -o {output.sorted_bam} {input}
    samtools index {output.sorted_bam} {output.index}
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
  shell:
    """
    bcftools mpileup -Ou --max-depth 10000 --max-idepth 10000 --annotate FORMAT/AD -f {input.ref} {input.sorted_bam} | bcftools call -mv -Oz -o {params.prefix}_unfilled_calls.vcf.gz
    bcftools index {params.prefix}_unfilled_calls.vcf.gz
    bcftools norm -f {input.ref} {params.prefix}_unfilled_calls.vcf.gz -Ob -o {params.prefix}_calls.norm.bcf
    bcftools filter --IndelGap 5 {params.prefix}_calls.norm.bcf -Ob -o {params.prefix}_calls.norm.flt-indels.bcf
    bcftools +fill-tags {params.prefix}_unfilled_calls.vcf.gz -Oz -o {output.compressed_vcf} -- -t VAF
    bcftools index {output.compressed_vcf}
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

