###############################################
###### Snakefile for ICRA and SGV Finder ######
###############################################

from os.path import join, isfile
from os import listdir
import sys
import glob

# Specify project directories

DATA_DIR = config("preprocessed_reads_directory")
PROJECT_DIR = config("output_directory")
SGV_DIR = config("SGV_code_directory")

# Get file names

FILES = [f for f in listdir(DATA_DIR) if (isfile(join(DATA_DIR,f)) and f.endswith('.fq.gz') and not f.startswith('._'))]

SAMPLE_PREFIX = [f.replace('_1.fq.gz','') for f in FILES]
SAMPLE_PREFIX = [f.replace('_2.fq.gz','') for f in SAMPLE_PREFIX]

#################################################

rule all:
	input:
		expand(join(PROJECT_DIR, "00_unzip/{sample}_{read}.fastq",sample=SAMPLE_PREFIX, read=['1','2'])),
		expand(join(PROJECT_DIR, "01_icra/{sample}.jsdel",sample=SAMPLE_PREFIX),
		expand(join(PROJECT_DIR, "01_icra/{sample}.pmp",sample=SAMPLE_PREFIX),
		expand(join(PROJECT_DIR, "02_sgv_perfile/{sample}.txt", sample=SAMPLE_PREFIX),
		join(PROJECT_DIR, "03_sgv/deletion-sgv.csv"),
		join(PROJECT_DIR, "03_sgv/variable-sgv.csv")

#################################################

rule unzip:
	input: expand(join(DATA_DIR,"{sample}_{read}.fq.gz",sample=SAMPLE_PREFIX,read=['1','2']))
	output: join(PROJECT_DIR,"00_unzip/{sample}_{read}.fastq")
	params:
		unzipped = join(DATA_DIR,"{sample}_{read}.fq",sample=SAMPLE_PREFIX,read=['1','2'])),
		outdir = join(PROJECT_DIR,"00_unzip/")
	shell: """
		mkdir -p {params.outdir}
		gunzip {input}
		cp {params.unzipped} {params.outdir}
		for f in {params.outdir};do mv "$f" "${f/fq/fastq}";done
	"""

#################################################

rule icra:
	input: expand(join(PROJECT_DIR, "00_unzip/{sample}_{read}.fastq"))
	output: expand(join(PROJECT_DIR, "01_icra/{sample}.jsdel"), expand(join(PROJECT_DIR, "00_icra/{sample}.pmp"))
	params:
		outfol = join(PROJECT_DIR,"01_icra/"),
		fq = expand(join(PROJECT_DIR,"00_unzip/{sample}")),
		max_mismatch = config['icra']['max_mismatch'],
		ignore_lengths = config['icra']['ignore_lengths'],
		epsilon = config['icra']['epsilon'],
		max_iterations = config['icra']['max_iterations'],
		min_bins = config['icra']['min_bins'],
		max_bins = config['icra']['max_bins'],
		min_reads = config['icra']['min_reads'],
		dense_region_coverage = config['icra']['dense_region_coverage'],
		length_minimum = config['icra']['length_minimum'],
		length_maximum = config['icra']['length_maximum'],
		usage = config['icra']['usage'],
	threads: 4
	resources:
		time = 10,
		mem = 32
	shell: """
		mkdir {PROJECT_DIR}/01_icra
		python {SGV_DIR}/ICRA_cmd.py {params.outfol} {input} --pe --max_mismatch {params.max_mismatch} --ignore_lengths {params.ignore_lengths} --epsilon {params.epsilon} --max_iterations {params.max_iterations} --min_bins {params.min_bins} --max_bins {params.max_bins} --min_reads {params.min_reads} --dense_region_coverage {params.dense_region_coverage} --length_minimum {params.length_minimum} --length_maximum {params.length_maximum} --usage {params.usage} --use_theta 
	"""

#################################################

rule sgv_perfile:
	input: expand(join(PROJECT_DIR, "01_icra/{sample}.jsdel"))
	output: expand(join(PROJECT_DIR, "02_sgv_perfile/{sample}.txt"))
	params:
		x_coverage = config['sgv_perfile']['x_coverage'],
		rate_param = config['sgv_perfile']['rate_param']
	threads: 4
	resources:
		time = 10,
		mem = 32
	shell: """
		mkdir {PROJECT_DIR}/02_sgv_perfile
		AVG=$(awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("%f",m);}'  {sample}_{read}.fastq)
		python {SGV_DIR}/SGVF_PerFile_cmd.py {input} {output} $AVG --x_coverage {params.x_coverage} --rate_param {params.rate_param}
	"""

#################################################

rule sgv_finder
	input: join(PROJECT_DIR,"02_sgv_perfile/")
	output:
		de = join(PROJECT_DIR, "03_sgv/deletion-sgv.csv"),
		var = join(PROJECT_DIR, "03_sgv/deletion-sgv.csv")
	params:
		x_coverage = config['sgv_finder']['x_coverage'],
		rate_param = config['sgv_finder']['rate_param'],
		min_samp_cutoff = config['sgv_finder']['min_samp_cutoff'],
		dels_detect_thresh = config['sgv_finder']['dels_detect_thresh'],
		real_del_thresh = config['sgv_finder']['real_del_thresh'],
		vsgv_dissim_thresh = config['sgv_finder']['vsgv_dissim_thresh'],
		dels_cooc_thresh = config['sgv_finder']['dels_cooc_thresh'],
		vsgv_clip_quantile = config['sgv_finder']['vsgv_clip_quantile'],
		vsgv_fit_interval = config['sgv_finder']['vsgv_fit_interval'],
		vsgv_fit_method = config['sgv_finder']['vsgv_fit_method'],
		vsgv_dense_perc = config['sgv_finder']['vsgv_dense_perc']
	threads: 4 
	resources:
		time = 10,
		mem = 32
	shell: """
		mkdir {PROJECT_DIR}/03_sgv/
		python {SGV_DIR}/SGV_cmd.py {input}/*.txt {output.de} {output.var} --x_coverage {params.x_coverage} --rate_param {params.rate_param} --byorig --min_samp_cutoff {params.min_samp_cutoff} --dels_detect_thresh {params.dels_detect_thresh} --real_del_thresh {params.real_del_thresh} --vsgv_dissim_thresh {params.vsgv_dissim_thresh} --dels_cooc_thresh {params.dels_cooc_thresh} --vsgv_clip_quantile {params.vsgv_clip_quantile} --vsgv_fit_interval {params.vsgv_fit_interval} --vsgv_fit_method {params.vsgv_fit_method} --vsgv_dense_perc {params.vsgv_dense_perc} --csv_output
	"""














