rule all:
	input:"filtered_bam/{s}.filtered.bam".format(s = config['sample'])

rule bwa_index:
	input: config['reference']
	output:
		"{ref}.amb".format(ref=config['reference']),
		"{ref}.ann".format(ref=config['reference']),
		"{ref}.bwt".format(ref=config['reference']),
		"{ref}.pac".format(ref=config['reference']),
		"{ref}.sa".format(ref=config['reference'])
	resources:
		mem=2,
		time=1
	shell:
		"bwa index {input}"

# align reads to reference genome with bwa
rule bwa_align:
	input:
		ref_index = rules.bwa_index.output,
		r = config["reads"].split(",") #splits sample file names into list by comma
	output:
		"filtered_bam/{sample}.filtered.bam"
	resources:
		mem=32,
		time=6
	threads: 8
	params:
		ref = config['reference'],
		qual=config['mapq'],
		# nm=config['n_mismatches']
	shell:
		"bwa mem -t {threads} {params.ref} {input.r} | "\
		"samtools view -b -q {params.qual} | "\
		# "bamtools filter -tag 'NM:<={params.nm}' | "\
		"samtools sort --threads {threads} -o {output}"
