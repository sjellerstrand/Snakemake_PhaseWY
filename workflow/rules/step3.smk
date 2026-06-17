########################## Step 3 ##########################

# 3.1. Additional filtering of vcf. This is making a table 
# that will be used to identify contigs that have 2 or more 
# heterozygous sites. The table contains 
# the observed numbers of homozygotes and heterozygotes.
# Here, the original vcf is used.
rule summarize_vcf:
	input:
		input_vcf=config["input_vcf"],
		index=intermediate + "subsetted_genome.fasta.fai"
	output:
		bedfile=intermediate + "subsetted_genome.bed",
		allele_freqs=intermediate + "allele_freqs.tsv"
	params:
		outdir=intermediate,
		min_dp=config["min_dp"],
		missing=config["missing"],
		min_mean=config["min_mean"],
		max_mean=config["max_mean"],
		mask=config["mask"]
	conda:
		"../envs/bcfvcftools.yml"
	priority: 100
	shell:
		"""
		awk -F'\\t' -v OFS='\\t' '{{print $1, 0, $2}}' {input.index} > {output.bedfile}
		if [[ -n "{params.mask}" && -s "{params.mask}" && $(grep -vc '^#' {params.mask}) -gt 0 ]]
		then
			bcftools view {input.input_vcf} -R {output.bedfile} | \
			vcftools --vcf - --exclude-bed {params.mask} --minDP {params.min_dp} --max-missing {params.missing} --min-meanDP {params.min_mean} --max-meanDP {params.max_mean} \
			--hardy --stdout >{output.allele_freqs}
		else
			bcftools view {input.input_vcf} -R {output.bedfile} | \
			vcftools --vcf - --minDP {params.min_dp} --max-missing {params.missing} --min-meanDP {params.min_mean} --max-meanDP {params.max_mean} \
			--hardy --stdout >{output.allele_freqs}
		fi
		"""

# 3.2 Make a bedfile from the genome assembly and make
# a list of filtered contigs to be used for shapeit4.
rule make_bed_filt_contigs:
	input:
		allele_freqs=intermediate + "allele_freqs.tsv",
		bedfile=intermediate + "subsetted_genome.bed"
	output:
		filt_list=intermediate + "filtered_contigs_list.txt",
		nonphased=intermediate + "nonphased_contigs.bed"
	priority: 100
	shell:
		"""
		awk -F'\\t' 'NR>1 {{split($3,a,"/"); if(a[2]>=2) count[$1]++}} END {{for (c in count) if(count[c]>=2) print c}}' {input.allele_freqs} > {output.filt_list}
		if [ -s {output.filt_list} ]
		then
			grep -vFf {output.filt_list} {input.bedfile} > {output.nonphased} || touch {output.nonphased}
		else
			touch {output.nonphased}
		fi
		"""

# 3.3 Filter vcf 
rule format_vcf:
	input:
		input_vcf=config["input_vcf"],
		bedfile=intermediate + "subsetted_genome.bed"
	output:
		biallelic_vcf=intermediate + "biallelic.vcf.gz",
		done=intermediate + "biallelic_done.txt"
	params:
		min_dp=config["min_dp"],
		missing=config["missing"],
		min_mean=config["min_mean"],
		max_mean=config["max_mean"],
		mask=config["mask"]
	conda:
		"../envs/bcfvcftools.yml"
	priority: 100
	shell:
		"""
		if [[ -n "{params.mask}" && -s "{params.mask}" && $(grep -vc '^#' {params.mask}) -gt 0 ]]
		then
			bcftools view {input.input_vcf} -R {input.bedfile} | \
			vcftools --vcf - --exclude-bed {params.mask} --minDP {params.min_dp} --max-missing {params.missing} --min-meanDP {params.min_mean} --max-meanDP {params.max_mean} \
			--recode --recode-INFO-all --stdout | \
			bcftools annotate -x FORMAT,INFO -Ou | bcftools norm -m - -Ov | bgzip -c > {output.biallelic_vcf}
		else
			bcftools view {input.input_vcf} -R {input.bedfile} | \
			vcftools --vcf - --minDP {params.min_dp} --max-missing {params.missing} --min-meanDP {params.min_mean} --max-meanDP {params.max_mean} \
			--recode --recode-INFO-all --stdout | \
			bcftools annotate -x FORMAT,INFO | vcffixup - | bcftools norm -m - | bgzip -c > {output.biallelic_vcf}
		fi
		tabix {output.biallelic_vcf}
		>> {output.done}
		"""