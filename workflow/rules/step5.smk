########################## Step 5 ##########################

# 5.1 Run shapeit4 on the merged whatshap output. 
# This command takes a list of contigs and runs shapeit4 on 
# each contig at a time.
rule shapeit4:
	input:
		vcf = lambda wildcards: (
			intermediate + "biallelic.vcf.gz"
			if config.get("whatshap") == "OFF"
			else intermediate + "whatshap/whatshap_out_merged.vcf.gz"),
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		run_complete=intermediate + "shapeit4/shapeit4_complete.txt"
	params:
		log=logs_dir + "shapeit4/",
		outdir=intermediate + "shapeit4",
		use_ps = "" if config.get("disable_whatshap") == "OFF" else "--use-PS 0.0001"
	threads: 8
	conda:
		"../envs/shapeit4.yml"
	shell:
		"""
		mkdir -p {params.log}
		mkdir -p {params.outdir}
		while read line
		do
			contig=$(echo $line | cut -f1)
			shapeit4 --input {input.vcf} \
				--thread 1 \
				--seed 123456 \
				--region "$contig" \
				{params.use_ps} \
				--sequencing \
				--output {params.outdir}/"$contig"_shapeit4_out.vcf.gz \
				--log {params.log}/"$contig"_phased.log
			done <{input.filt_list}
		>{output.run_complete}
		"""

# 5.2 Modify vcf after shapeit4 run. There is a dependancy to ensure
# that shapeit4 completes first. 
rule modify_vcf:
	input:
		shapeit4_done=intermediate + "shapeit4/shapeit4_complete.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt"
	params:
		indir=intermediate + "shapeit4",
		outdir=intermediate + "shapeit4/modified"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		mkdir -p {params.outdir}
		while read line
		do
			contig=$(echo $line | cut -f1)
			bcftools view {params.indir}/"$contig"_shapeit4_out.vcf.gz -Ou | bcftools norm -m +any | vcffixup - | \
			awk -F'\\t' 'BEGIN {{OFS=FS}} /^[#]/ {{print; next}} {{for (i=10; i<=NF; i++) {{ gsub("/","|",$i)}} print}}' | \
			bgzip -c > {params.outdir}/"$contig"_phased_all_variants.vcf.gz
			tabix {params.outdir}/"$contig"_phased_all_variants.vcf.gz
			done <{input.filt_list}
		cp {input.shapeit4_done} {output.mod_done}
		"""

