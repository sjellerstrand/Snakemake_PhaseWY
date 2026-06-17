########################## Step 4 ##########################

# Check if whatshap is not disabled
if config.get("whatshap") != "OFF":

	# 4.1 Split vcf into one vcf per sample
	rule split_vcf:
		input:
			biallelic_vcf=intermediate + "biallelic.vcf.gz",
			bedfile=intermediate + "subsetted_genome.bed"
		output:
			fmt_vcf=intermediate + "whatshap/{sample}_formatted.vcf.gz"
		conda:
			"../envs/bcftools.yml"
		priority: 100
		shell:
			"""
			bcftools view {input.biallelic_vcf} -R {input.bedfile} -s {wildcards.sample} | bgzip -c > {output.fmt_vcf}
			tabix {output.fmt_vcf}
			"""

	# 4.2 Run whatshap on each bam using a sample specific vcf. 
	# Whatshap expects an index file in the same directory for each bam.
	# Note, it can run better if two cores are reserved (it uses
	# ~5-6GB of RAM) even though it runs only one core
	rule whatshap:
		input:
			fasta=intermediate + "subsetted_genome.fasta",
			fmt_vcf=intermediate + "whatshap/{sample}_formatted.vcf.gz",
			bam=sample_df.bamfile,
			idx=intermediate + "subsetted_genome.fasta.fai"
		output:
			whatshap_out=intermediate + "whatshap/{sample}_whatshap_out.vcf.gz"
		log:
			logs_dir + "whatshap/{sample}.log"
		conda:
			"../envs/whatshap.yml"
		priority: 100
		shell:
			"""
			whatshap phase -o {output.whatshap_out} --indels -r {input.fasta} {input.fmt_vcf} {input.bam} --sample {wildcards.sample} 2> {log}
			tabix {output.whatshap_out}
			"""

	# 4.3 Merge whatshap results then index the results.
	# Here, the list of files to merge is passed
	# to the shell command using params.files. The list of file names is
	# space separated.
	rule merge_whatshap:
		input:
			expand(intermediate + "whatshap/{sample}_whatshap_out.vcf.gz", sample=sample_df.index)
		output:
			intermediate + "whatshap/whatshap_out_merged.vcf.gz"
		params:
			files=lambda wildcards, input: ' '.join(input)
		conda:
			"../envs/bcftools.yml"
		priority: 100
		shell:
			"""
			bcftools merge {params.files} | bgzip -c > {output}
			tabix {output}
			"""

	# 4.4 Find heterogametic singletons
	# Find all singletons in heterogametic individuals
	rule find_singletons:
		input:
			sample_table=config["sample_table"],
			whatshap_vcf=intermediate + "whatshap/whatshap_out_merged.vcf.gz",
			bedfile=intermediate_clust_settings + "beds/unfiltered_sex_linked.bed",
		output:
			htgm_singletons_vcf=intermediate_clust_settings + "nonphased_singletons/whatshap_htgm_singletons.vcf.gz"
		conda:
			"../envs/bcfvcftools.yml"
		priority: 100
		shell:
			"""
			if [[ $(grep -vc '^#' {input.bedfile}) -gt 0 ]]
			then
				HETGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
				bcftools view {input.whatshap_vcf} -R {input.bedfile} | \
				vcffilter -f "AC = 1" | \
				bcftools view -s "$HETGAM" | \
				vcffilter -f "AC = 1" | \
				bgzip -c > {output.htgm_singletons_vcf}
			else
				bcftools view {input.whatshap_vcf} -h | bgzip -c > {output.htgm_singletons_vcf}
			fi
			tabix {output.htgm_singletons_vcf}
			"""

	# 4.5 List heterogametic singletons
	# List singletons in heterogametic individuals
	rule list_singletons:
		input:
			sample_table=config["sample_table"],
			htgm_singletons_vcf=intermediate_clust_settings + "nonphased_singletons/whatshap_htgm_singletons.vcf.gz"
		output:
			htgm_singletons_tsv=intermediate_clust_settings + "nonphased_singletons/whatshap_htgm_singletons.tsv",
			done=intermediate_clust_settings + "nonphased_singletons/whatshap_htgm_singletons.done"
		conda:
			"../envs/bcfvcftools.yml"
		priority: 100
		shell:
			"""
			HETGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
			echo -e "#CHROM\tPOS\tGT\tPS\tSample" > {output.htgm_singletons_tsv}

			echo $HETGAM | tr ',' '\n' | \
			while read IND
			do
				bcftools view {input.htgm_singletons_vcf} -s "$IND" | \
				vcffilter -f "AC = 1" | \
				bcftools query -f '%CHROM\t%POS[\t%GT\t%PS]\n' | \
				awk -F'\\t' -v OFS='\\t' -v ind=$IND '{{if($4 == ".") {{$4="None"}} print $1,$2,$3,$4,ind}}' \
				>> {output.htgm_singletons_tsv}
			done
			>> {output.done}
			"""

	# 4.6 Find heterogametic singletons
	# which are not part of a phase-set
	# and thus, not phased
	rule find_nonphased_singletons:
		input:
			sample_table=config["sample_table"],
			htgm_singletons_tsv=intermediate_clust_settings + "nonphased_singletons/whatshap_htgm_singletons.tsv",
			whatshap_vcf=intermediate + "whatshap/whatshap_out_merged.vcf.gz",
			done=intermediate_clust_settings + "nonphased_singletons/whatshap_htgm_singletons.done",
			bedfile=intermediate_clust_settings + "beds/unfiltered_sex_linked.bed"
		output:
			nonphased_htgm_singletons_bed=intermediate_clust_settings + "beds/nonphased_singletons_{sample}_singletons.bed",
			done=intermediate_clust_settings + "nonphased_singletons/whatshap_nonphased_{sample}_singleton_done.txt"
		params:
			outdir = intermediate_clust_settings + "nonphased_singletons/temp_singletons"
		conda:
			"../envs/bcfvcftools.yml"
		priority: 100
		shell:
			"""
			IND={wildcards.sample}
			mkdir -p {params.outdir}
			echo -e "#CHROM\tSTOP\tGT\tPS\tSample" > {output.nonphased_htgm_singletons_bed}
			if [[ $(grep -vc '^#' {input.bedfile}) -gt 0 && $(tail -n+2 {input.htgm_singletons_tsv} | awk -F'\\t' -v IND=$IND '{{if($5==IND) print}}' | wc -l) -gt 0 ]]
			then
				tail -n+2 {input.htgm_singletons_tsv} | awk -F'\\t' -v IND=$IND '{{if($5==IND) print}}' | \
				while IFS=$'\\t' read -r CHROM POS GT PS SAMPLE
				do
					if [ "$PS" != "None" ]
					then
						RANGE=$(echo $POS | awk '{{print $1-2500"-"$1+2500}}')
						bcftools view {input.whatshap_vcf} -s "$SAMPLE" -r "$CHROM":"$RANGE" -H | \
						awk -F'\\t' -v OFS='\\t' -v PS="$PS" '{{split($10, FORMAT, ":"); if(FORMAT[2] == PS) {{print $1,$2}}}}' \
						> {params.outdir}/"$SAMPLE"_PS_"$PS"_singleton_"$POS".tsv

						N_INF_VAR=$(bcftools view {input.whatshap_vcf} -R {params.outdir}/"$SAMPLE"_PS_"$PS"_singleton_"$POS".tsv | \
						vcffilter -f "AC > 1" | \
						bcftools view -H | wc -l)

						rm {params.outdir}/"$SAMPLE"_PS_"$PS"_singleton_"$POS".tsv

						if [ "$N_INF_VAR" -eq 0 ]
						then
							printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "$CHROM" "$POS" "$GT" "$PS" "$SAMPLE"
						fi

					else
						printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "$CHROM" "$POS" "$GT" "$PS" "$SAMPLE"
					fi

				done | awk -F'\\t' -v OFS="\\t" '{{if($1 ~ /^#CHROM/) {{col="START"}} else {{col=$2-1}}; print $1, col, $2, $3, $4, $5}}' >> {output.nonphased_htgm_singletons_bed}
			fi
			>> {output.done}
			"""

# If whatshap is disabled
else:

	# 4.7 Find all heterogametic singletons
	# since they cannot be accurately phased
	# to the correct sex chromosome without
	# Phase-set information
	rule find_hetgam_all_singletons:
		input:
			biallelic_vcf=intermediate + "biallelic.vcf.gz",
			bedfile=intermediate_clust_settings + "beds/unfiltered_sex_linked.bed",
			sample_table=config["sample_table"]
		output:
			nonphased_htgm_singletons_bed=results_clust_settings + "beds/nonphased_singletons_htgm_singletons.bed",
			done=intermediate_clust_settings + "beds/save_singletons.done"
		conda:
			"../envs/bcfvcftools.yml"
		shell:
			"""
			HETGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
			echo -e #CHROM\tSTART\tSTOP" > {output.nonphased_htgm_singletons_bed}
			bcftools view {input.biallelic_vcf} -R {input.bedfile} | \
			vcffilter -f "AC = 1" | \
			bcftools view -s "$HETGAM" | \
			vcffilter -f "AC = 1" | \
			bcftools view -H | \
			awk -F'\\t' -v OFS='\\t' '{{print $1,$2-1,$2}}' \
			>> {output.nonphased_htgm_singletons_bed}
			>> {output.done}
			"""
