########################### Step 6 ##########################

# 6.1 Generate a file filtered by minor allele count. The number
# of genomes is currently set to the sum of the ploidy column
# in the samples.tsv file.
rule make_mac_vcf:
	input:
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		mac_done=intermediate_clust_settings + "sex_linkage/mac_vcfs/mac_done.txt"
	params:
		indir=intermediate + "shapeit4/modified",
		outdir=intermediate_clust_settings + "sex_linkage/mac_vcfs",
		mac=mac,
		num_genomes=sum(sample_df.ploidy)
	log:
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		num_genomes="$(({params.num_genomes}-{params.mac}))"
		mkdir -p {params.outdir}
		while read line
		do
			contig=$(echo $line | cut -f1)
			bcftools view {params.indir}/"$contig"_phased_all_variants.vcf.gz | \
			vcftools --vcf - \
			--max-alleles 2 \
			--recode --recode-INFO-all --stdout | \
			vcffilter -f "AC > {params.mac} & AC < $num_genomes" | \
			bgzip -c > {params.outdir}/"$contig"_phased_all_variants_mac.vcf.gz
		done <{input.filt_list}
		cp {input.mod_done} {output.mac_done}
		"""

# 6.2 Run R script sex_linkage.r to determine sex-linkage of 
# each haplotype. Takes the mac VCFs from make_mac_vcf as input
# along with various parameters.

# 6.2.1 Large contigs
rule sex_linkage_large:
	input:
		sex_depth_windows_done=intermediate_clust_settings + "sex_depth_windows/{large_contigs}/{large_contigs}_sex_depth_windows_clump.done",
		mac_done=intermediate_clust_settings + "sex_linkage/mac_vcfs/mac_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		index=intermediate + "subsetted_genome.fasta.fai",
		sample_table=config["sample_table"]
	output:
		done=intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt"
	params:
		sex_linkage="workflow/scripts/sex_linkage.r",
		indir1=intermediate_clust_settings + "sex_linkage/mac_vcfs",
		indir2=intermediate_clust_settings + "sex_depth_windows",
		outdir=intermediate_clust_settings + "sex_linkage",
		depth_diff=sex_depth_threshold,
		model=dist_model
	log:
	conda:
		"../envs/r-sexlinkage.yml"
	shell:
		"""
		contig="{wildcards.large_contigs}"
		if grep -qx "$contig" {input.filt_list}
		then
			contiglength=$(grep -w "$contig" {input.index} | cut -f2)
			Rscript {params.sex_linkage} --no-save --args INDS={input.sample_table} \
			CONTIG="$contig" CONTIG_LENGTH="$contiglength" INDIR1={params.indir1} INDIR2={params.indir2} \
			OUTDIR={params.outdir}/"$contig" SEX_DEPTH_THRESH={params.depth_diff} MODEL={params.model}
		fi
		>> {output.done}
		"""

# 6.2.2  Clumped contigs
rule sex_linkage_clump:
	input:
		sex_depth_windows_done=intermediate_clust_settings + "sex_depth_windows/sex_depth_windows_clump.done",
		mac_done=intermediate_clust_settings + "sex_linkage/mac_vcfs/mac_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		index=intermediate + "subsetted_genome.fasta.fai",
		sample_table=config["sample_table"],
		clump_list=intermediate + "clump_contigs.txt"
	output:
		done=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt"
	params:
		sex_linkage="workflow/scripts/sex_linkage.r",
		indir1=intermediate_clust_settings + "sex_linkage/mac_vcfs",
		indir2=intermediate_clust_settings + "sex_depth_windows",
		outdir=intermediate_clust_settings + "sex_linkage",
		depth_diff=sex_depth_threshold,
		model=dist_model
	log:
	conda:
		"../envs/r-sexlinkage.yml"
	shell:
		"""
		if [ -s {input.clump_list} ]
		then
			grep -Ff {input.filt_list} {input.clump_list} > {params.outdir}/filt_clump_contigs.txt || true
			if [[ -s {params.outdir}/filt_clump_contigs.txt && $(grep -vc '^#' {params.outdir}/filt_clump_contigs.txt) -gt 0 ]]
			then
				while read line
				do
					contig=$(echo $line | cut -f1)
					mkdir -p {params.outdir}/"$contig"
					contiglength=$(grep -w "$contig" {input.index} | cut -f2 )
					Rscript {params.sex_linkage} --no-save --args INDS={input.sample_table} \
					CONTIG="$contig" CONTIG_LENGTH="$contiglength" INDIR1={params.indir1} INDIR2={params.indir2} \
					OUTDIR={params.outdir}/"$contig" SEX_DEPTH_THRESH={params.depth_diff} MODEL={params.model}
				done < {params.outdir}/filt_clump_contigs.txt
			fi
		fi
		>> {output.done}
		"""

# 6.3 Summarise sex linked regions in bed file
rule sexlinked_bed:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sexlinked_bed_done=intermediate_clust_settings + "sex_linkage/sexlinked_bed_done.txt"
	params:
		indir=intermediate_clust_settings + "sex_linkage"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(cat {params.indir}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				cat {params.indir}/"$contig"/"$contig"_phase_windows.bed | grep -E "#|Sex-linked" | cut -f1,2,3 | bedtools merge \
				> {params.indir}/"$contig"/"$contig"_sex_linked.bed
			fi
		done < {input.filt_list}
		>> {output.sexlinked_bed_done}
		"""

# 6.4 Identify problematic sex-linked sites. If heterogametes
# show homozygotic genotypes for both alleles it is likely a sign 
# of heterogametic drop out, unreliable identification of sex-linkage,
#or incomplete lineage sorting.
rule identify_sexlinked_ILS1_homozygotes:
	input:
		sexlinked_bed_done=intermediate_clust_settings + "sex_linkage/sexlinked_bed_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		sexlinked_homozygotes_done=intermediate_clust_settings + "sex_linkage/sexlinked_homozygotes_done.txt"
	params:
		indir1=intermediate + "shapeit4/modified",
		indir2=intermediate_clust_settings + "sex_linkage"
	log:
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		HETGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [[ -s {params.indir2}/"$contig"/"$contig"_sex_linked.bed && $(grep -vc '^#' {params.indir2}/"$contig"/"$contig"_sex_linked.bed) -gt 0 ]]
			then
				bcftools view {params.indir1}/"$contig"_phased_all_variants.vcf.gz \
				-s $HETGAM -R {params.indir2}/"$contig"/"$contig"_sex_linked.bed | \
				vcffilter -f "AC > 0 & AF < 1" | \
				vcftools --vcf - --hardy --stdout | \
				tail -n+2 | cut -f1,2,3 | awk -F'\\t|/' -v OFS='\\t' '{{ if($3 > 0 && $5 > 0) print $1,$2-1,$2 }}' \
				> {params.indir2}/"$contig"/"$contig"_sex_linked_ILS1.bed
			fi
		done <{input.filt_list}
		>> {output.sexlinked_homozygotes_done}
		"""

