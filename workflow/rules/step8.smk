########################## Step 8 ##########################

# 8.1 Combine and sort callable regions
rule cat_regions:
	input:
		done1=lambda wildcards: expand(intermediate + "samtools/{large_contigs}/{large_contigs}_callable_regions.bed", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate + "samtools/clump_contigs_callable_regions.done"
	output:
		results_clust_settings + "beds/callable_regions.bed"
	params:
		indir=intermediate + "samtools",
		index=intermediate + "subsetted_genome.fasta.fai"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cat {params.indir}/*/*_callable_regions.bed | bedtools sort -g {params.index} > {output}
		"""

# 8.2 Concatenate regions of the genome 
# with hetrogametic drop-out due to 
# sex depth differences
rule cat_htgm_dropout:
	input:
		done1=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate + "sex_depth_difference/clump_hetgam_dropout.done"
	output:
		results_clust_settings + "beds/hetgam_dropout.bed"
	params:
		indir=intermediate + "sex_depth_difference",
		index=intermediate + "subsetted_genome.fasta.fai"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cat $(find {params.indir}/ -name "*_hetgam_dropout.bed") | sed '/^$/d' | bedtools sort -g {params.index} > {output}
		"""

# 8.3 This concatenates callable regions 
# of the genome that have been phased
rule callable_regions_phased:
	input:
		bed1=results_clust_settings + "beds/callable_regions.bed",
		bed2=intermediate + "nonphased_contigs.bed"
	output:
		bed1=results_clust_settings + "beds/callable_regions_phased.bed",
		bed2=results_clust_settings + "beds/nonphased_contigs.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cp {input.bed2} {output.bed2}
		if [ -s {output.bed2} ]
		then
			bedtools subtract -a {input.bed1} -b {output.bed2} > {output.bed1}
		else
			cp {input.bed1} {output.bed1}
		fi
		"""

# 8.4 Extract sex-linked regions due to sex specfic 
# clustering and phase info. 
rule extract_sex_linked_regions:
	input:
		index=intermediate + "subsetted_genome.fasta.fai",
		bed=results_clust_settings + "beds/callable_regions_phased.bed",
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
	output:
		phase_windows=results_clust_settings + "beds/phase_windows.bed",
		unfiltered_sex_linked_bed=intermediate_clust_settings + "beds/unfiltered_sex_linked.bed",
		phase_info_bed=results_clust_settings + "beds/phase_info.bed",
		phase_sex_linked_bed=results_clust_settings + "beds/phase_sex_linked.bed",
		phase_sex_linked_tmp_bed=intermediate_clust_settings + "beds/callable_phase_sex_linked_tmp.bed",
		callable_phase_sex_linked_bed=results_clust_settings + "beds/callable_phase_sex_linked.bed",
		done=intermediate_clust_settings + "extract_sex_linked_regions_done"
	params:
		indir1=intermediate_clust_settings + "sex_linkage"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		awk '/#chrom/ {{print; exit}}' $(find {params.indir1}/ -name "*_phase_windows.bed" | head -n1) > {output.phase_windows} || :
		awk '!/#chrom/ {{print}}' $(find {params.indir1}/ -name "*_phase_windows.bed") >> {output.phase_windows} || :
		if [ $(tail -n+2 {output.phase_windows} | grep "Sex-linked" | wc -l) -gt 0 ]
		then
			grep "Sex-linked" {output.phase_windows} | cut -f 1,2,3 | bedtools sort -g {input.index} | \
			bedtools merge > {output.unfiltered_sex_linked_bed}
		else
			>> {output.unfiltered_sex_linked_bed}
		fi
		awk '/#chrom/ {{print; exit}}' $(find {params.indir1}/ -name "*_phase_info.bed" | head -n1) > {output.phase_info_bed} || :
		awk '!/#chrom/ {{print}}' $(find {params.indir1}/ -name "*_phase_info.bed") >> {output.phase_info_bed} || :
		if [ $(tail -n+2 {output.phase_info_bed} | grep "Sex-linked" | wc -l) -gt 0 ]
		then
			grep "Sex-linked" {output.phase_info_bed} | cut -f 1,2,3 | bedtools sort -g {input.index} | \
			bedtools merge > {output.phase_sex_linked_bed}
			bedtools intersect -a {output.phase_sex_linked_bed} -b {input.bed} > {output.phase_sex_linked_tmp_bed}
			if [ $(tail -n+2 {output.phase_windows} | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				bedtools intersect -a {output.unfiltered_sex_linked_bed} -b {output.phase_sex_linked_tmp_bed} > {output.callable_phase_sex_linked_bed}
			else
				>> {output.callable_phase_sex_linked_bed}
			fi
		else
			>> {output.phase_sex_linked_bed}
			>> {output.callable_phase_sex_linked_bed}
		fi
		>> {output.done}
		"""

# 8.5 Concatenate sex-linked regions with unknown
# sex-linkage due to location in border region
rule cat_sexlinked_border_variants:
	input:
		index=intermediate + "subsetted_genome.fasta.fai",
		bed=results_clust_settings + "beds/hetgam_dropout.bed",
		phase_info_bed=results_clust_settings + "beds/phase_info.bed",
		callable=results_clust_settings + "beds/callable_regions_phased.bed"
	output:
		bed=results_clust_settings + "beds/border_variants.bed"
	params:
		indir1=intermediate_clust_settings + "sex_linkage"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		awk -F'\\t' -v OFS='\\t' '{{if($4=="Unknown") print $1,$2,$3}}' {input.phase_info_bed} | \
		bedtools intersect -a - -b {input.callable} | bedtools subtract -a - -b {input.bed} > {output.bed}
		"""

# 8.6 Remove border regions and create
# the final target region bed.
rule remove_border_variants:
	input:
		bed1=results_clust_settings + "beds/callable_regions_phased.bed",
		bed2=results_clust_settings + "beds/border_variants.bed"
	output:
		bed=results_clust_settings + "beds/target_region.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		bedtools subtract -a {input.bed1} -b {input.bed2} > {output.bed}
		"""

# 8.7 Make a bedfile of the missing regions.
rule make_missing_bed:
	input:
		bed1=results_clust_settings + "beds/target_region.bed",
		index=intermediate + "subsetted_genome.fasta.fai"
	output:
		bed1=results_clust_settings + "beds/missing_region.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		bedtools complement -i {input.bed1} -g {input.index} > {output.bed1}
		"""

# 8.8 Make a bed files corresponding to each classification of region,
# i.e. sexlinked, X/Z-linked, Y/W linked, and autosomal.
# Note that 'sex_linked.bed' and 'sexshared.bed' is the same.
rule make_regions_bed:
	input:
		bed1=results_clust_settings + "beds/hetgam_dropout.bed",
		bed2=results_clust_settings + "beds/target_region.bed",
		bed3=results_clust_settings + "beds/callable_phase_sex_linked.bed",
		index=intermediate + "subsetted_genome.fasta.fai"
	output:
		bed1=results_clust_settings + "beds/target_phase_hetgam_dropout.bed",
		bed2=results_clust_settings + "beds/sex_linked.bed",
		bed3=results_clust_settings + "beds/sexshared.bed",
		bed4=results_clust_settings + "beds/sexlimited.bed",
		bed5=results_clust_settings + "beds/autosomal.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		bedtools intersect -a {input.bed1} -b {input.bed2} > {output.bed1}
		cat {output.bed1} {input.bed3} | bedtools sort -g {input.index} | bedtools merge > {output.bed2}
		cp {output.bed2} {output.bed3}
		bedtools subtract -a {output.bed2} -b {output.bed1} > {output.bed4}
		bedtools subtract -a {input.bed2} -b {output.bed2} | bedtools subtract -a - -b {input.bed1} > {output.bed5}
		"""

# 8.9 Problematic sex-linked sites due to signs of incomplete linage 
# sorting (ILS1 and ILS2)
rule sex_linked_ILS_sites:
	input:
		sexlinked_homozygotes_done=intermediate_clust_settings + "sex_linkage/sexlinked_homozygotes_done.txt",
		sex_linked_ILS2_bed_done=intermediate_clust_settings + "hetgam_vcfs/sex_linked_ILS2_bed_done.txt",
		bed1=results_clust_settings + "beds/hetgam_dropout.bed"
	output:
		bed1=results_clust_settings + "beds/sex_linked_ILS1.bed",
		bed2=results_clust_settings + "beds/sex_linked_ILS1_not_hetgam_dropout.bed",
		bed3=results_clust_settings + "beds/sex_linked_ILS2.bed"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "hetgam_vcfs",
		tmp1=intermediate_clust_settings + "beds/sex_linked_ILS1_tmp.bed",
		tmp2=intermediate_clust_settings + "beds/sex_linked_ILS1_not_hetgam_dropout_tmp.bed",
		tmp3=intermediate_clust_settings + "beds/sex_linked_ILS2_tmp.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		echo "###Header" > {params.tmp1}
		echo "###Header" > {params.tmp2}
		echo "###Header" > {params.tmp3}
		cat $(find {params.indir1}/ -name "*_sex_linked_ILS1.bed") >> {params.tmp1}
		cat $(find {params.indir2}/ -name "*_sex_linked_ILS2.bed") >> {params.tmp3}
		bedtools subtract -a {params.tmp1} -b {input.bed1} >> {params.tmp2}
		mv {params.tmp1} {output.bed1}
		mv {params.tmp2} {output.bed2}
		mv {params.tmp3} {output.bed3}
		"""

# Check if whatshap is not disabled
if config.get("whatshap") != "OFF":

	# 8.10 Save all singletons
	rule save_singletons_bed:
		input:
			index=intermediate + "subsetted_genome.fasta.fai",
			nonphased_htgm_singletons_bed=expand(intermediate_clust_settings + "beds/nonphased_singletons_{sample}_singletons.bed", sample=sample_df.index),
			done=expand(intermediate_clust_settings + "nonphased_singletons/whatshap_nonphased_{sample}_singleton_done.txt", sample=sample_df.index)
		output:
			bed=results_clust_settings + "beds/nonphased_singletons_htgm_singletons.bed",
			done=intermediate_clust_settings + "beds/save_singletons.done"
		params:
			indir=intermediate_clust_settings + "beds"
		conda:
			"../envs/bedtools.yml"
		shell:
			"""
			echo -e "#CHROM\tSTOP\tGT\tPS\tSample" > {params.indir}/nonphased_singletons_htgm_singletons_combined.bed
			grep --no-filename -v "#CHROM.*STOP" $(find {params.indir}/ -name "*_singletons.bed") >> {params.indir}/nonphased_singletons_htgm_singletons_combined.bed || :
			if [ $(grep -vc '^#' {params.indir}/nonphased_singletons_htgm_singletons_combined.bed) -gt 0 ]
			then
				bedtools sort -i {params.indir}/nonphased_singletons_htgm_singletons_combined.bed -g {input.index} | bedtools merge > {output.bed}
			else
				cp {params.indir}/nonphased_singletons_htgm_singletons_combined.bed {output.bed}
			fi
			>> {output.done}
			"""

# 8.11 Summarise all problematic sites
rule save_all_unreliable_sites_bed:
	input:
		index=intermediate + "subsetted_genome.fasta.fai",
		bed1=results_clust_settings + "beds/border_variants.bed",
		bed2=results_clust_settings + "beds/nonphased_contigs.bed",
		bed3=results_clust_settings + "beds/sex_linked_ILS1_not_hetgam_dropout.bed",
		bed4=results_clust_settings + "beds/sex_linked_ILS2.bed",
		bed5=results_clust_settings + "beds/nonphased_singletons_htgm_singletons.bed",
		done=intermediate_clust_settings + "beds/save_singletons.done"
	output:
		bed=results_clust_settings + "beds/all_unreliable_sites.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		cat {input.bed1} {input.bed2} {input.bed3} {input.bed4} {input.bed5} | cut -f1-3 | \
		bedtools sort -g {input.index} | bedtools merge > {output.bed}
		"""
