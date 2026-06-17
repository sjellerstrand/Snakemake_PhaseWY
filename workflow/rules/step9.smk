################## Step 9 ############################

# Create final vcf files. 
# 9.1 Create a list of vcf files with all the phased variants
rule merge_phased_marker_vcfs:
	input:
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt"
	output:
		merge_done=intermediate_clust_settings + "vcfs/merge_phased_marker_vcfs_done.txt",
		vcf=results_clust_settings + "vcfs/phased_all_variants.vcf.gz"
	params:
		indir=intermediate + "shapeit4/modified",
		outdir=intermediate_clust_settings,
		temp_merge=intermediate_clust_settings + "temp_merge_phased_marker",
		outdir2=results_clust_settings + "vcfs",
		max_mem=lambda wildcards, resources: f"{int(resources.mem_mb * 0.8)}M"
	resources:
		mem_mb=16000
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		mkdir -p {params.temp_merge}
		find {params.indir} -name "*_phased_all_variants.vcf.gz" | sort -n > {params.outdir}phased_all_variants_vcf_list.txt
		bcftools concat -f {params.outdir}phased_all_variants_vcf_list.txt | bgzip -c > {params.outdir}vcfs/phase_all_variants_tmp.vcf.gz
		tabix {params.outdir}vcfs/phase_all_variants_tmp.vcf.gz
		bcftools sort --max-mem {params.max_mem} -T {params.temp_merge} {params.outdir}vcfs/phase_all_variants_tmp.vcf.gz | bgzip -c > {output.vcf}
		tabix {output.vcf}
		rm -r {params.temp_merge}*
		rm {params.outdir}vcfs/phase_all_variants_tmp.vcf.gz*
		cp {input.mod_done} {output.merge_done}
		"""

# 9.2 Create autosomal dataset: 
# Merge all autosomal vcfs back together.
rule merge_autosomal_vcfs:
	input:
		autosomal_contig_vcfs_done=intermediate_clust_settings + "autosomal/sex_depth_diff_autosomal_vcfs_done.txt",
		bed1=results_clust_settings + "beds/autosomal.bed",
		unreliable=results_clust_settings + "beds/all_unreliable_sites.bed",
		vcf=config["input_vcf"]
	output:
		merge_done=intermediate_clust_settings + "vcfs/merge_autosomal_vcfs_done.txt",
		vcf1=results_clust_settings + "vcfs/autosomal.vcf.gz",
		vcf2=results_clust_settings + "vcfs/autosomal_filtered.vcf.gz"
	params:
		indir=intermediate_clust_settings + "autosomal",
		outdir=intermediate_clust_settings,
		temp_merge=intermediate_clust_settings + "temp_merge_autosomal_vcfs"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		mkdir -p {params.temp_merge}
		if [ $(find {params.indir}/ -name "*_autosomal.vcf.gz" | wc -l) -gt 0 ]
		then
			find {params.indir}/ -name "*_autosomal.vcf.gz" > {params.outdir}/autosomal_vcf_list.txt
			bcftools concat -f {params.outdir}/autosomal_vcf_list.txt -Ou | bcftools sort -T {params.temp_merge} | bgzip -c > {output.vcf1}
			tabix {output.vcf1}
			if [ $(grep -vc '^#' {input.unreliable}) -gt 0  ]
			then
				bcftools view {output.vcf1} -R {input.bed1} --regions-overlap 0 | bcftools view -T ^{input.unreliable} | bgzip -c > {output.vcf2}
			else
				bcftools view {output.vcf1} -R {input.bed1} --regions-overlap 0 | bgzip -c > {output.vcf2}
			fi
			tabix {output.vcf2}
		else
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf1}
			tabix {output.vcf1}
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf2}
			tabix {output.vcf2}
		fi
		rm -r {params.temp_merge}
		>> {output.merge_done}
		"""

# 9.3 Create both sexshared and sexlimited datasets
rule make_sexshared_and_sexlimited_datasets:
	input:
		filter_sex_depth_diff_homgam_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/filter_sex_depth_diff_sexshared_vcfs_done.txt",
		remove_hetgam_dropout_done=intermediate_clust_settings + "hetgam_vcfs/remove_hetgam_dropout_done.txt",
		bed1=results_clust_settings + "beds/sexshared.bed",
		bed2=results_clust_settings + "beds/sexlimited.bed",
		unreliable=results_clust_settings + "beds/all_unreliable_sites.bed",
		vcf=config["input_vcf"]
	output:
		merge_done=intermediate_clust_settings + "vcfs/merge_sexshared_and_sexlimited_vcfs_done.txt",
		vcf_sexshared1=results_clust_settings + "vcfs/sexshared.vcf.gz",
		vcf_sexlimited1=results_clust_settings + "vcfs/sexlimited.vcf.gz",
		vcf_sexshared2=results_clust_settings + "vcfs/sexshared_filtered.vcf.gz",
		vcf_sexlimited2=results_clust_settings + "vcfs/sexlimited_filtered.vcf.gz"
	params:
		indir=intermediate_clust_settings + "hetgam_vcfs",
		outdir=intermediate_clust_settings,
		temp_merge_dir1=intermediate_clust_settings + "temp_merge_sexshared_vcfs",
		temp_merge_dir2=intermediate_clust_settings + "temp_merge_sexlimited_vcfs"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		mkdir -p {params.temp_merge_dir1}
		mkdir -p {params.temp_merge_dir2}
		if [ $(find {params.indir}/ -name "*_sexshared.vcf.gz" | grep -v HETGAM | wc -l) -gt 0 ]
		then
			find {params.indir}/ -name "*_sexshared.vcf.gz" | grep -v HETGAM > {params.outdir}/sexshared_vcf_list.txt
			bcftools concat -a -f {params.outdir}/sexshared_vcf_list.txt -Ou | bcftools sort -T {params.temp_merge_dir1} | bcftools norm --rm-dup all | bgzip -c > {output.vcf_sexshared1}
			tabix {output.vcf_sexshared1}
			if [ $(grep -vc '^#' {input.unreliable}) -gt 0  ]
			then
				bcftools view {output.vcf_sexshared1} -R {input.bed1} --regions-overlap 0 | bcftools view -T ^{input.unreliable} | bgzip -c > {output.vcf_sexshared2}
			else
				bcftools view {output.vcf_sexshared1} -R {input.bed1} --regions-overlap 0 | bgzip -c > {output.vcf_sexshared2}
			fi
			tabix {output.vcf_sexshared2}
		else
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf_sexshared1}
			tabix {output.vcf_sexshared1}
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf_sexshared2}
			tabix {output.vcf_sexshared2}
		fi
		if [ $(find {params.indir}/ -name "*_sexlimited.vcf.gz" | wc -l) -gt 0 ]
		then
			find {params.indir}/ -name "*_sexlimited.vcf.gz" > {params.outdir}/sexlimited_vcf_list.txt
			bcftools concat -f {params.outdir}/sexlimited_vcf_list.txt -Ou | bcftools sort -T {params.temp_merge_dir2} | bgzip -c > {output.vcf_sexlimited1}
			tabix {output.vcf_sexlimited1}
			if [ $(grep -vc '^#' {input.unreliable}) -gt 0  ]
			then
				bcftools view {output.vcf_sexlimited1} -R {input.bed2} --regions-overlap 0 | bcftools view -T ^{input.unreliable} | bgzip -c > {output.vcf_sexlimited2}
			else
				bcftools view {output.vcf_sexlimited1} -R {input.bed2} --regions-overlap 0 | bgzip -c > {output.vcf_sexlimited2}
			fi
			tabix {output.vcf_sexlimited2}
		else
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf_sexlimited1}
			tabix {output.vcf_sexlimited1}
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf_sexlimited2}
			tabix {output.vcf_sexlimited2}
		fi
		cat {input.filter_sex_depth_diff_homgam_vcfs_done} {input.remove_hetgam_dropout_done} >> {output.merge_done}
		"""

# 9.4 Save sites from contigs that could not be phased.
rule save_nonphased_sites:
	input:
		vcf=config["input_vcf"],
		callable=results_clust_settings + "beds/callable_regions.bed",
		nonphased=intermediate + "nonphased_contigs.bed"
	output:
		done=intermediate_clust_settings + "vcfs/nonphased_variants_done.txt",
		vcf=results_clust_settings + "vcfs/nonphased_variants.vcf.gz"
	params:
		min_dp=config["min_dp"],
		missing=config["missing"],
		min_mean=config["min_mean"],
		max_mean=config["max_mean"]
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		if [ $(grep -vc '^#' {input.nonphased}) -gt 0 ]
		then
			bcftools view {input.vcf} -R {input.nonphased} | vcftools --vcf - --bed {input.callable} --minDP {params.min_dp} --max-missing {params.missing} --min-meanDP {params.min_mean} --max-meanDP {params.max_mean} \
			--recode --recode-INFO-all --stdout | bgzip -c > {output.vcf}
		else
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf}
		fi
		tabix {output.vcf}
		>> {output.done}
		"""

# 9.5 Save all problematic sex-linked sites removed
# due to location in border region
rule save_problematic_sex_linked_sites:
	input:
		vcf=config["input_vcf"],
		callable=results_clust_settings + "beds/callable_regions.bed",
		border_variants=results_clust_settings + "beds/border_variants.bed"
	output:
		vcf=results_clust_settings + "vcfs/border_variants.vcf.gz",
		done=intermediate_clust_settings + "vcfs/border_variants_done.txt"
	params:
		min_dp=config["min_dp"],
		missing=config["missing"],
		min_mean=config["min_mean"],
		max_mean=config["max_mean"],
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		if [ $(grep -vc '^#' {input.border_variants}) -gt 0 ]
		then
			bcftools view {input.vcf} -R {input.border_variants} | vcftools --vcf - --bed {input.callable} --minDP {params.min_dp} --max-missing {params.missing} --min-meanDP {params.min_mean} --max-meanDP {params.max_mean} \
			--recode --recode-INFO-all --stdout | bgzip -c > {output.vcf}
		else
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf}
		fi
		tabix {output.vcf}
		>> {output.done}
		"""

# 9.6 Save all problematic sex-linked sites (ILS1 and ILS2).
rule save_problem_sex_linked_incomplete_sorting:
	input:
		vcf=config["input_vcf"],
		merge_done=intermediate_clust_settings + "vcfs/merge_sexshared_and_sexlimited_vcfs_done.txt",
		vcf_sexshared=results_clust_settings + "vcfs/sexshared.vcf.gz",
		vcf_sexlimited=results_clust_settings + "vcfs/sexlimited.vcf.gz",
		bed1=results_clust_settings + "beds/sex_linked_ILS1_not_hetgam_dropout.bed",
		bed2=results_clust_settings + "beds/sex_linked_ILS2.bed"
	output:
		vcf1=results_clust_settings + "vcfs/all_inds_ILS1_not_hetgam_dropout.vcf.gz",
		vcf2=results_clust_settings + "vcfs/sexshared_ILS1_not_hetgam_dropout.vcf.gz",
		vcf3=results_clust_settings + "vcfs/sexlimited_ILS1_not_hetgam_dropout.vcf.gz",
		vcf4=results_clust_settings + "vcfs/all_inds_ILS2.vcf.gz",
		vcf5=results_clust_settings + "vcfs/sexshared_ILS2.vcf.gz",
		vcf6=results_clust_settings + "vcfs/sexlimited_ILS2.vcf.gz",
		done=intermediate_clust_settings + "vcfs/ILS_done.txt"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		if [ $(grep -vc '^#' {input.bed1}) -gt 0  ]
		then
			bcftools view {input.vcf} -R {input.bed1} | bgzip -c > {output.vcf1}
			bcftools view {input.vcf_sexshared} -R {input.bed1} | bgzip -c > {output.vcf2}
			bcftools view {input.vcf_sexlimited} -R {input.bed1} | bgzip -c > {output.vcf3}
		else
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf1}
			bcftools view {input.vcf_sexshared} -h | bgzip -c > {output.vcf2}
			bcftools view {input.vcf_sexlimited} -h | bgzip -c > {output.vcf3}
		fi
		tabix {output.vcf1}
		tabix {output.vcf2}
		tabix {output.vcf3}
		if [ $(grep -vc '^#' {input.bed2}) -gt 0  ]
		then
			bcftools view {input.vcf} -R {input.bed2} | bgzip -c > {output.vcf4}
			bcftools view {input.vcf_sexshared} -R {input.bed2} | bgzip -c > {output.vcf5}
			bcftools view {input.vcf_sexlimited} -R {input.bed2} | bgzip -c > {output.vcf6}
		else
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf4}
			bcftools view {input.vcf_sexshared} -h | bgzip -c > {output.vcf5}
			bcftools view {input.vcf_sexlimited} -h | bgzip -c > {output.vcf6}
		fi
		tabix {output.vcf4}
		tabix {output.vcf5}
		tabix {output.vcf6}
		>> {output.done}
		"""

# 9.7 Save all singletons
rule save_singletons:
	input:
		nonphased_htgm_singletons_bed=results_clust_settings + "beds/nonphased_singletons_htgm_singletons.bed",
		vcf=config["input_vcf"],
		done=intermediate_clust_settings + "beds/save_singletons.done"
	output:
		vcf=results_clust_settings + "vcfs/heterogametic_sexlinked_singletons.vcf.gz",
		done=intermediate_clust_settings + "vcfs/singletons_done.txt"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		if [ $(grep -vc '^#' {input.nonphased_htgm_singletons_bed}) -gt 0 ]
		then
			bcftools view {input.vcf} -R {input.nonphased_htgm_singletons_bed} | bgzip -c > {output.vcf}
		else
			bcftools view {input.vcf} -h | bgzip -c > {output.vcf}
		fi
		tabix {output.vcf}
		>> {output.done}
		"""

# 9.8 Summarise all problematic sites
rule save_all_unreliable_sites:
	input:
		bed=results_clust_settings + "beds/all_unreliable_sites.bed",
		vcf1=results_clust_settings + "vcfs/nonphased_variants.vcf.gz",
		vcf2=results_clust_settings + "vcfs/border_variants.vcf.gz",
		vcf3=results_clust_settings + "vcfs/all_inds_ILS1_not_hetgam_dropout.vcf.gz",
		vcf4=results_clust_settings + "vcfs/all_inds_ILS2.vcf.gz",
		vcf5=results_clust_settings + "vcfs/heterogametic_sexlinked_singletons.vcf.gz",
		done1=intermediate_clust_settings + "vcfs/nonphased_variants_done.txt",
		done3=intermediate_clust_settings + "vcfs/border_variants_done.txt",
		done4=intermediate_clust_settings + "vcfs/ILS_done.txt",
		done5=intermediate_clust_settings + "vcfs/singletons_done.txt"
	output:
		done=intermediate_clust_settings + "vcfs/all_unreliable_sites_done.txt",
		vcf=results_clust_settings + "vcfs/all_unreliable_sites.vcf.gz"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		bcftools concat -a {input.vcf1} {input.vcf2} {input.vcf3} {input.vcf4} {input.vcf5} | \
		bcftools norm --rm-dup all | \
		bgzip -c > {output.vcf}
		tabix {output.vcf}
		>> {output.done}
		"""
