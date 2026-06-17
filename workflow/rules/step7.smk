########################## Step 7 ##########################

# 7.1 Split phased contig with sex-linked regions identified 
# based on phased genotypes into autosomal and 
# sex-linked regions, and sexshared and sexlimited genomic 
# regions.

# 7.1.1 Combine sex-linked sites identified from phased genotypes
# (*phase_windows.bed) and from depth differences between sexes 
# (*hetgam_dropout.bed). 
rule sex_link_hetgam_dropout_bed:
	input:
		sexlinked_bed_done=intermediate_clust_settings + "sex_linkage/sexlinked_bed_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		index=intermediate + "subsetted_genome.fasta.fai",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sex_link_hetgam_dropout_bed_done=intermediate_clust_settings + "sex_link_hetgam_dropout/sex_link_hetgam_dropout_bed_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate + "sex_depth_difference",
		outdir=intermediate_clust_settings + "sex_link_hetgam_dropout"
	conda: "../envs/bedtools.yml"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				if [[ -s {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
					then
					mkdir -p {params.outdir}/"$contig"
					# Create a sex-linked bed file based on sex-depth differences
					cat {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed \
					{params.indir1}/"$contig"/"$contig"_sex_linked.bed | \
					bedtools sort -g {input.index} | bedtools merge \
					> {params.outdir}/"$contig"/"$contig"_sex_link_hetgam_dropout.bed 
				fi
			fi
		done < {input.filt_list}
		>> {output.sex_link_hetgam_dropout_bed_done}
		"""

# 7.1.2 Extract sex-linked sites from the phased vcf files
# for homgam individuals
rule sex_link_homgam_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sex_link_hetgam_dropout_bed_done=intermediate_clust_settings + "sex_link_hetgam_dropout/sex_link_hetgam_dropout_bed_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		sex_link_homgam_vcfs_done=intermediate_clust_settings + "homgam_vcfs/sex_link_homgam_vcfs_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate + "sex_depth_difference",
		indir3=intermediate + "shapeit4/modified",
		indir4=intermediate_clust_settings + "sex_link_hetgam_dropout",
		outdir=intermediate_clust_settings + "homgam_vcfs"
	conda:
		"../envs/bcftools.yml"
	shell:
		"""
		HOMGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HOMGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
				then
				if [[ -s {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
					then
					bcftools view {params.indir3}/"$contig"_phased_all_variants.vcf.gz \
					-s $HOMGAM -R {params.indir4}/"$contig"/"$contig"_sex_link_hetgam_dropout.bed | \
					bgzip -c > {params.outdir}/"$contig"_HOMGAM.vcf.gz;
				else
					bcftools view {params.indir3}/"$contig"_phased_all_variants.vcf.gz \
					-s $HOMGAM -R {params.indir1}/"$contig"/"$contig"_sex_linked.bed | \
					bgzip -c > {params.outdir}/"$contig"_HOMGAM.vcf.gz
				fi;
				tabix {params.outdir}/"$contig"_HOMGAM.vcf.gz
			fi
		done < {input.filt_list}
		>> {output.sex_link_homgam_vcfs_done}
		"""

# 7.1.3 Extract sex-linked sites from the phased vcf files
# for hetgam individuals
rule sex_link_hetgam_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sex_link_hetgam_dropout_bed_done=intermediate_clust_settings + "sex_link_hetgam_dropout/sex_link_hetgam_dropout_bed_done.txt",
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		sex_link_hetgam_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/sex_link_hetgam_vcfs_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate + "shapeit4/modified",
		indir3=intermediate_clust_settings + "sex_link_hetgam_dropout",
		outdir=intermediate_clust_settings + "hetgam_vcfs"
	conda:
		"../envs/bcftools.yml"
	shell:
		"""
		HETGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
		while read line
		do
			contig=$(echo $line | cut -f1)
			mkdir -p {params.outdir}/"$contig"
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				bcftools view {params.indir2}/"$contig"_phased_all_variants.vcf.gz \
				-s $HETGAM -R {params.indir1}/"$contig"/"$contig"_sex_linked.bed | \
				bgzip -c > {params.outdir}/"$contig"/"$contig"_HETGAM.vcf.gz
			fi
		done < {input.filt_list}
		>> {output.sex_link_hetgam_vcfs_done}
		"""

# 7.2 Split heterogamets into sexlimited and sexshared genotypes
rule split_hetgam_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sex_link_hetgam_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/sex_link_hetgam_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		split_hetgam_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/split_hetgam_vcfs_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		split_phase="workflow/scripts/split_phase.py",
		outdir=intermediate_clust_settings + "hetgam_vcfs"
	conda:
		"../envs/python_bcftools.yml"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				python3 {params.split_phase} -i {params.outdir}/"$contig"/"$contig"_HETGAM.vcf.gz -p left -o {params.outdir}/"$contig"/"$contig"_HETGAM_left.vcf
				bgzip -f {params.outdir}/"$contig"/"$contig"_HETGAM_left.vcf
				tabix {params.outdir}/"$contig"/"$contig"_HETGAM_left.vcf.gz
				python3 {params.split_phase} -i {params.outdir}/"$contig"/"$contig"_HETGAM.vcf.gz -p right -o {params.outdir}/"$contig"/"$contig"_HETGAM_right.vcf
				bgzip -f {params.outdir}/"$contig"/"$contig"_HETGAM_right.vcf
				tabix {params.outdir}/"$contig"/"$contig"_HETGAM_right.vcf.gz
			fi
		done < {input.filt_list}
		>> {output.split_hetgam_vcfs_done}
		"""

# 7.2.1 Concatenate individual heterogametes into sexlimited and sexshared phase (part 1)
rule hetgam_sexlimited_vcf:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		split_hetgam_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/split_hetgam_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		hetgam_sexlimited_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/all_HETGAM_{hetgam_sample}_sexlimited_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "hetgam_vcfs"
	conda:
		"../envs/bcftools.yml"
	shell:
		"""
		mkdir -p {params.indir2}/{wildcards.hetgam_sample}_sexlimited
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				mkdir -p {params.indir2}/{wildcards.hetgam_sample}_sexlimited/"$contig"_temp_merge

				if [ $(cat {params.indir1}/"$contig"/"$contig"_{wildcards.hetgam_sample}_het_left.bed | grep -v "#" | awk '{{ print; exit }}' | wc -l) -gt 0 ]
				then
					###### Extract sexlimited sites with left phase
					bcftools view {params.indir2}/"$contig"/"$contig"_HETGAM_left.vcf.gz -s {wildcards.hetgam_sample} \
					-R {params.indir1}/"$contig"/"$contig"_{wildcards.hetgam_sample}_het_left.bed | \
					bgzip -c > {params.indir2}/{wildcards.hetgam_sample}_sexlimited/"$contig"_temp_merge/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexlimited1.vcf.gz
				fi

				if [ $(cat {params.indir1}/"$contig"/"$contig"_{wildcards.hetgam_sample}_het_right.bed | grep -v "#" | awk '{{ print; exit }}' | wc -l) -gt 0 ]
				then
					###### Extract sexlimited sites with right phase
					bcftools view {params.indir2}/"$contig"/"$contig"_HETGAM_right.vcf.gz -s {wildcards.hetgam_sample} \
					-R {params.indir1}/"$contig"/"$contig"_{wildcards.hetgam_sample}_het_right.bed  | \
					bgzip -c > {params.indir2}/{wildcards.hetgam_sample}_sexlimited/"$contig"_temp_merge/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexlimited2.vcf.gz
				fi

				###### Concatenate sexlimited genotypes
				if [ $(find {params.indir2}/{wildcards.hetgam_sample}_sexlimited/"$contig"_temp_merge/ -name ""$contig"_HETGAM_{wildcards.hetgam_sample}_sexlimited?.vcf.gz" | wc -l) -gt 0 ]
				then
					concat_files=$(find {params.indir2}/{wildcards.hetgam_sample}_sexlimited/"$contig"_temp_merge/ -name ""$contig"_HETGAM_{wildcards.hetgam_sample}_sexlimited?.vcf.gz" | tr '\\n' ' ')
					bcftools concat $concat_files -Ou | \
					bcftools sort -T {params.indir2}/{wildcards.hetgam_sample}_sexlimited/"$contig"_temp_merge/ | \
					bgzip -c > {params.indir2}/"$contig"/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexlim.vcf.gz
					tabix {params.indir2}/"$contig"/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexlim.vcf.gz
				else
					bcftools view {params.indir2}/"$contig"/"$contig"_HETGAM_right.vcf.gz -s {wildcards.hetgam_sample} -h | \
					bgzip -c > {params.indir2}/"$contig"/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexlim.vcf.gz
					tabix {params.indir2}/"$contig"/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexlim.vcf.gz
				fi
				rm -r {params.indir2}/{wildcards.hetgam_sample}_sexlimited/"$contig"_temp_merge
			fi
		done < {input.filt_list}
		rm -r {params.indir2}/{wildcards.hetgam_sample}_sexlimited
		>> {output.hetgam_sexlimited_vcfs_done}
		"""

# 7.2.2 Concatenate individual heterogametes into sexlimited and sexshared phase (part 2)
rule hetgam_sexshared_vcf:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		split_hetgam_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/split_hetgam_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		hetgam_sexshared_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/all_HETGAM_{hetgam_sample}_sexshared_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "hetgam_vcfs"
	conda:
		"../envs/bcftools.yml"
	shell:
		"""
		mkdir -p {params.indir2}/{wildcards.hetgam_sample}_sexshared
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				mkdir -p {params.indir2}/{wildcards.hetgam_sample}_sexshared/"$contig"_temp_merge

				if [ $(cat {params.indir1}/"$contig"/"$contig"_{wildcards.hetgam_sample}_het_left.bed | grep -v "#" | awk '{{ print; exit }}' | wc -l) -gt 0 ]
				then
					###### Extract sexshared sites with right phase
					bcftools view {params.indir2}/"$contig"/"$contig"_HETGAM_right.vcf.gz -s {wildcards.hetgam_sample} \
					-R {params.indir1}/"$contig"/"$contig"_{wildcards.hetgam_sample}_het_left.bed | \
					bgzip -c > {params.indir2}/{wildcards.hetgam_sample}_sexshared/"$contig"_temp_merge/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexshared1.vcf.gz
				fi

				if [ $(cat {params.indir1}/"$contig"/"$contig"_{wildcards.hetgam_sample}_het_right.bed | grep -v "#" | awk '{{ print; exit }}' | wc -l) -gt 0 ]
				then
					###### Extract sexshared sites with left phase
					bcftools view {params.indir2}/"$contig"/"$contig"_HETGAM_left.vcf.gz -s {wildcards.hetgam_sample} \
					-R {params.indir1}/"$contig"/"$contig"_{wildcards.hetgam_sample}_het_right.bed  | \
					bgzip -c > {params.indir2}/{wildcards.hetgam_sample}_sexshared/"$contig"_temp_merge/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexshared2.vcf.gz
				fi

				###### Concatenate sexshared genotypes
				if [ $(find {params.indir2}/{wildcards.hetgam_sample}_sexshared/"$contig"_temp_merge/ -name ""$contig"_HETGAM_{wildcards.hetgam_sample}_sexshared?.vcf.gz" | wc -l) -gt 0 ]
				then
					concat_files=$(find {params.indir2}/{wildcards.hetgam_sample}_sexshared/"$contig"_temp_merge/ -name ""$contig"_HETGAM_{wildcards.hetgam_sample}_sexshared?.vcf.gz" | tr '\\n' ' ')
					bcftools concat $concat_files -Ou | \
					bcftools sort -T {params.indir2}/{wildcards.hetgam_sample}_sexshared/"$contig"_temp_merge/ | \
					bgzip -c > {params.indir2}/"$contig"/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexshared.vcf.gz
					tabix {params.indir2}/"$contig"/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexshared.vcf.gz
				else
					bcftools view {params.indir2}/"$contig"/"$contig"_HETGAM_left.vcf.gz -s {wildcards.hetgam_sample} -h | \
					bgzip -c > {params.indir2}/"$contig"/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexshared.vcf.gz
					tabix {params.indir2}/"$contig"/"$contig"_HETGAM_{wildcards.hetgam_sample}_sexshared.vcf.gz
				fi
				rm -r {params.indir2}/{wildcards.hetgam_sample}_sexshared/"$contig"_temp_merge
			fi
		done < {input.filt_list}
		rm -r {params.indir2}/{wildcards.hetgam_sample}_sexshared
		>> {output.hetgam_sexshared_vcfs_done}
		"""

# 7.2.3 Merge heterogametic haplotypes
rule merge_sexlimited_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlimited_vcfs_done=lambda wildcards: expand(intermediate_clust_settings + "hetgam_vcfs/all_HETGAM_{hetgam_sample}_sexlimited_done.txt", hetgam_sample=hetgam_samples),
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		merge_sexlimited_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/merge_sexlimited_vcfs_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "hetgam_vcfs"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		INDS_ORDER_HETS=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				merge_files=$(ls {params.indir2}/"$contig"/"$contig"_HETGAM_*_sexlim.vcf.gz);
				bcftools merge $merge_files --force-single -Ou | bcftools view -s $INDS_ORDER_HETS | \
				vcffixup - | \
				bgzip -c > {params.indir2}/"$contig"/"$contig"_sexlimited_temp1.vcf.gz;
			fi
		done < {input.filt_list}
		>> {output.merge_sexlimited_vcfs_done}
		"""

# 7.2.4 If there were heterogametic drop-out regions, remove these in the sexlimited file
rule remove_hetgam_dropout:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		merge_sexlimited_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/merge_sexlimited_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		remove_hetgam_dropout_done=intermediate_clust_settings + "hetgam_vcfs/remove_hetgam_dropout_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate + "sex_depth_difference",
		indir3=intermediate_clust_settings + "hetgam_vcfs",
		outdir=intermediate_clust_settings + "hetgam_vcfs"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				if [[ -s {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					vcftools --gzvcf {params.indir3}/"$contig"/"$contig"_sexlimited_temp1.vcf.gz \
					--exclude-bed {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed \
					--recode --recode-INFO-all --stdout | \
					bgzip -c > {params.outdir}/"$contig"/"$contig"_sexlimited.vcf.gz;
				else
					cp {params.indir3}/"$contig"/"$contig"_sexlimited_temp1.vcf.gz {params.outdir}/"$contig"/"$contig"_sexlimited.vcf.gz;
				fi
				tabix {params.outdir}/"$contig"/"$contig"_sexlimited.vcf.gz;
			fi
		done < {input.filt_list}
		>> {output.remove_hetgam_dropout_done}
		"""

# 7.3.1 Merge sexshared haplotypes
rule merge_sexshared_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexshared_vcfs_done=lambda wildcards: expand(intermediate_clust_settings + "hetgam_vcfs/all_HETGAM_{hetgam_sample}_sexshared_done.txt", hetgam_sample=hetgam_samples),
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		merge_sexshared_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/merge_sexshared_vcfs_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "hetgam_vcfs"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		INDS_ORDER_HETS=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)	
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				merge_files=$(ls {params.indir2}/"$contig"/"$contig"_HETGAM_*_sexshared.vcf.gz);
				bcftools merge $merge_files --force-single -Ou | \
				bcftools view -s $INDS_ORDER_HETS | \
				bgzip -c > {params.indir2}/"$contig"/"$contig"_sexshared_temp1.vcf.gz;
				tabix {params.indir2}/"$contig"/"$contig"_sexshared_temp1.vcf.gz;
			fi
		done < {input.filt_list}
		>> {output.merge_sexshared_vcfs_done}
		"""

# 7.3.2 If there were heterogametic drop-out variants, 
# add these to the sexshared file for the heterogametic sex
rule hetgam_sexshared_dropout_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt",
		sex_link_hetgam_dropout_bed_done=intermediate_clust_settings + "sex_link_hetgam_dropout/sex_link_hetgam_dropout_bed_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		hetgam_sexshared_dropout_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/hetgam_sexshared_dropout_vcfs_done.txt"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate + "sex_depth_difference",
		indir3=intermediate + "shapeit4/modified",
		indir4=intermediate_clust_settings + "sex_link_hetgam_dropout",
		outdir=intermediate_clust_settings + "hetgam_vcfs",
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		HETGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1);
		while read line
		do
			contig=$(echo $line | cut -f1)
			mkdir -p {params.outdir}/"$contig"
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				if [[ -s {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					bcftools view {params.indir3}/"$contig"_phased_all_variants.vcf.gz \
					-s $HETGAM -R {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed | \
					vcftools --vcf - \
					--exclude-bed {params.indir1}/"$contig"/"$contig"_sex_linked.bed \
					--recode --recode-INFO-all --stdout | \
					bgzip -c > {params.outdir}/"$contig"/"$contig"_sexshared_temp2.vcf.gz
				fi
			fi
		done < {input.filt_list}
		>> {output.hetgam_sexshared_dropout_vcfs_done}
		"""

# 7.3.3 Re-code homozygose genotypes to haploid for the heterogametic sex
rule sex_link_diploid2haploid:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		hetgam_sexshared_dropout_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/hetgam_sexshared_dropout_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sex_link_diploid2haploid_done=intermediate_clust_settings + "hetgam_vcfs/sex_link_diploid2haploid_done.txt"
	conda:
		"../envs/python_bcftools.yml"
	params:
		diploid2haploid="workflow/scripts/diploid2haploid.py",
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate + "sex_depth_difference",
		indir3=intermediate_clust_settings + "hetgam_vcfs"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				if [[ -s {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					if [ $(bcftools view -H {params.indir3}/"$contig"/"$contig"_sexshared_temp2.vcf.gz | awk '{{ print; exit }}' | wc -l) -gt 0 ]
					then
						###### Re-code homozygose genotypes to haploid for the heterogametic sex
						python3 {params.diploid2haploid} -i {params.indir3}/"$contig"/"$contig"_sexshared_temp2.vcf.gz \
						-o {params.indir3}/"$contig"/"$contig"_HETGAM_sexshared.vcf
						bgzip -f {params.indir3}/"$contig"/"$contig"_HETGAM_sexshared.vcf
						tabix {params.indir3}/"$contig"/"$contig"_HETGAM_sexshared.vcf.gz
					fi
				fi
			fi
		done < {input.filt_list}
		>> {output.sex_link_diploid2haploid_done}
		"""

# 7.3.4 Concatenate sexshared sites for the heterogametic sex
rule concatenate_hom_het_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		sex_link_diploid2haploid_done=intermediate_clust_settings + "hetgam_vcfs/sex_link_diploid2haploid_done.txt",
		merge_sexshared_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/merge_sexshared_vcfs_done.txt",
		hetgam_sexshared_dropout_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/hetgam_sexshared_dropout_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		concatenate_hom_het_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/concatenate_hom_het_vcfs_done.txt"
	conda:
		"../envs/bcftools.yml"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate + "sex_depth_difference",
		indir3=intermediate_clust_settings + "hetgam_vcfs"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				if [[ -s {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					if [ $(bcftools view -H {params.indir3}/"$contig"/"$contig"_sexshared_temp2.vcf.gz | awk '{{ print; exit }}' | wc -l) -gt 0 ]
					then
						##### Concatenate sexshared sites for the heterogametic sex
						bcftools concat {params.indir3}/"$contig"/"$contig"_sexshared_temp1.vcf.gz \
						{params.indir3}/"$contig"/"$contig"_HETGAM_sexshared.vcf.gz -Ou | \
						bcftools sort -T {params.indir3}/"$contig"/"$contig"_temp_merge | \
						bgzip -c > {params.indir3}/"$contig"/"$contig"_sexshared_temp3.vcf.gz
					else
						cp {params.indir3}/"$contig"/"$contig"_sexshared_temp1.vcf.gz \
						{params.indir3}/"$contig"/"$contig"_sexshared_temp3.vcf.gz
					fi
				else
					cp {params.indir3}/"$contig"/"$contig"_sexshared_temp1.vcf.gz \
					{params.indir3}/"$contig"/"$contig"_sexshared_temp3.vcf.gz
				fi
				tabix {params.indir3}/"$contig"/"$contig"_sexshared_temp3.vcf.gz
			fi
		done < {input.filt_list}
		>> {output.concatenate_hom_het_vcfs_done}
		"""

# 7.3.5 Merge sexshared haplotypes from the homogametic and heterogametic sex
rule merge_sex_link_sexshared_hom_het_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sex_link_homgam_vcfs_done=intermediate_clust_settings + "homgam_vcfs/sex_link_homgam_vcfs_done.txt",
		concatenate_hom_het_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/concatenate_hom_het_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		merge_sex_link_sexshared_hom_het_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/merge_sex_link_sexshared_hom_het_vcfs_done.txt"
	conda:
		"../envs/bcftools.yml"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "homgam_vcfs",
		indir3=intermediate_clust_settings + "hetgam_vcfs"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				bcftools merge {params.indir2}/"$contig"_HOMGAM.vcf.gz \
				{params.indir3}/"$contig"/"$contig"_sexshared_temp3.vcf.gz | \
				bgzip -c > {params.indir3}/"$contig"/"$contig"_sexshared_temp4.vcf.gz
				tabix {params.indir3}/"$contig"/"$contig"_sexshared_temp4.vcf.gz
			fi
		done < {input.filt_list}
		>> {output.merge_sex_link_sexshared_hom_het_vcfs_done}
		"""

# 7.3.6 Order individuals and filter sites
rule filter_sex_link_sexshared_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		merge_sex_link_sexshared_hom_het_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/merge_sex_link_sexshared_hom_het_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		filter_sex_link_sexshared_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/filter_sex_link_sexshared_vcfs_done.txt"
	conda:
		"../envs/bcfvcftools.yml"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "hetgam_vcfs"
	shell:
		"""
		INDS_ORDER=$(cat {input.sample_table} | awk -F'\\t' '$1!="sample_name" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				bcftools view {params.indir2}/"$contig"/"$contig"_sexshared_temp4.vcf.gz -s $INDS_ORDER | \
				vcffixup - | \
				bgzip -c > {params.indir2}/"$contig"/"$contig"_sexshared.vcf.gz
				tabix {params.indir2}/"$contig"/"$contig"_sexshared.vcf.gz
			fi
		done < {input.filt_list}
		>> {output.filter_sex_link_sexshared_vcfs_done}
		"""

# 7.4.1 Identify problematic sex-linked sites
rule overlaps_het_sexshared_bed:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		remove_hetgam_dropout_done=intermediate_clust_settings + "hetgam_vcfs/remove_hetgam_dropout_done.txt",
		filter_sex_link_sexshared_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/filter_sex_link_sexshared_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		overlaps_het_sexshared_bed_done=intermediate_clust_settings + "hetgam_vcfs/overlaps_het_sexshared_bed_done.txt"
	conda:
		"../envs/bedtools.yml"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "hetgam_vcfs"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				bedtools intersect -a {params.indir2}/"$contig"/"$contig"_sexlimited.vcf.gz \
				-b {params.indir2}/"$contig"/"$contig"_sexshared.vcf.gz | \
				awk 'NR==1 {{
				info1_col=0; info2_col=0; gt_count=0;
				for(i=1;i<=NF;i++){{
				if($i == "GT" || $i ~ /^GT(:|$)/){{
				gt_count++;
				if(gt_count==1) info1_col=i-1; \
				else if(gt_count==2) info2_col=i-1;
				}}
				}}
				split($(info1_col), arr1, ";"); af1_index=0; \
				for(i=1;i<=length(arr1);i++) if(arr1[i] ~ /^AF=/) af1_index=i;
				split($(info2_col), arr2, ";"); af2_index=0;
				for(i=1;i<=length(arr2);i++) if(arr2[i] ~ /^AF=/) af2_index=i;
				}}
				{{
				split($(info1_col), arr1, ";");
				split($(info2_col), arr2, ";");
				af1 = (af1_index>0) ? substr(arr1[af1_index],4) : "NA";
				af2 = (af2_index>0) ? substr(arr2[af2_index],4) : "NA";
				if(af1 != "NA" && af2 != "NA" && af1 > 0 && af1 < 1 && af2 > 0 && af2 < 1){{
				print $1"\t"$2-1"\t"$2;
				}} 
				}}'	> {params.indir2}/"$contig"/"$contig"_overlaps.bed
			fi
		done < {input.filt_list}
		>> {output.overlaps_het_sexshared_bed_done}
		"""

# 7.4.2 Check if both alleles occur on both sex chromosomes.
# This is likely a sign of unsuccessful phasing or incomplete lineage sorting
rule sex_linked_ILS2_bed:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		overlaps_het_sexshared_bed_done=intermediate_clust_settings + "hetgam_vcfs/overlaps_het_sexshared_bed_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sex_linked_ILS2_bed_done=intermediate_clust_settings + "hetgam_vcfs/sex_linked_ILS2_bed_done.txt"
	conda:
		"../envs/bcftools.yml"
	params:
		indir1=intermediate_clust_settings + "sex_linkage",
		indir2=intermediate_clust_settings + "hetgam_vcfs"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir1}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -gt 0 ]
			then
				if [ $(cat {params.indir2}/"$contig"/"$contig"_overlaps.bed | awk 'NR>0 {{ print; exit }}' | wc -l) -gt 0 ]
				then
					bcftools query {params.indir2}/"$contig"/"$contig"_sexlimited.vcf.gz \
					-R {params.indir2}/"$contig"/"$contig"_overlaps.bed -f '%CHROM\\t%POS\\t%INFO/AN\\t%INFO/AC\\n' | \
					awk -F'\\t|,' '{{ sum=0; for(i=4; i==NF; ++i) sum+=$i }} \
					{{ if((sum==0 || $4/$3==1)) print $1"\\t"$2"\\tmonomorphic\\t"; \
					else if(NF>4) print $1"\\t"$2"\\tmultiallelic"; \
					else print $1"\\t"$2"\\tpolymorphic" }}' | \
					sort -k2,2 > {params.indir2}/"$contig"/"$contig"_sexlimited_var.txt;

					bcftools query {params.indir2}/"$contig"/"$contig"_sexshared.vcf.gz \
					-R {params.indir2}/"$contig"/"$contig"_overlaps.bed -f '%CHROM\\t%POS\\t%INFO/AN\\t%INFO/AC\\n' | \
					awk -F'\\\t|,' '{{ sum=0; for(i=4; i==NF; ++i) sum+=$i }} \
					{{ if((sum==0 || $4/$3==1)) print $1"\\t"$2"\\tmonomorphic\\t"; \
					else if(NF>4) print $1"\\t"$2"\\tmultiallelic"; \
					else print $1"\\t"$2"\\tpolymorphic" }}' | \
					sort -k2,2 > {params.indir2}/"$contig"/"$contig"_sexshared_var.txt;

					join -j 2 -o 1.1,1.2,1.3,2.3 {params.indir2}/"$contig"/"$contig"_sexlimited_var.txt {params.indir2}/"$contig"/"$contig"_sexshared_var.txt | \
					awk '{{ if(($3 == "polymorphic" && $3 == $4) || \
					($3 == "multiallelic" && $4 != "monomorphic") || \
					($4 == "multiallelic" && $3 != "monomorphic")) \
					print $1"\\t"$2-1"\\t"$2"\\t"$3"\\t"$4 }}' \
					> {params.indir2}/"$contig"/"$contig"_sex_linked_ILS2.bed;
				fi
			fi
		done < {input.filt_list}
		>> {output.sex_linked_ILS2_bed_done}
		"""

# 7.5 Split phased contig with sex-linked regions identified 
# due to sex depth differences into sex-linked regions,
# and sexshared and sexlimited 
# genomic regions.

# 7.5.1 Split sex-linked regions into heterogametes and homogametes (part 1)
rule sex_depth_diff_sexshared_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		sex_depth_diff_sexshared_vcfs_done=intermediate_clust_settings + "temp_sex_link/sex_depth_diff_sexshared_vcfs_done.txt"
	params:
		indir1=intermediate + "sex_depth_difference",
		indir2=intermediate_clust_settings + "sex_linkage",
		indir3=intermediate + "shapeit4/modified",
		outdir=intermediate_clust_settings + "temp_sex_link",
	conda:
		"../envs/bcftools.yml"
	shell:
		"""
		HOMGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HOMGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)			
		while read line
		do
			contig=$(echo $line | cut -f1)
			mkdir -p {params.outdir}/"$contig"
			if [ $(tail -n+2 {params.indir2}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -eq 0 ]
			then
				if [[ -s {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					bcftools view {params.indir3}/"$contig"_phased_all_variants.vcf.gz \
					-s $HOMGAM -R {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed > {params.outdir}/"$contig"/"$contig"_HOMGAM.vcf
					bgzip -f {params.outdir}/"$contig"/"$contig"_HOMGAM.vcf
					tabix {params.outdir}/"$contig"/"$contig"_HOMGAM.vcf.gz
				fi
			fi
		done < {input.filt_list}
		>> {output.sex_depth_diff_sexshared_vcfs_done}
		"""

# 7.5.2 Split sex-linked regions into heterogametes and homogametes (part 2)
rule sex_depth_diff_hetgam_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		sex_depth_diff_hetgam_vcfs_done=intermediate_clust_settings + "temp_sex_link/sex_depth_diff_hetgam_vcfs_done.txt"
	params:
		indir1=intermediate + "sex_depth_difference",
		indir2=intermediate_clust_settings + "sex_linkage",
		indir3=intermediate + "shapeit4/modified",
		outdir=intermediate_clust_settings + "temp_sex_link"
	conda:
		"../envs/bcftools.yml"
	shell:
		"""
		HETGAM=$(cat {input.sample_table} | awk -F'\\t' '$3=="HETGAM" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
		while read line
		do
			contig=$(echo $line | cut -f1)
			mkdir -p {params.outdir}/"$contig"
			if [ $(tail -n+2 {params.indir2}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -eq 0 ]
			then
				if [[ -s {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					bcftools view {params.indir3}/"$contig"_phased_all_variants.vcf.gz \
					-s $HETGAM -R {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed \
					> {params.outdir}/"$contig"/"$contig"_HETGAM.vcf;
					bgzip -f {params.outdir}/"$contig"/"$contig"_HETGAM.vcf
					tabix {params.outdir}/"$contig"/"$contig"_HETGAM.vcf.gz
				fi
			fi
		done < {input.filt_list}
		>> {output.sex_depth_diff_hetgam_vcfs_done}
		"""

# 7.5.3 Re-code homozygose genotypes to haploid for the heterogametic sex
rule sex_depth_diff_diploid2haploid:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		sex_depth_diff_hetgam_vcfs_done=intermediate_clust_settings + "temp_sex_link/sex_depth_diff_hetgam_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sex_depth_diff_diploid2haploid_done=intermediate_clust_settings + "temp_sex_link/sex_depth_diff_diploid2haploid_done.txt"
	conda:
		"../envs/python_bcftools.yml"
	params:
		diploid2haploid="workflow/scripts/diploid2haploid.py",
		indir1=intermediate + "sex_depth_difference",
		indir2=intermediate_clust_settings + "sex_linkage",
		indir3=intermediate_clust_settings + "temp_sex_link"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir2}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -eq 0 ]
			then
				if [[ -s {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					if [ $(bcftools view -H {params.indir3}/"$contig"/"$contig"_HETGAM.vcf.gz | awk '{{ print; exit }}' | wc -l) -gt 0 ]
					then
						###### Re-code homozygose genotypes to haploid for the heterogametic sex
						python3 {params.diploid2haploid} -i {params.indir3}/"$contig"/"$contig"_HETGAM.vcf.gz -o {params.indir3}/"$contig"/"$contig"_HETGAM_sexshared.vcf
						bgzip -f {params.indir3}/"$contig"/"$contig"_HETGAM_sexshared.vcf
						tabix {params.indir3}/"$contig"/"$contig"_HETGAM_sexshared.vcf.gz
					fi
				fi
			fi
		done < {input.filt_list}
		>> {output.sex_depth_diff_diploid2haploid_done}
		"""

# 7.5.4 Merge sexshared haplotypes from the homogametic and heterogametic sex
rule merge_sex_depth_diff_sexshared_hom_het_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		sex_depth_diff_hetgam_vcfs_done=intermediate_clust_settings + "temp_sex_link/sex_depth_diff_hetgam_vcfs_done.txt",
		sex_depth_diff_diploid2haploid_done=intermediate_clust_settings + "temp_sex_link/sex_depth_diff_diploid2haploid_done.txt",
		sex_depth_diff_sexshared_vcfs_done=intermediate_clust_settings + "temp_sex_link/sex_depth_diff_sexshared_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		merge_sex_depth_diff_sexshared_hom_het_vcfs_done=intermediate_clust_settings + "temp_sex_link/merge_sex_depth_diff_sexshared_hom_het_vcfs_done.txt"
	conda:
		"../envs/bcftools.yml"
	params:
		indir1=intermediate + "sex_depth_difference",
		indir2=intermediate_clust_settings + "sex_linkage",
		indir3=intermediate_clust_settings + "temp_sex_link"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(tail -n+2 {params.indir2}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -eq 0 ]
			then
				if [[ -s {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					if [ $(bcftools view -H {params.indir3}/"$contig"/"$contig"_HETGAM.vcf.gz | awk '{{ print; exit }}' | wc -l) -gt 0 ]
					then
						bcftools merge {params.indir3}/"$contig"/"$contig"_HOMGAM.vcf.gz {params.indir3}/"$contig"/"$contig"_HETGAM_sexshared.vcf.gz > {params.indir3}/"$contig"/"$contig"_sexshared_temp4.vcf
						bgzip -f {params.indir3}/"$contig"/"$contig"_sexshared_temp4.vcf
						tabix {params.indir3}/"$contig"/"$contig"_sexshared_temp4.vcf.gz
					fi
				fi
			fi
		done < {input.filt_list}
		>> {output.merge_sex_depth_diff_sexshared_hom_het_vcfs_done}
		"""

# 7.5.5 Order individuals and filter sites
rule filter_sex_depth_diff_sexshared_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		sex_depth_diff_hetgam_vcfs_done=intermediate_clust_settings + "temp_sex_link/sex_depth_diff_hetgam_vcfs_done.txt",
		merge_sex_depth_diff_sexshared_hom_het_vcfs_done=intermediate_clust_settings + "temp_sex_link/merge_sex_depth_diff_sexshared_hom_het_vcfs_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt",
		sample_table=config["sample_table"]
	output:
		filter_sex_depth_diff_sexshared_vcfs_done=intermediate_clust_settings + "hetgam_vcfs/filter_sex_depth_diff_sexshared_vcfs_done.txt"
	conda:
		"../envs/bcfvcftools.yml"
	params:
		indir1=intermediate + "sex_depth_difference",
		indir2=intermediate_clust_settings + "sex_linkage",
		indir3=intermediate_clust_settings + "temp_sex_link",
		outdir=intermediate_clust_settings + "hetgam_vcfs"
	shell:
		"""
		INDS_ORDER=$(cat {input.sample_table} | awk -F'\\t' '$1!="sample_name" {{ print $1 }}' | tr '\\n' ',' | head -c-1)
		while read line
		do
			contig=$(echo $line | cut -f1)
			mkdir -p {params.outdir}/"$contig"
			if [ $(tail -n+2 {params.indir2}/"$contig"/"$contig"_phase_windows.bed | grep "Sex-linked" | wc -l) -eq 0 ]
			then
				if [[ -s {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					if [ $(bcftools view -H {params.indir3}/"$contig"/"$contig"_HETGAM.vcf.gz | awk '{{ print; exit }}' | wc -l) -gt 0 ]
					then
						bcftools view {params.indir3}/"$contig"/"$contig"_sexshared_temp4.vcf.gz -s $INDS_ORDER | \
						vcffixup - > {params.outdir}/"$contig"/"$contig"_sexshared.vcf
						bgzip -f {params.outdir}/"$contig"/"$contig"_sexshared.vcf
						tabix {params.outdir}/"$contig"/"$contig"_sexshared.vcf.gz
					fi
				fi
			fi
		done < {input.filt_list}
		>> {output.filter_sex_depth_diff_sexshared_vcfs_done}
		"""

# 7.6 Extract autosomal regions.

# 7.6.1 Extract autosomal regions (part 1)
rule sex_depth_diff_autosomal_bed:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sex_depth_diff_autosomal_bed_done=intermediate_clust_settings + "autosomal/sex_depth_diff_autosomal_bed_done.txt"
	params:
		indir1=intermediate + "sex_depth_difference",
		indir2=intermediate_clust_settings + "sex_linkage",
		outdir=intermediate_clust_settings + "autosomal"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			mkdir -p {params.outdir}/"$contig"
			if [ $(tail -n+2 {params.indir2}/"$contig"/"$contig"_phase_windows.bed | grep "Autosomal" | wc -l) -gt 0 ]
			then
				cat {params.indir2}/"$contig"/"$contig"_phase_windows.bed | grep -E "#|Autosomal" | cut -f1,2,3 \
				> {params.outdir}/"$contig"/pseudoautosomal.bed
				if [[ -s {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed && $(grep -vc '^#' {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed) -gt 0 ]]
				then
					bedtools subtract -a {params.outdir}/"$contig"/pseudoautosomal.bed \
					-b {params.indir1}/"$contig"/"$contig"_hetgam_dropout.bed \
					> {params.outdir}/"$contig"/autosomal.bed
				else
					cp {params.outdir}/"$contig"/pseudoautosomal.bed {params.outdir}/"$contig"/autosomal.bed
				fi
				rm {params.outdir}/"$contig"/pseudoautosomal.bed
			fi
		done < {input.filt_list}
		>> {output.sex_depth_diff_autosomal_bed_done}
		"""

# 7.6.2 Extract autosomal regions (part 2)
rule sex_depth_diff_autosomal_vcfs:
	input:
		done1=lambda wildcards: expand(intermediate_clust_settings + "sex_linkage/{large_contigs}/{large_contigs}_done.txt", large_contigs=get_large_contigs(wildcards)),
		done2=intermediate_clust_settings + "sex_linkage/clump_sex_linkage_done.txt",
		sex_depth_diff_autosomal_bed_done=intermediate_clust_settings + "autosomal/sex_depth_diff_autosomal_bed_done.txt",
		mod_done=intermediate + "shapeit4/modified/modifications_done.txt",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sex_depth_diff_autosomal_vcfs_done=intermediate_clust_settings + "autosomal/sex_depth_diff_autosomal_vcfs_done.txt"
	conda:
		"../envs/bcftools.yml"
	params:
		indir1=intermediate_clust_settings + "autosomal",
		indir2=intermediate + "shapeit4/modified"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			if [ $(cat {params.indir1}/"$contig"/autosomal.bed | wc -l) -gt 0 ]
			then
				bcftools view {params.indir2}/"$contig"_phased_all_variants.vcf.gz -R {params.indir1}/"$contig"/autosomal.bed > {params.indir1}/"$contig"/"$contig"_autosomal.vcf;
				bgzip -f {params.indir1}/"$contig"/"$contig"_autosomal.vcf
				tabix {params.indir1}/"$contig"/"$contig"_autosomal.vcf.gz
			fi
		done < {input.filt_list}
		>> {output.sex_depth_diff_autosomal_vcfs_done}
		"""