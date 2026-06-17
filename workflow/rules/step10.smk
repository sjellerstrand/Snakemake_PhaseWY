################## Step 10 ###########################
# 10.1 Print vcf statistics. Note, there are several 
# files that are listed as input that are not called 
# in the shell command. This is to ensure the previous 
# rules run first
rule vcf_stats_pt1:
	input:
		results_vcf=results_clust_settings + "vcfs/{results_vcf}.vcf.gz",
		all_phased=results_clust_settings + "vcfs/phased_all_variants.vcf.gz",
		nonphased=results_clust_settings + "vcfs/nonphased_variants.vcf.gz",
		border_variants=results_clust_settings + "vcfs/border_variants.vcf.gz",
		autosomal=results_clust_settings + "vcfs/autosomal.vcf.gz",
		sexshared=results_clust_settings + "vcfs/sexshared.vcf.gz",
		sexlimited=results_clust_settings + "vcfs/sexlimited.vcf.gz",
		all_ILS1=results_clust_settings + "vcfs/all_inds_ILS1_not_hetgam_dropout.vcf.gz",
		homgam_ILS1=results_clust_settings + "vcfs/sexshared_ILS1_not_hetgam_dropout.vcf.gz",
		sexlimited_ILS1=results_clust_settings + "vcfs/sexlimited_ILS1_not_hetgam_dropout.vcf.gz",
		all_ILS2=results_clust_settings + "vcfs/all_inds_ILS2.vcf.gz",
		homgam_ILS2=results_clust_settings + "vcfs/sexshared_ILS2.vcf.gz",
		sexlimited_ILS2=results_clust_settings + "vcfs/sexlimited_ILS2.vcf.gz",
		singletons=results_clust_settings + "vcfs/heterogametic_sexlinked_singletons.vcf.gz",
		unreliable=results_clust_settings + "vcfs/all_unreliable_sites.vcf.gz",
		all_phased_done=intermediate_clust_settings + "vcfs/merge_phased_marker_vcfs_done.txt",
		nonphased_done=intermediate_clust_settings + "vcfs/nonphased_variants_done.txt",
		border_variants_done=intermediate_clust_settings + "vcfs/border_variants_done.txt",
		autosomal_done=intermediate_clust_settings + "vcfs/merge_autosomal_vcfs_done.txt",
		homgam_sexlimited_done=intermediate_clust_settings + "vcfs/merge_sexshared_and_sexlimited_vcfs_done.txt",
		ILS_done=intermediate_clust_settings + "vcfs/ILS_done.txt",
		singletons_done=intermediate_clust_settings + "vcfs/singletons_done.txt",
		unreliable_done=intermediate_clust_settings + "vcfs/all_unreliable_sites_done.txt"
	output:
		num_var=results_clust_settings + "stats/{results_vcf}/{results_vcf}_number_variants.txt",
		missing=results_clust_settings + "stats/{results_vcf}/{results_vcf}.lmiss",
		freq2=results_clust_settings + "stats/{results_vcf}/{results_vcf}.frq",
		het=results_clust_settings + "stats/{results_vcf}/{results_vcf}.het",
		miss_indiv=results_clust_settings + "stats/{results_vcf}/{results_vcf}.imiss",
		stats_done=intermediate_clust_settings + "stats/{results_vcf}/{results_vcf}.done"
	params:
		logfile=results_clust_settings + "stats/{results_vcf}/{results_vcf}.log"
	conda:
		"../envs/bcfvcftools.yml"
	shell:
		"""
		vcftools --gzvcf {input.results_vcf} --missing-site --stdout >{output.missing} 2>>{params.logfile}
		vcftools --gzvcf {input.results_vcf} --freq2 --max-alleles 2 --stdout >{output.freq2} 2>>{params.logfile}
		vcftools --gzvcf {input.results_vcf} --het --stdout >{output.het} 2>>{params.logfile}
		vcftools --gzvcf {input.results_vcf} --missing-indv --stdout >{output.miss_indiv} 2>>{params.logfile}
 		wc -l {output.missing} | awk '{{print $1-1}}' > {output.num_var} 2>>{params.logfile}
		>> {output.stats_done}
		"""

# 10.2 Run 'filter_stats.r' script.
rule vcf_stats_pt2:
	input:
		stats_done=intermediate_clust_settings + "stats/{results_vcf}/{results_vcf}.done",
		num_var=results_clust_settings + "stats/{results_vcf}/{results_vcf}_number_variants.txt"
	output:
		filter_stats_r_done=intermediate_clust_settings + "stats/{results_vcf}/{results_vcf}_filter_stats_rscript.done"
	params:
		filter_stats_r="workflow/scripts/filter_stats.r",
		vcf_out=results_clust_settings + "stats/{results_vcf}/{results_vcf}",
		logfile=results_clust_settings + "stats/{results_vcf}/{results_vcf}_filter_stats_rscript.log"
	conda:
		"../envs/r-plots.yml"
	shell:
		"""
		num_var=$(cut -f1 {input.num_var})
		if (($num_var > 0))
		then
			Rscript {params.filter_stats_r} --no-save --args VCF_OUT={params.vcf_out} 2>{params.logfile}
			echo 'filter_stats r-script done' > {output.filter_stats_r_done}
		else
			echo 'filter_stats r-script not run. Too few variants' > {output.filter_stats_r_done}
		fi
		"""

# 10.3 Calculate principal components if the number
# of variants is sufficient. Note the stats_done dependancy is 
# included to ensure that the number of variants is known
rule calc_pca:
	input:
		results_vcf=results_clust_settings + "vcfs/{results_vcf}.vcf.gz",
		stats_done=intermediate_clust_settings + "stats/{results_vcf}/{results_vcf}.done",
		num_var=results_clust_settings + "stats/{results_vcf}/{results_vcf}_number_variants.txt"
	output:
		pca_done=intermediate_clust_settings + "pca/{results_vcf}/{results_vcf}_calc_pca.done"
	params:
		outdir=results_clust_settings + "pca/{results_vcf}",
		outfiles=results_clust_settings + "pca/{results_vcf}/{results_vcf}",
		prune_in=results_clust_settings + "pca/{results_vcf}/{results_vcf}.prune.in"
	conda:
		"../envs/plink.yml"
	shell:
		"""
		mkdir -p {params.outdir}
		num_var=$(cut -f1 {input.num_var})
		if (($num_var > 1))
		then
			zcat {input.results_vcf} | awk '/^#/ {{print; next}} {{site = $1":"$2":"$4":"$5; if(site != last) print; last=site}}' | \
			plink --vcf /dev/stdin --double-id --allow-extra-chr --set-missing-var-ids @:# --vcf-half-call missing --indep-pairwise 50 10 0.1 --geno 0.1 -out {params.outfiles}
			if [[ -s "{params.prune_in}" ]]
				then
				prune_count=$(wc -l < {params.prune_in})
				if (( $prune_count > 1 ))
				then
					zcat {input.results_vcf} | awk '/^#/ {{print; next}} {{site = $1":"$2":"$4":"$5; if(site != last) print; last=site}}' | \
					plink --vcf /dev/stdin --double-id --allow-extra-chr --set-missing-var-ids @:# --vcf-half-call missing --extract {params.prune_in} --pca 20 --out {params.outfiles}
					echo 'PCA calculation done' > {output.pca_done}
				else
					echo 'PCA calculation not done. Too few variants' > {output.pca_done}
				fi
			else
				echo 'PCA calculation not done. Too few individuals' > {output.pca_done}
			fi
		else
			echo 'PCA calculation not done. Too few variants' > {output.pca_done}
		fi
		"""

# 10.4 Plot principal component 1 and 2 if there are at least 2 eigenvectors
rule plot_pc1_pc2:
	input:
		calc_done=intermediate_clust_settings + "pca/{results_vcf}/{results_vcf}_calc_pca.done",
		sample_table=config["sample_table"]
	output:
		pca_done=intermediate_clust_settings + "pca/{results_vcf}/{results_vcf}_pca.done"
	params:
		eigenvec=results_clust_settings + "pca/{results_vcf}/{results_vcf}.eigenvec",
		pca_r="workflow/scripts/pca.r",
		vcf_out=results_clust_settings + "pca/{results_vcf}/{results_vcf}",
		logfile=results_clust_settings + "pca/{results_vcf}/{results_vcf}_pca.log"
	conda:
		"../envs/r-plots.yml"
	shell:
		"""
		if (( $(cat {params.eigenvec} | head -n1 | awk '{{print NF}}') > 1 ))
		then
			Rscript {params.pca_r} --no-save --args VCF_OUT={params.vcf_out} INDS={input.sample_table} 2>{params.logfile}
			echo 'plotting complete' > {output.pca_done}
		else
			echo 'plotting not done. Too few variants' > {output.pca_done}
		fi
		"""

# 10.5 Plot haplotype clustering info as manhattan plots
rule plot_haplotype_cluster:
	input:
		phase_windows=results_clust_settings + "beds/phase_windows.bed",
		index=intermediate + "subsetted_genome.fasta.fai"
	output:
		hap_clust_done=intermediate_clust_settings + "genome_summary/haplotype_cluster_summary_done.txt"
	params:
		haplotype_cluster_r="workflow/scripts/haplotype_cluster_plot.r",
		width=config["width"],
		height=config["height"],
		min_len=config["min_len"],
		logfile=results_clust_settings + "genome_summary/haplotype_cluster.log",
		outdir=results_clust_settings + "genome_summary/"
	conda:
		"../envs/r-plots.yml"
	shell:
		"""
		mkdir -p {params.outdir}
		Rscript {params.haplotype_cluster_r} --no-save --args DATA={input.phase_windows} INDEX={input.index} OUT={params.outdir} WIDTH={params.width} HEIGHT={params.height} MIN_LEN={params.min_len} 2>{params.logfile}
		cp {input.index} {params.outdir}
		echo 'plotting complete' > {output.hap_clust_done}
		"""

# 10.6 Genome summary
rule make_genome_summary_beds:
	input:
		bed1=results_clust_settings + "beds/target_region.bed",
		bed2=results_clust_settings + "beds/autosomal.bed",
		bed3=results_clust_settings + "beds/target_phase_hetgam_dropout.bed",
		bed4=results_clust_settings + "beds/callable_phase_sex_linked.bed",
		idx=intermediate + "subsetted_genome.fasta.fai"
	output:
		bed1=intermediate_clust_settings + "genome_summary/Missing_data.bed",
		bed2=intermediate_clust_settings + "genome_summary/Autosomal.bed",
		bed3=intermediate_clust_settings + "genome_summary/Sex_phase_depth_difference.bed",
		bed4=intermediate_clust_settings + "genome_summary/Sex_depth_difference.bed",
		bed5=intermediate_clust_settings + "genome_summary/Sex_phase_difference.bed",
		bed6=results_clust_settings + "genome_summary/Genome_summary.bed"
	params:
		bed3tmp=intermediate_clust_settings + "genome_summary/tmp_Sex_phase_depth_difference.bed"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		bedtools complement -i {input.bed1} -g {input.idx} | awk '{{print $0"\tMissing data"}}' > {output.bed1}
		bedtools intersect -a {input.bed1} -b {input.bed2} | awk '{{print $0"\tAutosomal"}}' > {output.bed2}
		bedtools intersect -a {input.bed3}  -b {input.bed4} | awk '{{print $0"\tSex haplotype clustering & depth difference"}}' > {params.bed3tmp}
		bedtools subtract -a {input.bed3}  -b {params.bed3tmp} | awk '{{print $0"\tSex depth difference"}}' > {output.bed4}
		bedtools subtract -a {input.bed4} -b {params.bed3tmp} | awk '{{print $0"\tSex haplotype clustering"}}'  > {output.bed5}
		mv {params.bed3tmp} {output.bed3}
		cat {output.bed1} {output.bed2} {output.bed3} {output.bed4} {output.bed5} | bedtools sort -g {input.idx} \
		> {output.bed6}
		"""

# 10.7 Plot genome summary
rule plot_genome_summary:
	input:
		genome_summary=results_clust_settings + "genome_summary/Genome_summary.bed",
		index=intermediate + "subsetted_genome.fasta.fai"
	output:
		genome_summary_done=intermediate_clust_settings + "genome_summary/genome_summary.done"
	params:
		genome_summary_r="workflow/scripts/genome_summary_plot.r",
		out_plot=results_clust_settings + "genome_summary/Genome_summary",
		width=config["width"],
		height=config["height"],
		min_len=config["min_len"],
		logfile=results_clust_settings + "genome_summary/Genome_summary.log"
	conda:
		"../envs/r-plots.yml"
	shell:
		"""
		Rscript {params.genome_summary_r} --no-save --args DATA={input.genome_summary} INDEX={input.index} OUT={params.out_plot} WIDTH={params.width} HEIGHT={params.height} MIN_LEN={params.min_len} 2>{params.logfile}
		echo 'plotting complete' > {output.genome_summary_done}
		"""

# 10.8 Create sex depth difference summary
rule sexdiff_genome_summary_bed:
	input:
		idx=intermediate + "subsetted_genome.fasta.fai",
		sex_diff_subsample_done=intermediate + "sex_depth_difference/clump_sex_diff_subsample_" + str(subsample) + "_percent.done",
		bed1=results_clust_settings + "beds/callable_regions_phased.bed",
		bed2=intermediate_clust_settings + "genome_summary/Autosomal.bed",
		bed3=intermediate_clust_settings + "genome_summary/Sex_phase_depth_difference.bed",
		bed4=intermediate_clust_settings + "genome_summary/Sex_depth_difference.bed",
		bed5=intermediate_clust_settings + "genome_summary/Sex_phase_difference.bed"
	output:
		bed6=results_clust_settings + "genome_summary/sexdiff_genome_summary_subset_" +  str(subsample) + "_percent.bed"
	params:
		bed3tmp=intermediate_clust_settings + "genome_summary/tmp_Sex_phase_depth_difference.bed",
		indir=intermediate + "sex_depth_difference/",
		outdir=intermediate_clust_settings + "genome_summary/"
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		awk -F'\\t' -v OFS='\\t' '{{print $1, $2, $3, $4}}' $(find {params.indir}/ -name "*_sex_diff_subsample.txt") | bedtools sort -g {input.idx} | \
		bedtools intersect -a - -b {input.bed1} > {params.outdir}sex_diff_subset.bed
		bedtools intersect -a {params.outdir}sex_diff_subset.bed -b {input.bed2} | awk -F'\\t' -v OFS='\\t' '{{print $4, "Autosomal"}}' > {output.bed6}
		bedtools intersect -a {params.outdir}sex_diff_subset.bed -b {input.bed3} | awk -F'\\t' -v OFS='\\t' '{{print $4, "Sex haplotype clustering & depth difference"}}' >> {output.bed6}
		bedtools intersect -a {params.outdir}sex_diff_subset.bed -b {input.bed4} | awk -F'\\t' -v OFS='\\t' '{{print $4, "Sex depth difference"}}' >> {output.bed6}
		bedtools intersect -a {params.outdir}sex_diff_subset.bed -b {input.bed5} | awk -F'\\t' -v OFS='\\t' '{{print $4, "Sex haplotype clustering"}}' >> {output.bed6}
		"""

# 10.9 Plot distribution of depth differences.
# A fraction of all sites have been subsampled.
rule plot_sex_depth_distirbution:
	input:
		sexdiff=results_clust_settings + "genome_summary/sexdiff_genome_summary_subset_" +  str(subsample) + "_percent.bed"
	output:
		done=intermediate_clust_settings + "genome_summary/sexdiff_genome_summary_subset_" +  str(subsample) + "_percent.done"
	params:
		sex_depth_distirbution_r="workflow/scripts/sex_depth_distirbution_plot.r",
		out_plot=results_clust_settings + "genome_summary/sexdiff_genome_summary_subset_" +  str(subsample) + "_percent",
		depth_diff=sex_depth_threshold,
		logfile=results_clust_settings + "genome_summary/sex_depth_distirbution_plot.log"
	conda:
		"../envs/r-plots.yml"
	shell:
		"""
		Rscript {params.sex_depth_distirbution_r} --no-save --args DATA={input.sexdiff} OUT={params.out_plot} THRESHOLD={params.depth_diff} 2>{params.logfile}
		>> {output.done}
		"""

# 10.10 Perform genome alignment to synteny genome 
if config.get("synteny_genome") and os.path.exists(config["synteny_genome"]):
	rule synteny_alignment:
		input:
			genome=config["genome"],
			synteny=config["synteny_genome"]
		output:
			synteny_fwd=results_synteny + "fwd_" + synteny_name + ".chain",
			synteny_rev=results_synteny + "rev_" + synteny_name + ".chain",
			synteny_alignment_done=results_synteny + synteny_name + "_synteny_alignment.done"
		params:
			synteny_name=synteny_name,
			outdir=results_synteny,
			crossmap_delta_to_chain="workflow/scripts/crossmap_delta_to_chain.pl"
		threads: 20 # default if nothing else is set in slurm/config.yaml
		conda:
			"../envs/synteny_liftover.yml"
		shell:
			"""
			nucmer -l 20 -c 30 -b 1000 -t {threads} {input.genome} {input.synteny} -p {params.outdir}{params.synteny_name}
			delta-filter -1 {params.outdir}{params.synteny_name}.delta > {params.outdir}{params.synteny_name}.filtered.delta
			perl {params.crossmap_delta_to_chain} -fwd {output.synteny_fwd} -rev {output.synteny_rev} {params.outdir}{params.synteny_name}.filtered.delta
			echo 'synteny alignment complete' > {output.synteny_alignment_done}
			"""

	# 10.11 Perform lift over to synteny genome 
	rule liftover:
		input:
			synteny_alignment_done=results_synteny + synteny_name + "_synteny_alignment.done",
			synteny_fwd=results_synteny + "fwd_" + synteny_name + ".chain",
			genome=config["genome"],
			synteny=config["synteny_genome"],
			genome_summary=results_clust_settings + "genome_summary/Genome_summary.bed"
		output:
			liftover_done=intermediate_clust_settings + "genome_summary/" + synteny_name + "_liftover.done",
			liftover_genome_summary=intermediate_clust_settings + "genome_summary/" + synteny_name + "_liftover_Genome_summary_unprocessed.bed",
			liftover_genome_summary2=results_clust_settings + "genome_summary/" + synteny_name + "_liftover_Genome_summary.bed"
		params:
			synteny_name=synteny_name,
			outdir=intermediate_clust_settings + "genome_summary/",
			index=config["synteny_genome"] + ".fai"
		conda:
			"../envs/synteny_liftover.yml"
		shell:
			"""
			cat {input.genome_summary} | tr ' ' ':' > {params.outdir}Genome_summary_modified.bed
			CrossMap bed {input.synteny_fwd} {params.outdir}Genome_summary_modified.bed | grep -v Unmap | cut -f6- | uniq | tr ':' ' ' | bedtools sort -g {params.index} > {output.liftover_genome_summary}
			bedtools complement -i {output.liftover_genome_summary} -g {params.index} | awk -F'\\t' -v OFS='\\t' '{{print $1,$2,$3,"No alignment"}}' > {params.outdir}NoAlign.bed
			awk -F'\\t' '{{if($4=="Missing data") print}}' {output.liftover_genome_summary}  | bedtools merge -c 4 -o distinct > {params.outdir}Miss.bed
			awk -F'\\t' '{{if($4=="Autosomal") print}}' {output.liftover_genome_summary} | bedtools merge -c 4 -o distinct > {params.outdir}Auto.bed
			awk -F'\\t' '{{if($4=="Sex haplotype clustering & depth difference") print}}' {output.liftover_genome_summary} | bedtools merge -c 4 -o distinct > {params.outdir}SexClustDepth.bed
			awk -F'\\t' '{{if($4=="Sex depth difference") print}}' {output.liftover_genome_summary} | bedtools merge -c 4 -o distinct > {params.outdir}SexDepth.bed
			awk -F'\\t' '{{if($4=="Sex haplotype clustering") print}}' {output.liftover_genome_summary} | bedtools merge -c 4 -o distinct > {params.outdir}SexClust.bed
			cat {params.outdir}NoAlign.bed {params.outdir}Miss.bed {params.outdir}Auto.bed {params.outdir}SexClustDepth.bed {params.outdir}SexDepth.bed {params.outdir}SexClust.bed > {output.liftover_genome_summary2}
			echo 'liftover complete' > {output.liftover_done}
			"""

	# 10.12 Plot genome summary
	rule plot_liftover_genome_summary:
		input:
			liftover_done=intermediate_clust_settings + "genome_summary/" + synteny_name + "_liftover.done",
			liftover_genome_summary=results_clust_settings + "genome_summary/" + synteny_name + "_liftover_Genome_summary.bed"
		output:
			genome_summary_done=intermediate_clust_settings + "genome_summary/" + synteny_name + "_liftover_genome_summary.done"
		params:
			genome_summary_r="workflow/scripts/genome_summary_plot.r",
			out_plot=results_clust_settings + "genome_summary/" + synteny_name + "_liftover_genome_summary",
			width=config["width"],
			height=config["height"],
			min_len=config["min_len"],
			logfile=results_clust_settings + "genome_summary/" + synteny_name + "_liftover_genome_summary.log",
			index=config["synteny_genome"] + ".fai",
			outdir=results_clust_settings + "genome_summary/" + config["synteny_species"] + ".fasta.fai"
		conda:
			"../envs/r-plots.yml"
		shell:
			"""
			Rscript {params.genome_summary_r} --no-save --args DATA={input.liftover_genome_summary} INDEX={params.index} OUT={params.out_plot} WIDTH={params.width} HEIGHT={params.height} MIN_LEN={params.min_len} 2>{params.logfile}
			cp {params.index} {params.outdir}
			echo 'plotting complete' > {output.genome_summary_done}
			"""
