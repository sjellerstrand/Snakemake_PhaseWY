########################## Step 1 ##########################

# 1.1 Subset genome.fasta based on user supplied
# list of contigs.
rule subset_genome:
	input:
		genome=config["genome"],
		contig_list=config["contig_list"]
	output:
		fasta=intermediate + "subsetted_genome.fasta"
	conda:
		"../envs/seqtk.yml"
	shell:
		"""
		if [[ "{input.contig_list}" == *.fai ]]
		then
			cp {input.genome} {output.fasta}
		else
			seqtk subseq {input.genome} {input.contig_list} > {output.fasta}
		fi
		"""

# 1.2 index genome and define which contig
# should be run as independent jobs, and 
# which should be clumped
checkpoint index_genome:
	input:
		fasta=intermediate + "subsetted_genome.fasta"
	output:
		index=intermediate + "subsetted_genome.fasta.fai",
		clump_contig=intermediate + "clump_contigs.txt",
		large_contig=intermediate + "large_contigs.txt"
	params:
		clump_max_len=clump_max_len
	conda:
		"../envs/samtools.yml"
	shell:
		"""
		samtools faidx {input.fasta} -o {output.index}
		awk -F'\\t' -v max_len={params.clump_max_len} '{{if($2 < (max_len)*1000000) print $1}}' {output.index} \
		> {output.clump_contig}
		awk -F'\\t' -v max_len={params.clump_max_len} '{{if($2 >= (max_len)*1000000) print $1}}' {output.index} \
		> {output.large_contig}
		"""

# Define which contig should be run as independent jobs
def get_large_contigs(wildcards):
	cp_output = checkpoints.index_genome.get()
	with open(cp_output.output['large_contig']) as f:
		large_contigs = [
			line.strip() for line in f if line.strip()
		]
	return large_contigs

# 1.3 Use samtools depth to get the depth for each site
# from all of the bams. Here, we are using the user specified list of contigs. 
# One output file is generated per contig and includes all of the bams.
# ('-H' flag so there is a header)

# 1.3.1 Large contigs
rule get_depth_large:
	input:
		bams=expand("{bamfile}", bamfile=sample_df.bamfile),
	output:
		table=intermediate + "samtools/{large_contigs}/{large_contigs}_depth_table.txt",
		done=intermediate + "samtools/{large_contigs}/{large_contigs}_depth_table_done.txt"
	params:
		files=lambda wildcards, input: ' '.join(input.bams),
		samples=lambda wildcards, input: '\t'.join(sample_df.sample_name)
	threads: 8
	conda:
		"../envs/samtools.yml"
	shell:
		"""
		samtools depth -aa -J -H -s --threads {threads} -r {wildcards.large_contigs} {params.files} -o {output.table}
		awk 'NR==1 {{print "#CHROM\tPOS\t{params.samples}"}} NR>1 {{print}}' {output.table} > {output.table}.tmp && mv {output.table}.tmp {output.table}
		>> {output.done}
		"""

# 1.3.2 Clumped contigs
rule get_depth_clump:
	input:
		bams=expand("{bamfile}", bamfile=sample_df.bamfile),
		clump_list=intermediate + "clump_contigs.txt",
	output:
		done=intermediate + "samtools/clump_contigs_depth_table_done.txt"
	params:
		files=lambda wildcards, input: ' '.join(input.bams),
		samples=lambda wildcards, input: '\t'.join(sample_df.sample_name),
		outdir=intermediate + "samtools/"
	threads: 8
	conda:
		"../envs/samtools.yml"
	shell:
		"""
		if [ -s {input.clump_list} ]
		then
			while read line
			do
				contig=$(echo $line | cut -f1)
				mkdir -p {params.outdir}"$contig"
				samtools depth -aa -J -H -s --threads {threads} -r "$contig" {params.files} -o {params.outdir}"$contig"/"$contig"_depth_table.txt
				awk 'NR==1 {{print "#CHROM\tPOS\t{params.samples}"}} NR>1 {{print}}' {params.outdir}"$contig"/"$contig"_depth_table.txt \
				> {params.outdir}"$contig"/"$contig"_depth_table.tmp && mv {params.outdir}"$contig"/"$contig"_depth_table.tmp {params.outdir}"$contig"/"$contig"_depth_table.txt
			done < {input.clump_list}
		fi
		>> {output.done}
		"""

# 1.4 Parse the samtools depth files and filter
# sites based on min depth, missingness, min_mean and max_mean.
# This rule outputs a tmp file for each contig that will be
# filtered further in the next rule (filt_mask). The logfile
# contains info about sites removed etc.

# 1.4.1 Large contigs
rule parse_depth_large:
	input:
		depth_table=intermediate + "samtools/{large_contigs}/{large_contigs}_depth_table.txt",
		done=intermediate + "samtools/{large_contigs}/{large_contigs}_depth_table_done.txt"
	output:
		bed=intermediate + "samtools/{large_contigs}/{large_contigs}_depth_table_filt_tmp.bed"
	params:
		min_dp=config["min_dp"],
		missing=config["missing"],
		min_mean=config["min_mean"],
		max_mean=config["max_mean"]
	shell:
		"""
		awk -F'\\t' -v MIN_DP={params.min_dp} -v OFS='\\t' 'NR > 1 {{sum=0; missing=0; for(i=3; i<=NF; i++) {{sum+=$i; if($i < MIN_DP) {{missing++}}}} print $0, sum/(NF-2), 1-(missing/(NF-2))}}' {input.depth_table} | \
		awk -F'\\t' -v MIN_MEAN={params.min_mean} -v MAX_MEAN={params.max_mean} -v MISSING={params.missing} -v OFS='\\t' \
		'{{if($(NF-1) >= MIN_MEAN && $(NF-1) <= MAX_MEAN && $NF >= MISSING) print $1, ($2-1), $2}}' \
		> {output.bed}
		"""

# 1.4.2 Clumped contigs
rule parse_depth_clump:
	input:
		clump_list=intermediate + "clump_contigs.txt",
		done=intermediate + "samtools/clump_contigs_depth_table_done.txt"
	output:
		done=intermediate + "samtools/clump_contigs_depth_table_filt_tmp.done"
	params:
		min_dp=config["min_dp"],
		missing=config["missing"],
		min_mean=config["min_mean"],
		max_mean=config["max_mean"],
		outdir=intermediate + "samtools/"
	shell:
		"""
		if [ -s {input.clump_list} ]
		then
			while read line
			do
				contig=$(echo $line | cut -f1)
				awk -F'\\t' -v MIN_DP={params.min_dp} -v OFS='\\t' 'NR > 1 {{sum=0; missing=0; for(i=3; i<=NF; i++) {{sum+=$i; if($i < MIN_DP) {{missing++}}}} print $0, sum/(NF-2), 1-(missing/(NF-2))}}' {params.outdir}"$contig"/"$contig"_depth_table.txt | \
				awk -F'\\t' -v MIN_MEAN={params.min_mean} -v MAX_MEAN={params.max_mean} -v MISSING={params.missing} -v OFS='\\t' \
				'{{if($(NF-1) >= MIN_MEAN && $(NF-1) <= MAX_MEAN && $NF >= MISSING) print $1, ($2-1), $2}}' \
				> {params.outdir}"$contig"/"$contig"_depth_table_filt_tmp.bed
			done < {input.clump_list}
		fi
		>> {output.done}
		"""

# 1.5 Filter sites that were annotated
# as masked. Then, reformat the bedfile.

# 1.5.1 Large contigs
rule filt_mask_large:
	input:
		bed=intermediate + "samtools/{large_contigs}/{large_contigs}_depth_table_filt_tmp.bed"
	output:
		bed=intermediate + "samtools/{large_contigs}/{large_contigs}_callable_regions.bed"
	params:
		mask=config["mask"]
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		if [[ -n "{params.mask}" && -s "{params.mask}" && $(grep -vc '^#' {params.mask}) -gt 0 ]]
		then
			bedtools merge -i {input.bed} | bedtools subtract -a stdin -b {params.mask} > {output.bed}
		else
			bedtools merge -i {input.bed} > {output.bed}
		fi
		"""

# 1.5.2 Clumped contigs
rule filt_mask_clump:
	input:
		clump_list=intermediate + "clump_contigs.txt",
		done=intermediate + "samtools/clump_contigs_depth_table_filt_tmp.done"
	output:
		done=intermediate + "samtools/clump_contigs_callable_regions.done"
	conda:
		"../envs/bedtools.yml"
	params:
		outdir=intermediate + "samtools/",
		mask=config["mask"]
	shell:
		"""
		if [ -s {input.clump_list} ]
		then
			if [[ -n "{params.mask}" && -s "{params.mask}" && $(grep -vc '^#' {params.mask}) -gt 0 ]]; then
				while read line
					do
					contig=$(echo $line | cut -f1)
					bedtools merge -i {params.outdir}"$contig"/"$contig"_depth_table_filt_tmp.bed | \
					bedtools subtract -a stdin -b {params.mask} > {params.outdir}"$contig"/${{contig}}_callable_regions.bed
				done < {input.clump_list}
			else
				while read line
					do
					contig=$(echo $line | cut -f1)
					bedtools merge -i {params.outdir}"$contig"/"$contig"_depth_table_filt_tmp.bed > {params.outdir}"$contig"/${{contig}}_callable_regions.bed
				done < {input.clump_list}
			fi
		fi
		>> {output.done}
		"""