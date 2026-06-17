########################## Step 2 ##########################

# 2.1 Calculate the normalised depth ratio per site between 
# sexes using the script 'sex_depth_diff_contig.r'

# 2.1.1 Large contigs
rule sex_depth_difference_large:
	input:
		depth=intermediate + "samtools/{large_contigs}/{large_contigs}_depth_table.txt",
		sample_table=config["sample_table"]
	output:
		sex_depth_difference=intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_sex_diff.txt"
	params:
		outdir=intermediate + "sex_depth_difference/{large_contigs}"
	shell:
		"""
		cut -f1,2 {input.depth} | tail -n+2 > {params.outdir}/coords.txt
		while read line
		do
			sample=$(echo $line | cut -d" " -f1)
			ind_dp=$(echo $line | cut -d" " -f2)
			sex=$(echo $line | cut -d" " -f3)
			col=$(head -n1 {input.depth} | tr '\\t' '\\n' | awk -v sample=$sample '{{if($1==sample) print NR}}')
			awk -F'\\t' -v col=$col -v ind_dp=$ind_dp '{{if(NR > 1) print $col/ind_dp}}' {input.depth} > {params.outdir}/"$sample"_"$sex"_depth.txt
		done < <(tail -n+2 {input.sample_table})

		files1=$(ls {params.outdir}/*_HETGAM_depth.txt)
		files2=$(ls {params.outdir}/*_HOMGAM_depth.txt)
		paste {params.outdir}/coords.txt $files1 | \
		awk -F'\\t' -v OFS='\\t' '{{sum=0; for(i=3; i<=NF; i++) sum+=$i; avg = sum/(NF-2); print $1, $2, avg}}' > {params.outdir}/HETGAM_avg_depth.txt
		paste {params.outdir}/coords.txt $files2 | \
		awk -F'\\t' -v OFS='\\t' '{{sum=0; for(i=3; i<=NF; i++) sum+=$i; avg = sum/(NF-2); print $1, $2, avg}}' > {params.outdir}/HOMGAM_avg_depth.txt

		paste {params.outdir}/HETGAM_avg_depth.txt {params.outdir}/HOMGAM_avg_depth.txt | \
		awk -F'\\t' -v OFS='\\t' '{{ if($6 == 0) {{ratio=1}} else {{ratio=$3/$6}}; print $1, $2, ratio }}' \
		> {output.sex_depth_difference}
		rm $files1 $files2 {params.outdir}/HETGAM_avg_depth.txt {params.outdir}/HOMGAM_avg_depth.txt {params.outdir}/coords.txt 
		"""

# 2.1.2 Clumped contigs
rule sex_depth_difference_clump:
	input:
		clump_list=intermediate + "clump_contigs.txt",
		done=intermediate + "samtools/clump_contigs_depth_table_done.txt",
		sample_table=config["sample_table"]
	output:
		done=intermediate + "sex_depth_difference/clump_sex_diff.done"
	params:
		indir=intermediate + "samtools",
		outdir=intermediate + "sex_depth_difference"
	shell:
		"""
		if [ -s {input.clump_list} ]
		then
			while read line
			do
				contig=$(echo $line | cut -f1)
				mkdir -p {params.outdir}/"$contig"
				cut -f1,2 {params.indir}/"$contig"/"$contig"_depth_table.txt | tail -n+2 > {params.outdir}/"$contig"/coords.txt
				while read line
					do
					sample=$(echo $line | cut -d" " -f1)
					ind_dp=$(echo $line | cut -d" " -f2)
					sex=$(echo $line | cut -d" " -f3)
					col=$(head -n1 {params.indir}/"$contig"/"$contig"_depth_table.txt | tr '\\t' '\\n' | awk -v sample=$sample '{{if($1==sample) print NR}}')
					awk -F'\\t' -v col=$col -v ind_dp=$ind_dp '{{if(NR > 1) print $col/ind_dp}}' {params.indir}/"$contig"/"$contig"_depth_table.txt \
					> {params.outdir}/"$contig"/"$sample"_"$sex"_depth.txt
				done < <(tail -n+2 {input.sample_table})

				files1=$(ls {params.outdir}/"$contig"/*_HETGAM_depth.txt)
				files2=$(ls {params.outdir}/"$contig"/*_HOMGAM_depth.txt)		

				paste {params.outdir}/"$contig"/coords.txt $files1 | \
				awk -F'\\t' -v OFS='\\t' '{{sum=0; for(i=3; i<=NF; i++) sum+=$i; avg = sum/(NF-2); print $1, $2, avg}}' > {params.outdir}/"$contig"/HETGAM_avg_depth.txt
				paste {params.outdir}/"$contig"/coords.txt $files2 | \
				awk -F'\\t' -v OFS='\\t' '{{sum=0; for(i=3; i<=NF; i++) sum+=$i; avg = sum/(NF-2); print $1, $2, avg}}' > {params.outdir}/"$contig"/HOMGAM_avg_depth.txt

				paste {params.outdir}/"$contig"/HETGAM_avg_depth.txt {params.outdir}/"$contig"/HOMGAM_avg_depth.txt | \
				awk -F'\\t' -v OFS='\\t' '{{ if($6 == 0) {{ratio=1}} else {{ratio=$3/$6}}; print $1, $2, ratio }}' \
				> {params.outdir}/"$contig"/"$contig"_sex_diff.txt
				rm $files1 $files2 {params.outdir}/"$contig"/HETGAM_avg_depth.txt {params.outdir}/"$contig"/HOMGAM_avg_depth.txt {params.outdir}/"$contig"/coords.txt
			done < {input.clump_list}
		fi
		>> {output.done}
		"""

# 2.2 Extract sex-linked regions due to sex depth differences.
# These are intersected with callable regions from Step 1.

# 2.2.1 Large contigs
rule extract_sexlinked_regions_large:
	input:
		sex_depth_difference=intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_sex_diff.txt",
		bed=intermediate + "samtools/{large_contigs}/{large_contigs}_callable_regions.bed"
	output:
		sexlinked_regions=intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed"
	params:
		indir=intermediate + "sex_depth_difference",
		depth_diff=sex_depth_threshold
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		awk -F'\\t' -v OFS='\\t' -v DEPTH_DIFF={params.depth_diff} '{{ if($3 < DEPTH_DIFF) print $1, $2-1, $2 }}' {input.sex_depth_difference} | \
		bedtools merge | bedtools intersect -a - -b {input.bed} \
		> {output.sexlinked_regions}
		"""

# 2.2.2 Clumped contigs
rule extract_sexlinked_regions_clump:
	input:
		clump_list=intermediate + "clump_contigs.txt",
		done=intermediate + "sex_depth_difference/clump_sex_diff.done",
		bed=intermediate + "samtools/clump_contigs_callable_regions.done"
	output:
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done"
	params:
		indir1=intermediate + "samtools",
		indir2=intermediate + "sex_depth_difference",
		depth_diff=sex_depth_threshold
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		if [ -s {input.clump_list} ]
		then
			while read line
				do
				contig=$(echo $line | cut -f1)
				awk -F'\\t' -v OFS='\\t' -v DEPTH_DIFF={params.depth_diff} '{{ if($3 < DEPTH_DIFF) print $1, $2-1, $2 }}' {params.indir2}/"$contig"/"$contig"_sex_diff.txt | \
				bedtools merge | bedtools intersect -a - -b  {params.indir1}/"$contig"/"$contig"_callable_regions.bed \
				> {params.indir2}/"$contig"/"$contig"_hetgam_dropout.bed
			done < {input.clump_list}
		fi
		>> {output.sexlinked_regions_done}
		"""

# 2.3 Calculate windowed avarages

# 2.3.1 Large
rule sex_depth_windows_large:
	input:
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		bed1=lambda wildcards: expand(intermediate + "samtools/{large_contigs}/{large_contigs}_callable_regions.bed", large_contigs=get_large_contigs(wildcards)),
		index=intermediate + "subsetted_genome.fasta.fai",
	output:
		sex_depth_windows_done=intermediate_clust_settings + "sex_depth_windows/{large_contigs}/{large_contigs}_sex_depth_windows_clump.done"
	params:
		indir1=intermediate + "samtools",
		indir2=intermediate + "sex_depth_difference",
		outdir=intermediate_clust_settings + "sex_depth_windows",
		depth_diff=sex_depth_threshold,
		window=window,
		step=step
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		contig="{wildcards.large_contigs}"
		mkdir -p {params.outdir}/"$contig"
		contiglength=$(grep -w "$contig" {input.index} | cut -f2)
		WINDOW={params.window}
		STEP={params.step}

		# Calculate windows
		if (( contiglength < $WINDOW * ($WINDOW / $STEP) ))
		then
			WINDOW0=$WINDOW
			STEP0=$STEP
			WINDOW=$(( ($contiglength + ($WINDOW/$STEP0)-1) / ($WINDOW/$STEP0) ))
			STEP=$(( (($STEP0 * $contiglength) + ((($WINDOW0/$STEP0)-1)* $WINDOW0) - 1) / (($WINDOW0/$STEP0) * $WINDOW0) ))
		fi
		if (( $contiglength == 1 ))
		then
			startpos="0"
			endpos="1"
		else
			startpos=$(seq 1 $STEP $(($contiglength - $WINDOW)))
			echo $startpos | tr ' ' '\n' | awk -v c="$contig" -v OFS='\\t' -v W=$WINDOW '{{print c, $1-1, $1+W-1}}' \
			> {params.outdir}/"$contig"/"$contig"_windows.bed
			if (( $(tail -n1 {params.outdir}/"$contig"/"$contig"_windows.bed | cut -f3) < $contiglength ))
			then
				printf "%s\t%s\t%s\n" "$contig" "$(( $contiglength - $WINDOW ))" "$contiglength" \
				>> {params.outdir}/"$contig"/"$contig"_windows.bed
			fi
		fi
		if [[ -s "{params.indir1}/"$contig"/"$contig"_callable_regions.bed" && $(grep -vc '^#' {params.indir1}/"$contig"/"$contig"_callable_regions.bed) -gt 0 ]]
		then
			awk -F'\\t' -v OFS='\\t' '{{ print $1, $2-1, $2, $3 }}' {params.indir2}/"$contig"/"$contig"_sex_diff.txt | \
			bedtools intersect -a - -b  {params.indir1}/"$contig"/"$contig"_callable_regions.bed | \
			bedtools map -a {params.outdir}/"$contig"/"$contig"_windows.bed -b - -c 4 -o mean -null "NA" \
			> {params.outdir}/"$contig"/"$contig"_sex_depth_windows.bed
		else
			awk -F'\\t' -v OFS='\\t' '{{ print $1, $2, $3, "NA" }}' {params.outdir}/"$contig"/"$contig"_windows.bed \
			> {params.outdir}/"$contig"/"$contig"_sex_depth_windows.bed
		fi
		>> {output.sex_depth_windows_done}
		"""

# 2.3.2 Clumped
rule sex_depth_windows_clump:
	input:
		clump_list=intermediate + "clump_contigs.txt",
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		bed2=intermediate + "samtools/clump_contigs_callable_regions.done",
		index=intermediate + "subsetted_genome.fasta.fai",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sex_depth_windows_done=intermediate_clust_settings + "sex_depth_windows/sex_depth_windows_clump.done"
	params:
		indir1=intermediate + "samtools",
		indir2=intermediate + "sex_depth_difference",
		outdir=intermediate_clust_settings + "sex_depth_windows",
		depth_diff=sex_depth_threshold,
		window=window,
		step=step
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		if [ -s {input.clump_list} ]
		then
			while read line
			do
				contig=$(echo $line | cut -f1)
				mkdir -p {params.outdir}/"$contig"
				contiglength=$(grep -w "$contig" {input.index} | cut -f2)
				WINDOW={params.window}
				STEP={params.step}

				# Calculate windows
				if (( $contiglength < $WINDOW * ($WINDOW / $STEP) ))
				then
					WINDOW0=$WINDOW
					STEP0=$STEP
					WINDOW=$(( ($contiglength + ($WINDOW/$STEP0)-1) / ($WINDOW/$STEP0) ))
					STEP=$(( (($STEP0 * $contiglength) + ((($WINDOW0/$STEP0)-1)* $WINDOW0) - 1) / (($WINDOW0/$STEP0) * $WINDOW0) ))
				fi
				if (( $contiglength == 1 ))
				then
					startpos="0"
					endpos="1"
				else
					startpos=$(seq 1 $STEP $(($contiglength - $WINDOW)))
					echo $startpos | tr ' ' '\n' | awk -v c="$contig" -v OFS='\\t' -v W=$WINDOW '{{print c, $1-1, $1+W-1}}' \
					> {params.outdir}/"$contig"/"$contig"_windows.bed
					if (( $(tail -n1 {params.outdir}/"$contig"/"$contig"_windows.bed | cut -f3) < $contiglength ))
					then
						printf "%s\t%s\t%s\n" "$contig" "$(( $contiglength - $WINDOW ))" "$contiglength" \
						>> {params.outdir}/"$contig"/"$contig"_windows.bed
					fi
				fi
				if [[ -s "{params.indir1}/"$contig"/"$contig"_callable_regions.bed" && $(grep -vc '^#' {params.indir1}/"$contig"/"$contig"_callable_regions.bed) -gt 0 ]]
				then
					awk -F'\\t' -v OFS='\\t' '{{ print $1, $2-1, $2, $3 }}' {params.indir2}/"$contig"/"$contig"_sex_diff.txt | \
					bedtools intersect -a - -b  {params.indir1}/"$contig"/"$contig"_callable_regions.bed | \
					bedtools map -a {params.outdir}/"$contig"/"$contig"_windows.bed -b - -c 4 -o mean -null "NA" \
					> {params.outdir}/"$contig"/"$contig"_sex_depth_windows.bed
				else
					awk -F'\\t' -v OFS='\\t' '{{ print $1, $2, $3, "NA" }}' {params.outdir}/"$contig"/"$contig"_windows.bed \
					> {params.outdir}/"$contig"/"$contig"_sex_depth_windows.bed
				fi
			done < {input.clump_list}
		fi
		>> {output.sex_depth_windows_done}
		"""

# 2.4 Subsample variants for sex depth distribution plot
rule subsample_sexdiff_sites:
	input:
		sexlinked_regions=lambda wildcards: expand(intermediate + "sex_depth_difference/{large_contigs}/{large_contigs}_hetgam_dropout.bed", large_contigs=get_large_contigs(wildcards)),
		sexlinked_regions_done=intermediate + "sex_depth_difference/clump_hetgam_dropout.done",
		bed1=lambda wildcards: expand(intermediate + "samtools/{large_contigs}/{large_contigs}_callable_regions.bed", large_contigs=get_large_contigs(wildcards)),
		bed2=intermediate + "samtools/clump_contigs_callable_regions.done",
		filt_list=intermediate + "filtered_contigs_list.txt"
	output:
		sex_diff_subsample_done=intermediate + "sex_depth_difference/clump_sex_diff_subsample_" + str(subsample) + "_percent.done"
	params:
		indir1=intermediate + "samtools",
		indir2=intermediate + "sex_depth_difference",
		depth_diff=sex_depth_threshold,
		subsample=subsample
	conda:
		"../envs/bedtools.yml"
	shell:
		"""
		while read line
		do
			contig=$(echo $line | cut -f1)
			awk -F'\\t' -v OFS='\\t' '{{ print $1, $2-1, $2, $3 }}' {params.indir2}/"$contig"/"$contig"_sex_diff.txt | \
			bedtools intersect -a - -b  {params.indir1}/"$contig"/"$contig"_callable_regions.bed | \
			awk 'rand() < ({params.subsample}/100)' > {params.indir2}/"$contig"/"$contig"_sex_diff_subsample.txt
		done < {input.filt_list}
		>> {output.sex_diff_subsample_done}
		"""
