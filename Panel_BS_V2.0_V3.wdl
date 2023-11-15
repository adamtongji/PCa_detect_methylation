workflow panel{

	String Sample
	String Fastq1
	String Fastq2
	String Outdir
	String Root
	String Database

	call makedir{
		input:
			outdir=Outdir,
			sample=Sample
	}

	call trim{
		input:
			outdir=makedir.dir,
			root=Root,
			sample=Sample,
			fastq1=Fastq1,
			fastq2=Fastq2,
		R1_pair="R1.fq.gz",
		R1_unpair="R1.unpair.fq.gz",
		R2_pair="R2.fq.gz",
		R2_unpair="R2.unpair.fq.gz",
		threads=15

	}

	call cut{
		input:
			outdir=makedir.dir,
			root=Root,
			sample=Sample,
			R1_pair=trim.R1_pair_file,
			R2_pair=trim.R2_pair_file,
		R1_clean="clean_R1.fq.gz",
		R2_clean="clean_R2.fq.gz",
	}

	call fastqc{
		input:
			outdir=makedir.dir,
			root=Root,
			sample=Sample,
			R1_clean=cut.R1_clean_file,
			R2_clean=cut.R2_clean_file
	}

	call bismark{
		input:
			outdir=makedir.dir,
			root=Root,
			sample=Sample,
			database=Database,
			R1_clean=cut.R1_clean_file,
			R2_clean=cut.R2_clean_file,
		fname1="clean_R1_bismark_bt2_pe"
	}

	call report{
		input:
			outdir=makedir.dir,
			root=Root,
			sample=Sample,
			database=Database,
			bam=bismark.deduplicated_sort_bam,
		fname1="clean_R1_bismark_bt2_pe"
	} 
}


## task.0 Create directory ##
task makedir{

	String outdir
	String sample

	command{
		echo -e "\n------------------------------------------\n[START] ${sample} at [`date +%F` `date +%T`]" > ${outdir}/"Panel.wdl.log" 
		mkdir -p ${outdir}
		mkdir -p ${outdir}/1.trim
		mkdir -p ${outdir}/2.cut
		mkdir -p ${outdir}/3.fastqc
		mkdir -p ${outdir}/4.bismark
		mkdir -p ${outdir}/5.report
		mkdir -p ${outdir}/symbol
	}	

	output{
		String dir="${outdir}"
	}
}

## task.1 trim ##
task trim{

	String outdir
	String root
	String sample
	String fastq1
	String fastq2
	String R1_pair
	String R1_unpair
	String R2_pair
	String R2_unpair
	Int threads=15

	command<<<
		if [ -f ${outdir}/symbol/1.trim.txt ];then
		echo "1.trim task success"
		else

		echo -e "\n\t[START] 1.trim at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		java -jar ${root}/software/trimmomatic-0.36.jar \
		PE -phred33 ${fastq1} ${fastq2} \
		${outdir}/1.trim/${sample}_${R1_pair} ${outdir}/1.trim/${sample}_${R1_unpair} \
		${outdir}/1.trim/${sample}_${R2_pair} ${outdir}/1.trim/${sample}_${R2_unpair} \
		-threads ${threads} \
		ILLUMINACLIP:"${root}/config/adapters/illumina.fa":2:30:10 \
		LEADING:3 \
		TRAILING:3 \
		SLIDINGWINDOW:4:15 \
		MINLEN:51 &&

		echo -e "\t[END]   1.trim at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&

		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/1.trim.txt
		fi
	>>>

	output{
    	File R1_pair_file="${outdir}/1.trim/${sample}_${R1_pair}"
		File R2_pair_file="${outdir}/1.trim/${sample}_${R2_pair}"
	}
}

## task.2 cut ##
task cut{
	String outdir
	String root
	String sample
	String R1_clean
	String R2_clean
	File R1_pair
	File R2_pair

	command<<<
		if [ -f ${outdir}/symbol/2.cut.txt ];then
		echo "2.cut task success"
		else

		# 2.1 Cut 25 bases of Read 1 fq 3' end
		echo -e "\n\t[START] 2.cut:read1 at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		${root}/software/seqtk/seqtk trimfq -b 10 -e 15 ${R1_pair} | gzip > ${outdir}/2.cut/${sample}_${R1_clean} &&

		echo -e "\t[END]   2.cut:read1 at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&

		# 2.2 Cut 25 bases of Read 2 fq 5' end
		echo -e "\n\t[START] 2.cut:read2 at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		${root}/software/seqtk/seqtk trimfq -b 15 -e 10 ${R2_pair} | gzip > ${outdir}/2.cut/${sample}_${R2_clean} &&

		echo -e "\t[END]   2.cut:read2 at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&

		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/2.cut.txt
		fi
	>>>

	output{
		File R1_clean_file="${outdir}/2.cut/${sample}_${R1_clean}"
		File R2_clean_file="${outdir}/2.cut/${sample}_${R2_clean}"
	}

}

## task.3 fastqc ##
task fastqc{
	String outdir
	String root
	String sample
	File R1_clean
	File R2_clean

	command<<<
		if [ -f ${outdir}/symbol/3.fastqc.txt ];then
		echo "3.fastqc task success"
		else

		echo -e "\n\t[START] 3.fastqc at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&
		
		${root}/software/fastqc -t 8 -o ${outdir}/3.fastqc ${R1_clean} ${R2_clean} &&
		
		echo -e "\t[END]   3.fastqc at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&

		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/3.fastqc.txt
		fi
	>>>


}


## task.4 bismark ##
task bismark{
	String outdir
	String root
	String sample
	String database
	String fname1
	File R1_clean
	File R2_clean

	command<<<
		if [ -f ${outdir}/symbol/4.bismark.txt ];then
		echo "4.bismark task success"
		else

		# 4.1 bismark bowtie2 alignment + samtools sort -n
		# output: 
		# 1. bt2.bam: contains all alignments plus methylation call strings
		# 2. PE_report.txt: contains alignment and methylation summary
		echo -e "\n\t[START] 4.1 bismark bowtie2 alignment at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		bismark ${database}/ \
		-1 ${R1_clean} \
		-2 ${R2_clean} \
		--parallel 4 \
		-o ${outdir}/4.bismark &&

		echo -e "\t[END]   4.1 bismark bowtie2 alignment at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&	


		# 4.2 bismark deduplicate_bismark
		# deduplicate the Bismark alignment BAM file and remove all reads but one which align to the the very same position and in the same orientation. 
		# output:
		# 1. sort-n.bam: sorted the data by chromosomes/contigs/scaffoldsz
		# 2. sort-n.deduplicated.bam
		# 3. sort-n.deduplication_report.txt
		# 4. sort-n.deduplicated.sort.bam
		echo -e "\n\t[START] 4.2 bismark deduplicate_bismark at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		samtools sort \
		-n ${outdir}/4.bismark/${sample}_${fname1}.bam \
		-o ${outdir}/4.bismark/${sample}_${fname1}.sort-n.bam &&

		deduplicate_bismark ${outdir}/4.bismark/${sample}_${fname1}.sort-n.bam --output_dir ${outdir}/4.bismark &&

		samtools sort \
		${outdir}/4.bismark/${sample}_${fname1}.sort-n.deduplicated.bam  \
		-o ${outdir}/4.bismark/${sample}_${fname1}.sort-n.deduplicated.sort.bam &&

		echo -e "\t[END]   4.2 bismark deduplicate_bismark at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&
		

		# 4.3 picard.jar CollectInsertSizeMetrics 
		# provides useful metrics for validating library construction including the insert size distribution and read orientation of paired-end libraries.
		# output:
		# 1. sort.bam
		# 2. insert_size_metrics.txt 
		# 3. insert_size_histogram.pdf
		echo -e "\n\t[START] 4.3 picard.jar CollectInsertSizeMetrics at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		samtools sort \
		${outdir}/4.bismark/${sample}_${fname1}.bam \
		-o ${outdir}/4.bismark/${sample}_${fname1}.sort.bam &&

		# java -jar ${root}/software/picard.jar CollectInsertSizeMetrics \
		# I=${outdir}/4.bismark/${sample}_${fname1}.sort.bam \
		# O=${outdir}/4.bismark/${sample}_${fname1}.sort.insert_size_metrics.txt \
		# H=${outdir}/4.bismark/${sample}_${fname1}.sort.insert_size_histogram.pdf \
		# M=0.5 &&

		# echo -e "\t[END]   4.3 picard.jar CollectInsertSizeMetrics at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&	


		# 4.4 picard.jar MarkDuplicates 
		# locates and tags duplicate reads in a BAM or SAM file
		# output:
		# 1. sort.dupMarked.bam: a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. 
		#                   Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. 
		# 2. sort.dupMarked.txt: a metrics file indicating the numbers of duplicates for both single- and paired-end reads.
		echo -e "\n\t[START] 4.4 picard.jar MarkDuplicates at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		java -jar ${root}/software/picard.jar MarkDuplicates \
		I=${outdir}/4.bismark/${sample}_${fname1}.sort.bam  \
		O=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.bam \
		M=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.txt &&

		echo -e "\t[END]   4.4 picard.jar MarkDuplicates at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&	

		# 4.6 bismark bismark_methylation_extractor
		# !!! must sort the bam using 'samtools sort -n' (by read name)
		# extract context-dependent (CpG/CHG/CHH) methylation
		# output:
		# 1. CpG_context_*_sort-n.deduplicated.sort
		# 2. CHG_context_*_sort-n.deduplicated.sort
		# 3. CHH_context_*_sort-n.deduplicated.sort
		# 4. T1_sort-n.deduplicated.sort_splitting_report.txt
		# 5. T1_sort-n.deduplicated.sort.M-bias.txt
		# 6. T1_sort-n.deduplicated.sort.bismark.cov.gz
		# 7. T1_sort-n.deduplicated.sort.bedGraph.gz
		echo -e "\n\t[START] 4.8 bismark bismark_methylation_extractor at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		bismark_methylation_extractor --gzip --comprehensive \
		--bedGraph ${outdir}/4.bismark/${sample}_${fname1}.sort-n.bam \
		--output ${outdir}/4.bismark &&

		echo -e "\t[END]   4.8 bismark bismark_methylation_extractor at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&	


		# 4.7 zcat|grep Lambda|awk
		# Phage lambda DNA as an internal control for quality control, it can accurately evaluate the efficiency of BS processing for each sample
		# output:
		# 1. conversion_rate.txt
		echo -e "\n\t[START] 4.9 zcat|grep Lambda|awk at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		zcat ${outdir}/4.bismark/${sample}_${fname1}.sort-n.bismark.cov.gz | grep Lambda | awk '{a+=$5;b+=$6}END{print 1-a/(a+b),a,b}' > ${outdir}/4.bismark/conversion_rate.txt &&

		zcat ${outdir}/4.bismark/${sample}_${fname1}.sort-n.bismark.cov.gz | grep Lambda | awk '$3 < 7000' | awk '{a+=$5;b+=$6}END{print 1-a/(a+b),a,b}' > ${outdir}/4.bismark/conversion_rate_under_7000.txt &&

		echo -e "\t[END]   4.9 zcat|grep Lambda|awk at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&	


		# 4.8 bismark bam2nuc
		# The script bam2nuc reads BAM files and calculates the mono- and di-nucleotide coverage of the reads and compares it to the average genomic sequence composition. 
		# output:
		# 1. dupMarked.sort.deduplicated.nucleotide_stats.txt
		echo -e "\n\t[START] 4.10 bismark bam2nuc at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		bam2nuc --genome_folder ${database} ${outdir}/4.bismark/${sample}_${fname1}.sort-n.deduplicated.sort.bam --dir ${outdir}/4.bismark &&

		echo -e "\t[END]   4.10 bismark bam2nuc at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&	


		# 4.9 bismark bismark2report
		# find Bismark alignment, deduplication and methylation extraction (splitting) reports as well as M-bias files to generate a graphical HTML report
		# output:
		# 1. T1_clean_R1_bismark_bt2_PE_report.html
		echo -e "\n\t[START] 4.11 bismark bismark2report at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		bismark2report \
		--alignment_report ${outdir}/4.bismark/${sample}_${fname1}_PE_report.txt \
		--dedup_report ${outdir}/4.bismark/${sample}_${fname1}.sort-n.deduplication_report.txt \
		--dir ${outdir}/4.bismark &&

		echo -e "\t[END]   4.11 bismark bismark2report at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&
		echo -e "\n------------------------------------------\n[END]   ${sample} at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" 		
		
		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/4.bismark.txt

		fi
	>>>

	output{
		File deduplicated_sort_bam = "${outdir}/4.bismark/${sample}_${fname1}.sort-n.deduplicated.sort.bam"
	}

}

## task.5 report ##
task report{
	String outdir
	String root
	String sample
	String database
	String fname1
	File bam

	command<<<
		if [ -f ${outdir}/symbol/5.report.txt ];then
		echo "5.report task success"
		else

		# 4.5 picard.jar CollectHsMetrics
		# collects metrics that are specific for sequence datasets generated through hybrid-selection. 
		# Hybrid-selection (HS) is the most commonly used technique to capture exon-specific sequences for targeted sequencing experiments such as exome sequencing
		# NOTE: This tool requires an aligned SAM or BAM file as well as bait and target interval files in Picard interval_list format.
		# NOTE: Metrics labeled as percentages are actually expressed as fractions!
		# output:
		# 1. hs_metrics.txt
		# 2. PER_TARGET_COVERAGE:  getting G/C content and mean sequence depth information for every target interval
		# echo -e "\n\t[START] 4.5 picard.jar CollectHsMetrics at [`date +%F` `date +%T`]" >> ${outdir}/"Panel.wdl.log" &&

		# java -jar ${root}/software/picard.jar CollectHsMetrics \
		# I=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.bam  \
		# O=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.covered_PC1.hs_metrics.txt \
		# PER_TARGET_COVERAGE=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.covered_PC1.PER_TARGET_COVERAGE \
		# MINIMUM_MAPPING_QUALITY=0  \
		# R=${database}/genome.fasta \
		# BAIT_INTERVALS=${root}/config/panel_bed/covered_PC1.interval_list \
		# TARGET_INTERVALS=${root}/config/panel_bed/covered_PC1.interval_list &&

		# java -jar ${root}/software/picard.jar CollectHsMetrics \
		# I=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.bam  \
		# O=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.target-region.hs_metrics.txt \
		# PER_TARGET_COVERAGE=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.target-region.PER_TARGET_COVERAGE \
		# MINIMUM_MAPPING_QUALITY=0  \
		# R=${database}/genome.fasta \
		# BAIT_INTERVALS=${root}/config/panel_bed/target-region.interval_list \
		# TARGET_INTERVALS=${root}/config/panel_bed/target-region.interval_list &&

		# java -jar ${root}/software/picard.jar CollectHsMetrics \
		# I=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.bam  \
		# O=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.ldt_cand.hs_metrics.txt \
		# PER_TARGET_COVERAGE=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.ldt_cand.PER_TARGET_COVERAGE \
		# MINIMUM_MAPPING_QUALITY=0  \
		# R=${database}/genome.fasta \
		# BAIT_INTERVALS=${root}/config/panel_bed/ldt_cand.interval_list \
		# TARGET_INTERVALS=${root}/config/panel_bed/ldt_cand.interval_list &&

		# java -jar ${root}/software/picard.jar CollectHsMetrics \
		# I=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.bam  \
		# O=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.late_target-region.hs_metrics.txt \
		# PER_TARGET_COVERAGE=${outdir}/4.bismark/${sample}_${fname1}.sort.dupMarked.late_target-region.PER_TARGET_COVERAGE \
		# MINIMUM_MAPPING_QUALITY=0  \
		# R=${database}/genome.fasta \
		# BAIT_INTERVALS=${root}/config/panel_bed/late_target-region.interval_list \
		# TARGET_INTERVALS=${root}/config/panel_bed/late_target-region.interval_list &&

		# echo -e "\t[END]   4.5 picard.jar CollectHsMetrics at [`date +%F` `date +%T`]\n" >> ${outdir}/"Panel.wdl.log" &&	

		# # 5.1.1 calculate the score (old versiom)
		# source /home/wangsg/software/miniconda3/etc/profile.d/conda.sh 
		# conda activate frag_3090
		# /home/wangsg/software/miniconda3/envs/frag_3090/bin/python \
		# ${root}/scripts/late_ldt_pipe.py \
		# ${bam} ${sample} \
		# ${outdir}/5.report &&

		# # 5.1.2 score bar plot
		# /usr/local/bin/Rscript ${root}/scripts/percentage_score_bar_plot.R \
		# ${outdir}/5.report \
		# ${outdir}/5.report &&

		# # 5.2.1 calculate the ratio (new version)
		# /home/wangsg/software/miniconda3/envs/frag_3090/bin/python \
		# ${root}/scripts/ratio_preprocess.py \
		# ${bam} ${sample} \
		# ${outdir}/5.report &&

		# # 5.2.2 early_model_v1.R (only score for now)
		# /usr/local/bin/Rscript ${root}/scripts/early_model_v1.R \
		# ${outdir}/5.report/${sample}_bam2ratio.txt \
		# ${outdir}/5.report &&

		# # 5.2.3 late_model_v2.R (score and prob)
		# /usr/local/bin/Rscript ${root}/scripts/late_model_v2.R \
		# ${outdir}/5.report/${sample}_bam2ratio.txt \
		# ${outdir}/5.report &&

		echo "[`date +%F` `date +%T`] Nothing is True. Everything is permitted." > ${outdir}/symbol/5.report.txt
		fi
	>>>

	output{

	}

}