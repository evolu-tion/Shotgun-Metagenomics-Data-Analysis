#!/bin/bash
#SBATCH --job-name=coassembly_analysis
#SBATCH --partition=largemem
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=500G
#SBATCH --ntasks=1
#SBATCH --mail-user=nattawet.sriwichai@oist.jp
#SBATCH --mail-type=END
#SBATCH --input=none
#SBATCH --output=job_coassembly%j.out
#SBATCH --error=job_coassembly%j.err

# Shotgun metagenomic analysis pipeline including 6 step of analysis processes
# 0) Checking qualitiy of seqeuncing data
# 1) Filtering good quality reads of each metatrascriptomic
# 2) Assigning taxonomic of each metagenomic data
# 3) Coassembling read to contig and scaffold
# 4) Calling gene and extract protein seqeunce from assembled read
# 5) Annoting protein sequence
# 6) Comparing fuctional of each metagenomic

# Load configuration file
source configuration
dir_out_qc=${dir_out}/${dir_out_qc}
dir_out_fread=${dir_out}/${dir_out_fread}
dir_out_taxa=${dir_out}/${dir_out_taxa}
dir_out_co_assem=${dir_out}/${dir_out_co_assem}
dir_out_bining=${dir_out}/${dir_out_bining}
dir_out_cgene=${dir_out}/${dir_out_cgene}
dir_out_anno=${dir_out}/${dir_out_anno}
dir_out_reassembly=${dir_out}/${dir_out_reassembly}


list_fa_file=()
list_fa_short_name=()
echo "Please type short name"
for R1 in `cd ${dir_seq} && ls *_R1_*`; do
	fastq_paired_end=(${R1} ${R1//_R1_/_R2_})
	case "$R1" in
		*.gz ) 
			sample_name=${R1%%.*}
			;;
		* )
			sample_name=${R1%.*}
			;;
	esac
	sample_name=${sample_name//_R1_001/}
	sample_name=${sample_name//_/}
	sample_name=${sample_name//-/}
	read -e -p "$(tput setaf 1) + ${R1}: $(tput sgr 0)" -i $sample_name sample_name
	list_fa_short_name+=(${sample_name})
	list_fa_file+=(${fastq_paired_end[0]},${fastq_paired_end[1]})
done
echo "--------------------------------------------"

echo " ___________________________________________________________"
echo "|         FastQ file name          |      Short name        |"
echo "|__________________________________|________________________|"
num_list=${#list_fa_file[@]}
for (( i=0; i<${#list_fa_file[@]}; i++ )); do
	echo "  ${list_fa_file[$i]%,*}     ${list_fa_short_name[$i]}  "
done
echo "|__________________________________|________________________|"

echo "Shotgun metagenomics analyis pipeline including: "
if [ "${step_0_check_quality}" == "Y" ]; then
	echo "+ Run step 0: Check sequence quality"
fi
if [ "${step_1_filtering}" == "Y" ]; then
	echo "+ Run step 1: Filtering read"
fi
if [ "${step_2_taxonomy_assigning}" == "Y" ]; then
	echo "+ Run step 2: Taxonomic assigenment"
fi
if [ "${step_3_1_coassembling}" == "Y" ]; then
	echo "+ Run step 3.1: coassembly metagenomic sample"
fi
if [ "${step_3_2_reassembling}" == "Y" ]; then
	echo "+ Run step 3.2: reassembly read"
fi
if [ "${step_4_gene_calling}" == "Y" ]; then
	echo "+ Run step 4: Gene calling"
fi
if [ "${step_5_annotating}" == "Y" ]; then
	echo "+ Run step 5: gene annotation"
fi
if [ "${display_result}" == "Y" ]; then
	echo "+ Display result"
fi

# # STEP 0: Checking quality of seqeuncing data
if [ "${step_0_check_quality}" == "Y" ]; then
	echo "$(tput setab 7)$(tput setaf 1)========== STEP 0: Checking quality of seqeuncing data ==========$(tput sgr 0)"
	mkdir -p ${dir_out_qc}
	for entry in ${list_fa_file[@]}; do
		fastq_paired_end=(${entry//,/ })
		echo "$(tput setaf 1)Checking seqeuncing read quality: $(tput setaf 2)${fastq_paired_end[0]} and ${fastq_paired_end[1]}$(tput sgr 0)"
		${dir_fastqc} ${dir_seq}/${fastq_paired_end[0]} --threads ${cpu_max} --outdir ${dir_out_qc} &
		${dir_fastqc} ${dir_seq}/${fastq_paired_end[1]} --threads ${cpu_max} --outdir ${dir_out_qc}
	done
	rm ${dir_out_qc}/*.zip
fi

# # STEP 1: Filtering good quality reads of each metatrascriptomic
# # trimmed left and right read that low quality (>= 20 bp) and maximum N is one base on each read.
# # output in directory seq_filtered and sufix file is sample name
if [ "${step_1_filtering}" == "Y" ]; then
	echo "$(tput setab 7)$(tput setaf 1)========== STEP 1: Filtering good quality reads of each metatrascriptomic ==========$(tput sgr 0)"
	mkdir -p ${dir_out_fread}
	mkdir -p ${dir_out_fread}/log
	mkdir -p ${dir_out_fread}/qc
	num=0
	for entry in ${list_fa_file[@]}; do
		fastq_paired_end=(${entry//,/ })
		echo "$(tput setaf 1)Filtering read quality and trimmed: $(tput setaf 2)${list_fa_short_name[num]}$(tput sgr 0)"

		# Trimmomatic filtering read
		java -jar ${dir_trimmomatic} PE \
			-threads 24 \
			-trimlog ${dir_out_fread}/log/${list_fa_short_name[num]}.log \
			${dir_seq}/${fastq_paired_end[0]} \
			${dir_seq}/${fastq_paired_end[1]} \
			${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R1.unpaired.fq.gz \
			${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R2.unpaired.fq.gz \
			ILLUMINACLIP:/home/n/nattawet-sriwichai/bin/software/filter/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
			LEADING:${p_trimming_LEADING} \
			TRAILING:${p_trimming_TRAILING} \
			SLIDINGWINDOW:${p_trimming_SLIDINGWINDOW} \
			MINLEN:${p_trimming_MINLEN} > ${dir_out_fread}/log/${list_fa_short_name[num]}.summary.txt 2>&1

		# Merge filtered read
		cat ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R1.unpaired.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq.gz
		cat ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R2.unpaired.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq.gz

		# Quality checking
		cat ${dir_out_fread}/${list_fa_short_name[num]}_R1.unpaired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R2.unpaired.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}.unpaired.fq.gz

		echo "$(tput setaf 1)Read quality checking: $(tput setaf 2)${list_fa_short_name[num]}_R1.paired.fq.gz and ${list_fa_short_name[num]}_R2.paired.fq.gz$(tput sgr 0)"
		${dir_fastqc} ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz -t ${cpu_max} --outdir ${dir_out_fread}/qc &
		${dir_fastqc} ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz -t ${cpu_max} --outdir ${dir_out_fread}/qc

		num_read_R1=$(cat ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz | echo $((`wc -l`/4)))
		num_read_R2=$(cat ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz | echo $((`wc -l`/4)))
		num_read_unpaired=$(cat ${dir_out_fread}/${list_fa_short_name[num]}.unpaired.fq.gz | echo $((`wc -l`/4)))
		printf "${list_fa_short_name[num]}\t${num_read_R1}\t${num_read_R2}\t${num_read_unpaired}\n" >> num_read.txt
		
		num=$((num+1))
	done
	rm ${dir_out_fread}/qc/*.zip
fi


# # STEP 2: Assigning taxonomy of each metagenomic data
# # Output of this step is taxonomy tab-separation files of each metagenomic data
if [ "${step_2_taxonomy_assigning}" == "Y" ]; then
	echo "$(tput setab 7)$(tput setaf 1)========== STEP 2: Assigning taxonomy of each metagenomic data ==========$(tput sgr 0)"
	mkdir -p ${dir_out_taxa}
	dir_out_metaxa2=${dir_out_taxa}/SSU
	mkdir -p ${dir_out_metaxa2}
	mkdir -p ${dir_out_metaxa2}/visualize_taxonomy

	dir_out_metaphlan2=${dir_out_taxa}/specific_marker_gene
	mkdir -p ${dir_out_metaphlan2}
	mkdir -p ${dir_out_metaphlan2}/visualize_taxonomy

	num_lib=${#list_fa_file[@]}
	for (( num=0; num<${num_lib}; num++ )); do
		echo "$(tput setaf 1)Taxonomy assignment based on SSU by Metaxa2: $(tput setaf 2)${list_fa_short_name[num]}$(tput sgr 0)"
		
		gunzip -c ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq
		gunzip -c ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq
		
		${dir_metaxa} \
			-1 ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq \
			-2 ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq \
			-o ${dir_out_metaxa2}/${list_fa_short_name[num]} \
			-t ${p_taxonomy_organism} \
			-T ${p_taxonomy_identity} \
			-m ${p_taxonomy_match} \
			--cpu=${cpu_max}

		rm ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq
		
		mkdir -p ${dir_out_metaxa2}/${list_fa_short_name[num]}
		mv ${dir_out_metaxa2}/${list_fa_short_name[num]}.* ${dir_out_metaxa2}/${list_fa_short_name[num]}/
		mkdir -p ${dir_out_metaxa2}/${list_fa_short_name[num]}/fasta_marker_16s
		mv ${dir_out_metaxa2}/${list_fa_short_name[num]}/*.fasta ${dir_out_metaxa2}/${list_fa_short_name[num]}/fasta_marker_16s/

		# count the number of species, genera by Metaxa2 Taxonomic Traversal Tool (metaxa2_ttt)
		${dir_metaxa2_ttt} -i ${dir_out_metaxa2}/${list_fa_short_name[num]}/${list_fa_short_name[num]}.taxonomy.txt \
			-o ${dir_out_metaxa2}/${list_fa_short_name[num]}/${list_fa_short_name[num]}
		mkdir -p ${dir_out_metaxa2}/${list_fa_short_name[num]}/taxonomy
		mv ${dir_out_metaxa2}/${list_fa_short_name[num]}/*.level* ${dir_out_metaxa2}/${list_fa_short_name[num]}/taxonomy/

		# taxonomy visualization
		cp ${dir_out_metaxa2}/${list_fa_short_name[num]}/taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.txt ${dir_out_metaxa2}/visualize_taxonomy
		awk -F \\t '{print $2 FS $1}' ${dir_out_metaxa2}/visualize_taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.txt | \
			sed 's/;/\t/g' > ${dir_out_metaxa2}/visualize_taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.converted.txt
		ktImportText ${dir_out_metaxa2}/visualize_taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.converted.txt \
			-o ${dir_out_metaxa2}/visualize_taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.html
		
		# # STEP2.1 high confident of Taxonomy identification based on specific marker gene by metaphlAn2: 
		mkdir -p ${dir_out_metaphlan2}/temp
		echo "$(tput setaf 1)Taxonomy assignment based on specific marker gene by MetaPhlAn2: $(tput setaf 2)${list_fa_short_name[num]}$(tput sgr 0)"
		${dir_MetaPhlAn} ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq.gz,${dir_out_fread}/${list_fa_short_name[num]}_R2.fq.gz \
			--bowtie2out ${dir_out_metaphlan2}/temp/${list_fa_short_name[num]}.bowtie2.bz2 \
			--nproc ${cpu_max} \
			--input_type fastq \
			--ignore_viruses \
			--ignore_eukaryotes \
			--min_cu_len 2000 \
			-t ${p_MPA_analysis_type} \
			> ${dir_out_metaphlan2}/${list_fa_short_name[num]}.profile.txt
	done
	rm ${dir_out_metaxa2}/visualize_taxonomy/*.level_${p_taxonomy_level}.txt 
	python src/taxonomy_comparision.py ${dir_out_metaxa2}/visualize_taxonomy/*.converted.txt > ${dir_out_metaxa2}/taxonomy_comparision.txt
fi


# # STEP 3.1 Co-assembly between data
if [ "${step_3_1_coassembling}" == "Y" ]; then

	echo "$(tput setab 7)$(tput setaf 1)========== STEP 3: Cossembling all sample read to contig and/or scaffold ==========$(tput sgr 0)"
	mkdir -p ${dir_out_co_assem}
	num_data=${#list_fa_file[@]}

	text_R1=""
	text_R2=""
	text_unpaired=""
	for (( i=0; i<${num_data}; i++ )); do
		text_R1="${text_R1}${dir_out_fread}/${list_fa_short_name[i]}_R1.paired.fq.gz "
		text_R2="${text_R2}${dir_out_fread}/${list_fa_short_name[i]}_R2.paired.fq.gz "
		text_unpaired="${text_unpaired}${dir_out_fread}/${list_fa_short_name[i]}.unpaired.fq.gz "
	done

	echo ${text_R1}
	cat ${text_R1} > ${dir_out_fread}/merge_R1.paired.fq.gz

	echo ${text_R2} 
	cat ${text_R2} > ${dir_out_fread}/merge_R2.paired.fq.gz

	echo ${text_unpaired}
	cat ${text_unpaired} > ${dir_out_fread}/merge.unpaired.fq.gz

	echo "$(tput setaf 1)Coassembling by metaSPAdes $(tput sgr 0)"

	${dir_SPAdes} \
		--meta \
		--threads ${cpu_max} \
		--memory ${memory_max} \
		-k 21,33,55,77,99,127 \
		--pe1-1 ${dir_out_fread}/merge_R1.paired.fq.gz \
		--pe1-2 ${dir_out_fread}/merge_R2.paired.fq.gz \
		--pe1-s ${dir_out_fread}/merge.unpaired.fq.gz \
		-o ${dir_out_co_assem}

	python src/merge_contig_scaffold.py \
			--contig ${dir_out_co_assem}/contigs.fasta \
			--output ${dir_out_co_assem}/GOEL.fa \
			--sample "GOEL"

	# # STEP3.1: Check quality of metagenomic assembly by using metaQuast
	echo "Check quality of metagenomic assembly"
	${dir_metaQuast} ${dir_out_co_assem}/contigs.fasta -o ${dir_out_co_assem}/quality --threads=${cpu_max} --max-ref-number 1000

fi



# # convert fastq to paired-end fasta
if [ "${step_3_2_reassembling}" == "Y" ]; then
	echo "========== Reassembly =========="
	mkdir -p ${dir_out_reassembly}
	mkdir -p ${dir_out_reassembly}/temp
	file_name_binning=reassembled

	num_lib=${#list_fa_file[@]}
	for ((i=0;i<$num_lib;i++)) do
		fastq_paired_end=(${entry//,/ })

		echo "convert fastq to fa of : ${list_fa_short_name[$i]} to ${dir_out_reassembly}/temp/${list_fa_short_name[$i]}.fa"
		# python src/fq2fa.py \
		# 	--fq1 ${dir_out_fread}/${list_fa_short_name[$i]}_R1.paired.fq.gz \
		# 	--fq2 ${dir_out_fread}/${list_fa_short_name[$i]}_R2.paired.fq.gz \
		# 	--output ${dir_out_reassembly}/temp/${list_fa_short_name[$i]}.fa
		# echo -e "${dir_out_reassembly}/temp/${list_fa_short_name[$i]}.fa\n" >> ${dir_out_reassembly}/temp/reads_list.txt
	done

	echo "start to bining and reassembly read"
	# ${dir_MaxBin} \
	# 	-contig ${dir_out_co_assem}/contigs.fasta \
	# 	-reads_list ${dir_out_reassembly}/temp/reads_list.txt \
	# 	-out ${dir_out_reassembly}/${file_name_binning} \
	# 	-thread 24 \
	# 	-reassembly

	# ============== SKIP have problem ==============
	# # # assembly by metaSPAdes
	# # separate merage paired-end to single file
	# mkdir -p ${dir_out_reassembly}/temp_assembly
	# mkdir -p ${dir_out_reassembly}/temp_assembly/assembled
	# for bined_reads in `cd ${dir_out_reassembly}/${file_name_binning}.reassem && ls *.reads.*`; do
	# 	python src/fa2split_paired.py \
	# 		--input_single_fasta ${dir_out_reassembly}/${file_name_binning}.reassem/${bined_reads} \
	# 		--output ${dir_out_reassembly}/temp_assembly/sequence/${bined_reads}

	# 	echo "reassembly of bin: ${bined_reads}"
	# 	mkdir -p ${dir_out_reassembly}/temp_assembly/assembled/${bined_reads}
	# 	${dir_SPAdes} \
	# 	--meta \
	# 	--threads ${cpu_max} \
	# 	--memory ${memory_max} \
	# 	-k 21,33,55,77,99,127 \
	# 	--pe1-1 ${dir_out_reassembly}/temp_assembly/sequence/${bined_reads}_R1.paired.fa \
	# 	--pe1-2 ${dir_out_reassembly}/temp_assembly/sequence/${bined_reads}_R2.paired.fa \
	# 	-o ${dir_out_reassembly}/temp_assembly/assembled/${bined_reads}
	# done
	# # assembly (skip)

	# Change heading of each contig and/or merge contig with scaffold
	echo "$(tput setab 7)$(tput setaf 1)Convert bined contig name$(tput sgr 0)"
	mkdir -p ${dir_out_reassembly}/temp_rename_contig_in_bin
	text=""
	for bined_reads in `cd ${dir_out_reassembly} && ls ${file_name_binning}.*.fasta`; do
		echo "$(tput setaf 1)Convert: $(tput setaf 2)${bined_reads}$(tput sgr 0)"
		bined_num=${bined_reads//.fasta/}
		bined_num=${bined_num//${file_name_binning}./}
		# python src/merge_contig_scaffold.py \
		# 	--contig ${dir_out_reassembly}/${bined_reads} \
		# 	--output ${dir_out_reassembly}/temp_rename_contig_in_bin/${bined_num}.fa \
		# 	--sample GO_${bined_num}
		text="${text} ${dir_out_reassembly}/temp_rename_contig_in_bin/${bined_num}.fa"
	done
	echo "$(tput setaf 1)Convert: $(tput setaf 2)${file_name_binning}.noclass$(tput sgr 0)"
	# python src/merge_contig_scaffold.py \
	# 		--contig ${dir_out_reassembly}/${file_name_binning}.noclass \
	# 		--output ${dir_out_reassembly}/temp_rename_contig_in_bin/noclass.fa \
	# 		--sample GO_noclass
	text="${text} ${dir_out_reassembly}/temp_rename_contig_in_bin/noclass.fa"

	cat ${text} > ${dir_out_reassembly}/contig.fa
	echo "$(tput setaf 2)Finish binning method in: ${dir_out_reassembly}/contig.fa$(tput sgr 0)"

	# # STEP4: Calling gene and extract protein seqeunce from assembled read
	mkdir 
	echo "$(tput setab 7)$(tput setaf 1)Calling gene on reassembled contig$(tput sgr 0)"
	mkdir -p ${dir_out_reassembly}/temp_call_gene
	# ${dir_FragGeneScan} -s ${dir_out_reassembly}/contig.fa \
	# 	-o ${dir_out_reassembly}/temp_call_gene/gene \
	# 	-w 0 \
	# 	-t illumina_5 \
	# 	-p ${cpu_max}

	# Change fasta head name
	echo "$(tput setaf 1)Change protein fasta name$(tput sgr 0)"
	# python src/change_gene_name.py -i ${dir_out_reassembly}/temp_call_gene/gene.faa \
	# 	-a ${dir_out_reassembly}/temp_call_gene/gene.out \
	# 	-o ${dir_out_reassembly}/prot.faa
	# python src/change_gene_name.py -i ${dir_out_reassembly}/temp_call_gene/gene.fnn \
	# 	-a ${dir_out_reassembly}/temp_call_gene/gene.out \
	# 	-o ${dir_out_reassembly}/nucl.faa


	# # Annotating protein domain by InterProScan
	# echo "$(tput setab 7)$(tput setaf 2)Protein domain predicting by InterProScan5 $(tput sgr 0)"
	dir_out_anno=annotation
	mkdir -p ${dir_out_reassembly}/${dir_out_anno}/temp_fasta/
	mkdir -p ${dir_out_reassembly}/${dir_out_anno}/protein_domain_prediction/
	mkdir -p ${dir_out_reassembly}/${dir_out_anno}/log_interproscan

	num_lib=${#list_fa_file[@]}
	# Split fasta file
	echo "$(tput setaf 1)Protein-domain prediction by InterProScan5: $(tput setaf 2)${dir_out_reassembly}/prot.faa $(tput sgr 0)"
	sed -i 's/*//g' ${dir_out_reassembly}/prot.faa
	python ${dir_split_fasta} \
		-i ${dir_out_reassembly}/prot.faa \
		-o ${dir_out_reassembly}/${dir_out_anno}/temp_fasta \
		--node ${num_node} \
		--threads ${cpu_max} \
		--memory ${memory_max} \
		--log_out ${dir_out_reassembly}/${dir_out_anno}/log_interproscan/coassembly.prot.faa \
		--log_err ${dir_out_reassembly}/${dir_out_anno}/log_interproscan/coassembly.prot.faa \
		--interproscan ${dir_InterProScan} \
		--out_interproscan_result ${dir_out_reassembly}/${dir_out_anno}/protein_domain_prediction/
		
	# # sbatch runing on each node
	work_dir=${dir_out_reassembly}/${dir_out_anno}/temp_fasta/prot
	list_sh_file=();
	for sh_file in `find $work_dir -type f -name run_interpro*`; do
		sbatch ${sh_file}
	done



	#mapped read of example to binned contig
	# mkdir -p ${dir_out_reassembly}/mapped_read
	mkdir -p "${dir_out_reassembly}/mapped_read/ref"

	# cp ${dir_out_reassembly} ${dir_out_reassembly}/mapped_read/ref
	# bowtie2-build ${dir_out_reassembly}/mapped_read/ref/contig.fa contig --threads 24

	module load samtools
	echo "$(tput setaf 1)Mapping reads$(tput sgr 0)"
	num_lib=${#list_fa_file[@]}
	for (( i=0; i<${num_lib}; i++ )); do
		echo "$(tput setaf 1)Mapping read: $(tput setaf 2)${list_fa_short_name[i]}$(tput sgr 0)"
		# bowtie2 -x ${dir_out_reassembly}/mapped_read/ref/contig \
		# 	-1 ${dir_out_fread}/${list_fa_short_name[i]}_R1.paired.fq.gz \
		# 	-2 ${dir_out_fread}/${list_fa_short_name[i]}_R2.paired.fq.gz \
		# 	-U ${dir_out_fread}/${list_fa_short_name[i]}.unpaired.fq.gz \
		# 	-S ${dir_out_reassembly}/mapped_read/${list_fa_short_name[i]}.sam \
		# 	--threads ${cpu_max}
		# samtools view -bS -@ ${cpu_max} ${dir_out_reassembly}/mapped_read/${list_fa_short_name[i]}.sam | samtools sort -@ ${cpu_max} -m 24G - ${dir_out_reassembly}/mapped_read/${list_fa_short_name[i]}
		# samtools index ${dir_out_reassembly}/mapped_read/${list_fa_short_name[i]}.bam
	done
	echo "Generating visualization"
	# python src/gene_position_on_contig.py ${dir_out_reassembly}/prot.faa.log > ${dir_out_reassembly}/prot.faa.pos
	# anvi-gen-contigs-database -f ${dir_out_reassembly}/mapped_read/ref/contig.fa -o ${dir_out_reassembly}/contig.db --skip-gene-calling

fi

# # assembly
