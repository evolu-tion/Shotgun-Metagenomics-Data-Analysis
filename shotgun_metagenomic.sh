#!/bin/bash
#SBATCH --job-name=metagenomic_analysis
#SBATCH --partition=compute
#SBATCH --time=84:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=125G
#SBATCH --ntasks=1
#SBATCH --input=none
#SBATCH --output=metagenomic_analysis_%j.out
#SBATCH --error=metagenomic_analysis_%j.err

# Shotgun metagenomic analysis pipeline including 6 step of analysis processes
# 0) Checking qualitiy of seqeuncing data
# 1) Filtering good quality reads of each metatrascriptomic
# 2) Assigning taxonomic of each metagenomic data
# 3) Assembling read to contig and scaffold
# 4) Calling gene and extract protein seqeunce from assembled read
# 5) Annoting protein sequence
# ----6) Comparing fuctional of each metagenomic

echo -e "=================================================================== 
              Starting metagenomics data analysis 
==================================================================="

# Load configuration file
source configuration
dir_out_qc=${dir_out}/${dir_out_qc}
dir_out_fread=${dir_out}/${dir_out_fread}
dir_out_taxa=${dir_out}/${dir_out_taxa}
dir_out_assem=${dir_out}/${dir_out_assem}
dir_out_bining=${dir_out}/${dir_out_bining}
dir_out_cgene=${dir_out}/${dir_out_cgene}
dir_out_anno=${dir_out}/${dir_out_anno}



# # Load fastq file
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

# list_fa_file=("GO-EL-03_S3_L001_R1_001.fastq.gz,GO-EL-03_S3_L001_R1_001.fastq.gz")
# list_fa_short_name=("GOEL03S3L001")
# list_fa_file=("${list_fa_file[@]:1}")

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
if [ "${step_3_assembling}" == "Y" ]; then
	echo "+ Run step 3: Assembly metagneomic"
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
			ILLUMINACLIP:src/adapters/TruSeq3-PE.fa:2:30:10 \
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

# # STEP3: Assembling read to contig and scaffold
if [ "${step_3_assembling}" == "Y" ]; then
	echo "$(tput setab 7)$(tput setaf 1)========== STEP 3: Assembling read to contig and/or scaffold ==========$(tput sgr 0)"
	mkdir -p ${dir_out_assem}
	contig_file_name=""

	num_lib=${#list_fa_file[@]}
	for (( num=0; num<${num_lib}; num++ )); do
		echo "$(tput setaf 1)Assembling by metaSPAdes: $(tput setaf 2)${list_fa_short_name[num]}$(tput sgr 0)"
		${dir_SPAdes} \
			--meta \
			--threads ${cpu_max} \
			--pe1-1 ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz \
			--pe1-2 ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz \
			--pe1-s ${dir_out_fread}/${list_fa_short_name[num]}.unpaired.fq.gz \
			-o ${dir_out_assem}/${list_fa_short_name[num]}
		
		# Change heading of each contig and/or merge contig with scaffold
		python src/merge_contig_scaffold.py \
			--contig ${dir_out_assem}/${list_fa_short_name[num]}/contigs.fasta \
			--output ${dir_out_assem}/${list_fa_short_name[num]}.fa \
			--sample ${list_fa_short_name[num]} 

		contig_file_name="${contig_file_name}${dir_out_assem}/${list_fa_short_name[num]}.fa "
	done

	# # STEP3.1: Check quality of metagenomic assembly by using metaQuast
	echo "$(tput setaf)Check quality of metagenomic assembly by metaQUAST $(tput sgr 0)"
	${dir_metaQuast} ${contig_file_name} -o ${dir_out_assem}/quality --threads=${cpu_max} --max-ref-number 1000

fi


# # STEP4: Calling gene and extract protein seqeunce from assembled read
if [ "${step_4_gene_calling}" == "Y" ]; then
	echo "$(tput setab 7)$(tput setaf 1)=== STEP 4: Calling gene and extract protein seqeunce from assembled read ===$(tput sgr 0)"
	mkdir -p ${dir_out_cgene}

	num_lib=${#list_fa_file[@]}
	for (( num=0; num<${num_lib}; num++ )); do
		echo "$(tput setaf 1)Gene calling by FragGeneScan: $(tput setaf 2)${list_fa_short_name[num]}$(tput sgr 0)"
		${dir_FragGeneScan} -s ${dir_out_assem}/${list_fa_short_name[num]}.fa \
			-o ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene \
			-w 0 \
			-t illumina_5 \
			-p ${cpu_max}

		# Change fasta head name
		python src/change_gene_name.py -i ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene.faa \
			-a ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene.out \
			-o ${dir_out_cgene}/${list_fa_short_name[num]}.prot.faa &

		# Change nucleotide fasta name
		python src/change_gene_name.py -i ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene.ffn \
			-a ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene.out \
			-o ${dir_out_cgene}/${list_fa_short_name[num]}.nucl.fa
	done
fi

# # STEP5: Annotating protein sequence divied into 2 part
if [ "${step_5_annotating}" == "Y" ]; then
	echo "$(tput setab 7)$(tput setaf 1)========== STEP 5: Annotating protein sequence ==========$(tput sgr 0)"
	mkdir -p ${dir_out_anno}
	num_lib=${#list_fa_file[@]}
	# PART A Annotating of each protein by GhostKOALA (http://www.kegg.jp/ghostkoala/) by manually
	#        Input file is protein sequence from gene calling tool in directory "${dir_out_cgene}" in ${entry_name}-calling_gene-contigs.faa
	#        Output download from websites

	# PART A Annotating of each nucleotide sequece by used multiple database KEGG, Pfam, GO, and SEED sub system information from the COG annotations.
	#        by useing COGNIZER tool
	mkdir -p ${dir_out_anno}/cognizer
	for (( num=0; num<${num_lib}; num++ )); do
		echo "$(tput setaf 1)Gene functional annotation by COGNIZER: $(tput setaf 2)${list_fa_short_name[num]}$(tput sgr 0)"
		cd ${dir_cognizer}
		mkdir ${dir_out_anno}/cognizer/${list_fa_short_name[num]}
		./cognizer -t ${cpu_max} -e ${p_cognizer_evalue} -i ${dir_out_cgene}/${list_fa_short_name[num]}.nucl.fa -o ${dir_out_anno}/cognizer/${list_fa_short_name[num]}
	done

	# # PART B Annotating protein domain by InterProScan
	# echo "$(tput setab 7)$(tput setaf 2)Protein domain predicting by InterProScan5 $(tput sgr 0)"
	# mkdir -p ${dir_out_anno}/temp_fasta/
	# mkdir -p ${dir_out_anno}/protein_domain_prediction/
	# mkdir -p ${dir_out_anno}/log_interproscan

	num_lib=${#list_fa_file[@]}
	for (( num=0; num<${num_lib}; num++ )); do
		# Split fasta file
		echo "$(tput setaf 1)Protein-domain prediction by InterProScan5: $(tput setaf 2)${list_fa_short_name[num]}$(tput sgr 0)"
		sed -i 's/*//g' ${dir_out_cgene}/${list_fa_short_name[num]}.prot.faa
		python ${dir_split_fasta} \
			-i ${dir_out_cgene}/${list_fa_short_name[num]}.prot.faa \
			-o ${dir_out_anno}/temp_fasta \
			--node ${num_node} \
			--threads ${cpu_max} \
			--memory ${memory_max} \
			--log_out ${dir_out_anno}/log_interproscan/${list_fa_short_name[num]}.prot.faa \
			--log_err ${dir_out_anno}/log_interproscan/${list_fa_short_name[num]}.prot.faa \
			--interproscan ${dir_InterProScan} \
			--out_interproscan_result ${dir_out_anno}/protein_domain_prediction/
		
		# # sbatch runing on each node
		work_dir=${dir_out_anno}/temp_fasta/${list_fa_short_name[num]}
		list_sh_file=();
		for sh_file in `find $work_dir -type f -name run_interpro*`; do
			sbatch ${sh_file}
		done
	done
fi


# # STEP6: Functional comparative between microbial communities
# $ Be constructing

## Display the result
if [ "${display_result}" == "Y" ]; then
	declare -A info
	num_lib=${#list_fa_file[@]}
	info[0,0]="Libary_id"
	info[1,0]="# paired-end read"
	info[2,0]="# paired read good qualtiy"
	info[3,0]="# forward read good quality"
	info[4,0]="# backward read good quality"
	info[5,0]="# SSU read"
	info[6,0]="# + archea SSU read"
	info[7,0]="# + bacteria SSU read"
	info[8,0]="# contig"
	info[9,0]="# protein-coding gene"
	for ((i=1;i<$num_lib+1;i++)) do
		fastq_paired_end=(${list_fa_file[$(($i-1))]//,/ })
		short_name=${list_fa_short_name[$(($i-1))]}

		# display short name
		info[0,$i]=${short_name}
		echo "Loading libary: ${short_name}"
		# paired-end read in library
		case ${fastq_paired_end[0]} in
			*.gz ) 
				info[1,$i]=$(zcat ${dir_seq}/${fastq_paired_end[0]} | echo $((`wc -l`/4)))
				;;
			* )
				info[1,$i]=$(cat ${dir_seq}/${fastq_paired_end[0]} | echo $((`wc -l`/4)))
				;;
		esac

		# Filtering good quality read in library
		# - Paired good quality
		info[2,$i]=$(zcat ${dir_out_fread}/${short_name}_R1.paired.fq.gz | echo $((`wc -l`/4)))
		
		# - Single forward good quality
		info[3,$i]=$(zcat ${dir_out_fread}/${short_name}_R1.unpaired.fq.gz | echo $((`wc -l`/4)))
		info[4,$i]=$(zcat ${dir_out_fread}/${short_name}_R2.unpaired.fq.gz | echo $((`wc -l`/4)))
		
		# Number of SSU read
		info[5,$i]=$(sed -n -e 's/Total number of classifications made by Metaxa://p' 2_taxonomy/SSU/${short_name}/${short_name}.summary.txt)
	
		# - Archrea SSU read
		# - Bacteria SSU read
		# Number of contig
		info[8,$i]=$(cat ${dir_out_assem}/${short_name}.fa | echo $((`wc -l`/2)))
		
		# Number of protein-coding gene in all contig
		info[9,$i]=$(cat ${dir_out_cgene}/${short_name}.prot.faa | echo $((`wc -l`/2)))

	done

	num_list_name=9
	for ((j=0;j<=${num_list_name};j++)) do
		for ((i=0;i<$num_lib+1;i++)) do
			printf "${info[$j,$i]}\t"
		done
	done

fi
