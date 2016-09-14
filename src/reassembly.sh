#!/bin/bash
#SBATCH --job-name=CoAs2Merge
#SBATCH --partition=largemem
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=500G
#SBATCH --ntasks=1
#SBATCH --mail-user=nattawet.sriwichai@oist.jp
#SBATCH --mail-type=END
#SBATCH --input=none
#SBATCH --output=job_reassembly_%j.out
#SBATCH --error=job_reassembly_%j.err

# Shotgun metagenomic analysis pipeline including 6 step of analysis processes
# 0) Checking qualitiy of seqeuncing data
# 1) Filtering good quality reads of each metatrascriptomic
# 2) Assigning taxonomic of each metagenomic data
# 3) Assembling read to contig and scaffold
# 4) Calling gene and extract protein seqeunce from assembled read
# 5) Annoting protein sequence
# 6) Comparing fuctional of each metagenomic


/apps/unit/GoryaninU/software/bining/MaxBin-2.2/run_MaxBin.pl \
-reads /work/GoryaninU/golf/1_seq_filtered/merge_R1.paired.fq.gz \
-reads2 /work/GoryaninU/golf/1_seq_filtered/merge_R2.paired.fq.gz \
-contig /work/GoryaninU/golf/3_coassembled_2/contigs.fasta \
-thread 24 \
-out /work/GoryaninU/golf/3_bining \
-reassembly

