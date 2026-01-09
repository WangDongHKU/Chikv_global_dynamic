#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH --mail-user=dongw21@hku.hk     #    Email address to receive notification
#SBATCH --output=mafft.%j.out
#SBATCH --error=mafft.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=17         # 申请 17 个 CPU 线程
#SBATCH --mem=10G                  # 内存（根据情况调）
#SBATCH --time=2:00:00            # 最长运行时间 24 小时
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --output=/scr/u/dongw21/Chikungunya/Genebank/mafft%j.out           
#SBATCH --error=/scr/u/dongw21/Chikungunya/Genebank/mafft%j.err   


module load cmake
module load gcc
module load R/4.4.3
module load MAFFT 2>/dev/null || true
mafft --auto --thread 16 \
  /scr/u/dongw21/Chikungunya/chikv_I_WestAfrica/West_African.fasta \
  > /scr/u/dongw21/Chikungunya/chikv_I_WestAfrica/West_African_aln.fasta

mafft --auto --thread 16 \
  /scr/u/dongw21/Chikungunya/chikv_II_ECSA/ECSA.fasta \
  > /scr/u/dongw21/Chikungunya/chikv_II_ECSA/ECSA_aln.fasta


mafft --auto --thread 16 \
  /scr/u/dongw21/Chikungunya/chikv_III_Asian/Asian.fasta\
  > /scr/u/dongw21/Chikungunya/chikv_III_Asian/Asian_aln.fasta
