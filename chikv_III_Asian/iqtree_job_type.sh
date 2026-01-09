#!/bin/bash
#SBATCH --job-name=chikv_iqtree
#SBATCH --mail-user=dongw21@hku.hk     #    Email address to receive notification
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=33         # 申请 33 个 CPU 线程
#SBATCH --mem=32G                  # 内存（根据情况调）
#SBATCH --time=6-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --output=/scr/u/dongw21/Chikungunya/chikv_iqtree%j.out           
#SBATCH --error=/scr/u/dongw21/Chikungunya/chikv_iqtree%j.err           
cd /scr/u/dongw21/Chikungunya
# 如需 module，请按你们系统习惯加载
module load IQ-TREE

#iqtree2 -s chikv_aln.fasta -m MFP -nt 16
# 如果之后要 bootstrap，可以改成：
# iqtree2 -s chikv_aln.fasta -m MFP -nt 16 -B 100
# 3. I-WestAfrica
cd /scr/u/dongw21/Chikungunya/chikv_I_WestAfrica
iqtree2 \
  -s /scr/u/dongw21/Chikungunya/chikv_I_WestAfrica/West_African_aln.fasta \
  -m MFP \
  -nt AUTO  \
  -pre West_African \
  -bb 1000 \
  -alrt 1000

cd /scr/u/dongw21/Chikungunya/chikv_III_Asian
# 2. III-Asian
iqtree2 \
  -s /scr/u/dongw21/Chikungunya/chikv_III_Asian/Asian_aln.fasta \
  -m MFP \
  -nt AUTO  \
  -pre Asian \
  -bb 1000 \
  -alrt 1000

cd /scr/u/dongw21/Chikungunya/chikv_II_ECSA
# 1. II-ECSA
iqtree2 \
  -s /scr/u/dongw21/Chikungunya/chikv_II_ECSA/ECSA_aln.fasta \
  -m MFP \
  -nt AUTO \
  -pre ECSA \
  -bb 1000 \
  -alrt 1000

