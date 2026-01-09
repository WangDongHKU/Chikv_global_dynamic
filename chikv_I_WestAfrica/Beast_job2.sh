#!/bin/bash
#SBATCH --job-name=BEAST2_1
#SBATCH --mail-user=dongw21@hku.hk     #    Email address to receive notification
#SBATCH --cpus-per-task=8         # 申请 8 个 CPU 线程
#SBATCH --mem=32G                  # 内存（根据情况调）
#SBATCH --time=5-10:00:00            # 最长运行时间 24 小时
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --output=/scr/u/dongw21/Chikungunya/chikv_I_WestAfrica/Beast_job_%j.out           
#SBATCH --error=/scr/u/dongw21/Chikungunya/chikv_I_WestAfrica/Beast_job_%j.err           

cd /scr/u/dongw21/Chikungunya/chikv_I_WestAfrica/

# 如需 module，请按你们系统习惯加载
module load cuda
module load cuda-toolkit
module load BEAGLE/4.0.0/amd

export PATH=$HOME/beast/bin:$PATH
which beast
beast -version


beast -beagle_info

beast \
 -beagle \
 -threads -1 \
 -overwrite West_African_aln_2.xml



