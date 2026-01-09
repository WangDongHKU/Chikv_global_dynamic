#!/bin/bash
#SBATCH --job-name=BEAST1_2
#SBATCH --mail-user=dongw21@hku.hk     #    Email address to receive notification
#SBATCH --cpus-per-task=20         # 申请 20 个 CPU 线程
#SBATCH --mem=32G                  # 内存（根据情况调）
#SBATCH --time=5-10:00:00            # 最长运行时间 24 小时
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --output=/scr/u/dongw21/Chikungunya/chikv_II_ECSA/Beast_job_%j.out           
#SBATCH --error=/scr/u/dongw21/Chikungunya/chikv_II_ECSA/Beast_job_%j.err           

cd /scr/u/dongw21/Chikungunya/chikv_II_ECSA/
# 如需 module，请按你们系统习惯加载
module load cuda
module load cuda-toolkit
module load BEAGLE/4.0.0/amd
module load BEAST
beast -beagle_info

NTHREADS=20   # 和 --cpus-per-task 保持一致

beast \
 -beagle \
 -beagle_CPU \
 -threads $NTHREADS \
 -beagle_thread_count $NTHREADS\
 -overwrite ECSA_aln_subsampled_1.xml

