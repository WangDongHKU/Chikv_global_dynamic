#!/bin/bash
#SBATCH --job-name=BEAST2_3
#SBATCH --mail-user=dongw21@hku.hk     #    Email address to receive notification
#SBATCH --output=BEAST.%j.out
#SBATCH --error=BEAST.%j.err
#SBATCH --cpus-per-task=4         # 申请 8 个 CPU 线程
#SBATCH --mem=16G                  # 内存（根据情况调）
#SBATCH --time=1-10:00:00            # 最长运行时间 24 小时
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --output=/scr/u/dongw21/Chikungunya/chikv_III_Asian/Beast_job_%j.out           
#SBATCH --error=/scr/u/dongw21/Chikungunya/chikv_III_Asian/Beast_job_%j.err           

cd /scr/u/dongw21/Chikungunya/chikv_III_Asian

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
 -overwrite Asian_2.xml

#beast -beagle -overwrite Asian.xml
#module load BEAST
#beast -beagle_GPU -overwrite /scr/u/dongw21/Chikungunya/chikv_III_Asian/Asian1.xml

#beast -beagle -overwrite /scr/u/dongw21/Chikungunya/chikv_III_Asian/Asian1.xml
