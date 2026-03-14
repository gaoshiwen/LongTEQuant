#!/bin/bash
#SBATCH --job-name=miniQuant_analysis
#SBATCH --partition=hpxg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --output=miniquant_%j.out
#SBATCH --error=miniquant_%j.err

# 环境设置
cd $SLURM_SUBMIT_DIR
echo "开始时间: $(date)"
echo "作业ID: $SLURM_JOB_ID"
echo "节点: $SLURM_JOB_NODELIST"

# 激活conda环境 (使用您的路径)
echo "激活conda环境..."
eval "$(/home/liaojianshu/project/lijiaproject/miniconda3/bin/conda shell.bash hook)"
conda activate miniquant

# 验证环境
echo "当前Python路径: $(which python)"
echo "当前环境: $CONDA_DEFAULT_ENV"
python --version

# 运行miniQuant
echo "开始运行miniQuant..."
cd /project/liaojianshu/miniQuant_old/isoform_quantification/

python main.py quantify \
    -gtf /project/liaojianshu/test_data/gencode.v40.annotation.exon.only.gtf \
    -srsam /project/liaojianshu/test_data/5ml_40ms/40M.sorted.sam \
    -lrsam /project/liaojianshu/test_data/5ml_40ms/simulated_reads_cONT.minimap2.sorted.bam \
    -o /project/liaojianshu/test_data/miniQuant_output \
    --alpha adaptive \
    --EM_choice hybrid \
    --threads 8

echo "完成时间: $(date)"