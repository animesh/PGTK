#!/bin/bash
#SBATCH --account=nn9036k
#SBATCH --job-name=resultsTKvep
#SBATCH --time=168:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --mail-user=animesh.sharma@ntnu.no
#SBATCH --mail-type=ALL
#SBATCH --output=resultsTKvep-%j.log
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd -P)}"
cd "${PROJECT_DIR}"
module purge
module load dotNET-Core/8.0.411
dotnet --version
#8.0.411
#wget "https://maxquant.org/p/maxquant/MaxQuant_v2.8.1.0.zip?md5=73u6i7B82tlhSVoJuSo6hQ&expires=1785704833"
#unzip -o MaxQuant_v2.8.1.0.zip
#wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/"
dotnet /cluster/home/ash022/scripts/MaxQuant_v2.8.1.0/bin/MaxQuantCmd.dll  /cluster/home/ash022/scripts/pgtk/mqpar.xml
#/cluster/projects/nn9036k/dotnet80303/dotnet
tail -n 20 resultsTKvep-${SLURM_JOB_ID}.log
