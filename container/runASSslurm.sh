#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=output.txt
#SBATCH --job-name=ass

# docker run did not work, especially not on slurm
# docker run --user $(id -u):$(id -g) -v /nfs/proj/Sys_CARE/AS_Simulator/ensembl_data/Homo_sapiens.GRCh38.99:/input -v /nfs/proj/Sys_CARE/AS_Simulator/ensembl_data/Homo_sapiens.GRCh38.99_output:/output biomedbigdata/ass
# do singularity run instead but first create image file: 
# singularity build ass_from_docker.sif docker://quirinmanz/ass_docker:latest
singularity run -B /nfs/proj/Sys_CARE/AS_Simulator/ensembl_data/Homo_sapiens.GRCh38.99:/input -B /nfs/proj/Sys_CARE/AS_Simulator/ensembl_data/Homo_sapiens.GRCh38.99_output:/output ass_from_docker.sif
