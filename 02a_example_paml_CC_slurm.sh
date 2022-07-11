#!/bin/bash

#SBATCH --account=def-cptol
#SBATCH --time=96:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=1000M
#SBATCH --mail-user=wesley.gerelle@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=monocot_mht_paml_branch_f61
#SBATCH --output=%x-%j.out

module load python/2.7
module load nixpkgs/16.09
module load intel/2018.3
module load paml/4.9h
source /home/wesley/project/wesley/selection-tests/paml/paml_venv/bin/activate
source /home/wesley/.bashrc

scontrol show hostname ${SLURM_JOB_NODELIST} > ./node_list_${SLURM_JOB_ID}

env_parallel --jobs 48 -u --sshloginfile ./node_list_${SLURM_JOB_ID} --workdir $PWD --joblog monocot_mht_paml_branch_f61.log "cd {} && python ../../../scripts/paml.py *.fasta ../*.tre branch ../test_taxa_*.txt" ::: monocot-mht/*/*/
