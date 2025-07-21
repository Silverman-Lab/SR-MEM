#!/bin/bash

#SBATCH --account=open          # specify the account
#SBATCH --partition=open        # specify the partition
#SBATCH --job-name=deda
#SBATCH --ntasks=2
#SBATCH --time=7:00:00
#SBATCH --mem=24GB

Rscript run_mems.R $1 $2
