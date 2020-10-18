#!/bin/bash
#SBATCH --account=def-sheila
snakemake --profile slurm --config paired='True'
