#!/bin/bash
#SBATCH --account=def-sheila
snakemake metrics_all --profile slurm
