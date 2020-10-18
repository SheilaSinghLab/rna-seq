#!/bin/bash
#SBATCH --account=def-sheila
snakemake qc_all --profile slurm
