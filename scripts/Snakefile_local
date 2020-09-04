# getting command line args for the method of normalization (default = tmm)
# an example command for running snakemake with custom arguments: snakemake --config nm="tmm"
try:
    normMethod = config['nm'].lower()
    if normMethod != "tmm" and normMethod != "rpkm":
        normMethod = "tmm"
except KeyError:
    normMethod = "tmm"

# creating a global wildcard to account for all count_matrix.txt files
COUNTS = glob_wildcards('counts/{count}.txt').count

localrules: all, normalize

# current final output is a normalized count matrix for each .txt file
rule all:
    input:
         counts = expand('norm_counts/{count}_norm.txt', count = COUNTS)

# rule that takes all counts and produces normalized counts by TPM
rule normalize:
    input:
         counts = 'counts/{file}.txt',
         rFile = 'normalize.R'
    output: 'norm_counts/{file}_norm.txt'
    shell:
        '''
        echo normalizing {input.counts} to create {output}
        Rscript {input.rFile} {normMethod} {input.counts} {output}
        '''