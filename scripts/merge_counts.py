import sys
import re
import pandas as pd
import os
num_files = len(sys.argv)
files = sys.argv[2:]
pattern = 'data/aligned/(.*)_ReadsPerGene.out.tab'
path = './data/ref/annotations/protein_coding.txt'
counts = {}
genes = []
samples=[]
for f in files:
    match = re.search(pattern, f)
    sample = match.group(1)
    samples.append(sample)
    counts[sample] = []
    with open(f, 'r') as f:
        for line in f:
            vals = line.rstrip().split()
            counts[sample].append(vals[3])
            genes.append(vals[0])
        num_genes = len(counts[sample])
counts['genes'] = genes[:num_genes]
df = pd.DataFrame(counts)
df = df[['genes']+samples]
df_pc = pd.read_csv(path,sep='\t')
merged_inner = pd.merge(left=df, right=df_pc, how='right',left_on='genes', right_on='gene_id')
merged = merged_inner[['gene_id','gene_name']+samples]

if os.path.isfile(sys.argv[1]):
    existing = pd.read_csv(sys.argv1, sep='\t')
    merged = pd.merge(left=existing, right=merged[['gene_id']+samples], how='right',left_on='gene_id',right_on='gene_id')
    merged.to_csv(sys.argv[1],header=True,index=False,sep='\t',mode='a')
else:
    merged.to_csv(sys.argv[1], header=True, index=False, sep='\t', mode='a')