 import re
 import pandas as pd
 path = './data/ref/annotations/Homo_sapiens.GRCh38.Ensembl86.gtf'
 pattern ='gene_id "(ENSG[0-9]+)";.+ gene_name "(.+)";.+ gene_biotype "protein_coding";.+' 
 with open(path,'r') as f: 
          protein_coding = [] 
          for line in f: 
              m = re.search(pattern, line) 
              if m: 
                 protein_coding.append(','.join([m.group(1),m.group(2)])) 
    reduced = list(set(protein_coding))
id = [x.split(',')[0] for x in reduced] 
name = [x.split(',')[1] for x in reduced] 
d = {} 
d['gene_id'] = id 
d['gene_name'] = name 
df = pd.DataFrame(d) 
df.to_csv('./data/ref/annotations/protein_coding.txt', header=True, index=False, sep='\t', mode='a') 