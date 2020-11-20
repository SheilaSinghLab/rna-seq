## -----------------------------------------------------------------------------------------------------------------------
load_data = function(path){
  data = read.csv(path, sep='\t', header=T, stringsAsFactors = F)
  data = clean(data)
}
clean = function(data){
  #filter NA rows
  data = data[!apply(is.na(data[,-1:-2,drop=F]),1,all),]
} 

make_dge = function(data, design=NULL, filt=T, norm=T){
  require(edgeR)
  dge = DGEList(counts=data[,-1:-2], genes=data[,1:2])
  if (filt){
    keep = filterByExpr(dge, design)
    dge = dge[keep,,keep.lib.sizes=F]
  if (norm){
    dge = calcNormFactors(dge)
  }
  }
}

#if multiple replicates
run_dge_replicates = function(dge,design, num_genes, bayes=F){
  require(limma)
  v = voom(dge, design)
  fit = lmFit(v, design)
  if(bayes){
    fit_bayes = eBayes(fit)
    hits =topTable(fit_bayes, coef=ncol(design), number=num_genes)
  } else {
    fit_treat = treat(fit, lfc=log2(1), trend=T)
    hits = topTreat(fit_treat,coef=ncol(design),number = num_genes)
  }
  
 
}

args = commandArgs(trailingOnly=T)
data_path = args[1] #path to counts.txt file
design_path = args[2] #path to design matrix indicating replicates 
num_genes = as.numeric(args[3]) #num of top differential genes to display
bayes = as.logical(args[4]) #implement bayes DEG or not (Boolean)
output= args[5] #where to save DEG output

data = load_data(data_path)
design = scan(design_path,sep=',')
design = factor(design)
design = model.matrix(~design)

dge = make_dge(data, design)
hits = run_dge_replicates(dge, design, num_genes, bayes)
write.table(hits, output,sep='\t', row.names=F,quote=F)

