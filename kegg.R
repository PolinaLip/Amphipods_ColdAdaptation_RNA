
### MAPK signaling pathway, Gla: 
N = 4710 # number of genes, annotated by KEGG
K = 71 # number of genes found in pathway for the species
n = 1726 # number of DE genes (0.001)
k = 50 # number of DE genes in the pathway
phyper(k - 1, K,  N - K, n, lower.tail = F) # hypergeometric test
# 5.992426e-09 - p-value

### ribosome biogenesis
N = 4710 # number of genes, annotated by KEGG
K = 54 # number of genes found in pathway for the species
n = 1726 # number of DE genes (0.001)
k = 35 # number of DE genes in the pathway
phyper(k - 1, K,  N - K, n, lower.tail = F) # hypergeometric test
# 2.21319e-05

### regulation of actin cytoskeleton
N = 4710 # number of genes, annotated by KEGG
K = 51 # number of genes found in pathway for the species
n = 1726 # number of DE genes (0.001)
k = 27 
phyper(k - 1, K,  N - K, n, lower.tail = F)
# 0.01235189

