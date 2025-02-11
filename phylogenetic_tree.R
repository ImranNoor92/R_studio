library(ape)
library(phangorn)
alignment <- read.phyDat("all_proteins_aligned.fasta", format = "fasta", type = "AA")
dist_matrix <- dist.ml(alignment, model = "JTT")
tree_nj <- NJ(dist_matrix)
plot(tree_nj, main = "Simplified NJ Tree", show.tip.label = TRUE)
# 
# 1. Collapse Very Short Branches (Highly Similar Sequences)
# Use di2multi() from ape to remove very short branches (effectively merging them into a single node).
#library(ape)
tree_simplified <- di2multi(tree_nj, tol = 0.01)  # Adjust tol for merging sensitivity
plot(tree_simplified, cex = 0.7)  # Re-plot the simplified tree
# 2. Collapse Large Clades into One Label
# Use collapse.singles() to simplify clades into single nodes.

tree_collapsed <- collapse.singles(tree_nj)
plot(tree_collapsed, cex = 0.7)
#🔹 This reduces small branches that do not have strong separation.

install.packages("phytools")
library(phytools)
tree_clustered <- collapseTree(tree_nj, nodes = c(5, 10, 20))  # Choose key nodes
plot(tree_clustered, cex = 0.7)



library(ape)
dist_matrix <- cophenetic(tree_nj)  # Get pairwise distance matrix
keep <- colnames(dist_matrix)[apply(dist_matrix, 2, function(x) min(x[x > 0]) > 0.1)]
tree_pruned <- drop.tip(tree_nj, setdiff(tree_nj$tip.label, keep))  
plot(tree_pruned, cex = 0.7)


library(ggtree)
ggtree(tree_nj) + 
  geom_tiplab(size = 3) + 
  geom_cladelab(node = 10, label = "Capsid Proteins", offset = 0.02, align = TRUE) +
  theme_tree2()

plot(tree_nj, type = "fan", cex = 0.7)
ggtree(tree_nj, layout = "circular") + geom_tiplab(size = 3)
ggtree(tree_nj) + geom_tiplab(aes(color = group), size = 3)
ggtree(tree_nj) + geom_tiplab(aes(color = group), size = 3) + theme_tree2()
