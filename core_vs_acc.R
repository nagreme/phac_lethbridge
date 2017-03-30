#Make a heatmap of the accessory difference organized by core clustering (tree)

#Nadège Pulgar-Vidal
#March 6, 2107

# ****************
# Note: There are functions at the bottom of the file 
# ****************

library(ape)
library(dendextend)
library(dendextendRcpp)
library(phangorn)
library(gplots)
library(RColorBrewer)
library(prabclus)
library(Matrix)
library(magrittr)

#/home/nadege/Desktop/mist/108_and_97/core_json_out.csv
#/home/nadege/Desktop/mist/108_and_97/acc_json_out.csv
#/home/nadege/Desktop/merging/cd-hit-est_attempt/merge_process_100/control/gene_presence_absence_accessory_97_108.Rtab

#/home/nadege/Desktop/mist/big_set/core_json_out.csv
#/home/nadege/Desktop/mist/big_set/acc_json_out.csv
#/home/nadege/Desktop/mist/big_set/acc_presence_absence.Rtab

core_file <- "/home/nadege/Desktop/mist/big_set/core_json_out.csv"
acc_file <- "/home/nadege/Desktop/mist/big_set/acc_json_out.csv"
binary_acc_file <- "/home/nadege/Desktop/mist/big_set/acc_presence_absence.Rtab"


#prep core genes
cg_dendro <- read.table(file=core_file, row.names = 1, header=T, sep=",") %>%
dist.gene(method = "pairwise", pairwise.deletion = FALSE, variance = FALSE) %>%
hclust(method = "complete") %>%
as.dendrogram()

#prep binary acc genes
bin_acc_genes <- read.csv(file=binary_acc_file, row.names=1, header=T, sep=",") 
colnames(bin_acc_genes) <- sub("X","", colnames(bin_acc_genes)) #Can't really pipe these because of this step
dist_bin_acc <- jaccard(as.matrix(bin_acc_genes))
#hc_bin_acc <- hclust(as.dist(dist_bin_acc), method = "complete") #these two are needed for the tanglegram and other dendros
#bin_acc_dendro <- as.dendrogram(hc_bin_acc)

#prep allelic acc genes
acc_genes <- read.table(file=acc_file, row.names = 1, header=T, sep=",")
colnames(acc_genes) <- sub("X","", colnames(acc_genes))
acc_genes[acc_genes == 0] <- NA 
#dist_acc <- dist.gene(acc_genes, method = "pairwise", pairwise.deletion = T, variance = F)
#hc_acc <- hclust(dist_acc, method = "complete")
#acc_dendro <- as.dendrogram(hc_acc)


#make a heatmap that uses core clustering but accessory differences on heatmap
heatmap.2(as.matrix(dist_bin_acc),Rowv = cg_dendro, Colv = "Rowv", distfun=return_same, dendrogram = "both", symm = T,
          trace = "none", col = brewer.pal(9,"PuBu"), main = "Accessory Distance Matrix Clustered by Core",
          key.title = "Accessory Distance Distribution", key.xlab = "Accessory Distance")
#cexRow = 0.35, cexCol = 0.35, offsetRow = 0.001, offsetCol = 0.001)

#tanglegram
#cg_vs_acc_dend <- dendlist(cg_dendro, bin_acc_dendro)
#untangled_dendro <- untangle(cg_vs_acc_dend,method = "step2side")
#tangle <- round(entanglement(untangled_dendro),3)
#tanglegram(untangled_dendro,common_subtrees_color_branches = TRUE, margin_inner = 4.5, lab.cex = 0.8, margin_outer = 2,
#           main="Core vs. Accessory",main_left="Core",main_right = "Accessory (binary)",sub=paste("entanglement = ",tangle))


#Scatterplot of binary acc dist (x axis) vs allelic acc dist (y axis, pariwise deletion on)
#Convert distances to similarities
acc_sim <- 1 - dist.gene(acc_genes, method = "percentage", pairwise.deletion = T, variance = F) %>% as.matrix()
bin_acc_sim <- 1 - dist_bin_acc 

acc_sim <- tril(acc_sim, k=-1)
bin_acc_sim <- tril(bin_acc_sim, k=-1)
acc_sim[acc_sim == 0] <- NA 
bin_acc_sim[bin_acc_sim == 0] <- NA 
acc_sim[acc_sim == 0] <- NA 
bin_acc_sim[bin_acc_sim == 0] <- NA

plot(x = bin_acc_sim, y = acc_sim, type = "p", main = "Gene Presence Absence vs. Allele Similarity in Accessory",
     xlab = "Overlapping Loci (Jaccard similarity)", ylab = "Allele Similarity", xlim = c(0,1), ylim = c(0,1))


#Histogram of # of shared loci
acc_counts <- apply(bin_acc_genes,2,sum)
hist(acc_counts, xlim= c(600,1400), main = "Number of Accessory Genes per Genome", xlab="# of Accessory Genes",
     ylab = "# of Genomes")



return_same <- function(x)
{
  as.dist(x)
}




#Read in the data
#core_genes <- read.table(file=core_file, row.names = 1, header=T, sep=",")
#acc_genes <- read.table(file=acc_file, row.names = 1, header=T, sep=",")
#colnames(acc_genes) <- sub("X","", colnames(acc_genes))
#bin_acc_genes <- read.csv(file=binary_acc_file, row.names=1, header=T, sep=",") 
#colnames(bin_acc_genes) <- sub("X","", colnames(bin_acc_genes))
#bin_acc_genes <- bin_acc_genes[,-1] #remove metagenome column

#Cluster/distance matrix calculations
#dist_cg <- dist.gene(core_genes, method = "pairwise", pairwise.deletion = FALSE, variance = FALSE)

#For the accessory either 0->NA and pairwise deletion with dist.gene, or jaccard, or something similar?
#acc_genes[acc_genes == 0] <- NA 
#dist_acc <- dist.gene(acc_genes, method = "pairwise", pairwise.deletion = T, variance = F)
#dist_acc_no_del <- dist.gene(acc_genes, method = "pairwise", pairwise.deletion = F, variance = F)

#The binary accessory uses jaccard distance because it's a very sparse matrix
#dist_bin_acc <- jaccard(as.matrix(bin_acc_genes))

#Cluster the data
#hc_cg <- hclust(dist_cg, method = "complete")
#hc_acc <- hclust(dist_acc, method = "complete")
#hc_acc_no_del <- hclust(dist_acc, method = "complete")
#hc_bin_acc <- hclust(as.dist(dist_bin_acc), method = "complete")

#Make dendrograms out of the clustering
#cg_dendro <- as.dendrogram(hc_cg)
#acc_dendro <- as.dendrogram(hc_acc)
#acc_no_del_dendro <- as.dendrogram(hc_acc_no_del)
#bin_acc_dendro <- as.dendrogram(hc_bin_acc)







