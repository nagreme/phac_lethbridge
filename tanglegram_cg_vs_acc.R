#library('ctv') 
#install.views('Phylogenetics')
#update.views('Phylogenetics')
library(ape)
library(dendextend)
library(dendextendRcpp)
library(phangorn)
library(gplots)
library(RColorBrewer)
library(prabclus)

#"/home/nadege/Desktop/mist/108_genomes/json_out.csv" 
#"/home/nadege/Desktop/roary/corrected_prokka_db_full_set_exclude_CP017855/gene_presence_absence_accessory.Rtab" 

#"/home/nadege/Desktop/mist/108_genomes/json_out_exc_61s.csv" 
#"/home/nadege/Desktop/roary/corrected_prokka_db_full_set_exclude_CP017855/gene_presence_absence_accessory_exc_61s.Rtab" 

#"/home/nadege/Desktop/mist/old_97/json_out.csv"
#"/home/nadege/Desktop/roary/old_full_set_97_meta0/gene_presence_absence_accessory.Rtab"

#"/home/nadege/Desktop/mist/108_genomes/json_out.csv"
#"/home/nadege/Desktop/panseq/full_no_d_p/binary_table_fixed_labels_accessory.csv" 

#read/format data frames
core_genes <- read.table(file="/home/nadege/Desktop/mist/108_genomes/json_out.csv", row.names = 1, header=T, sep=",")
acc_genes_t <- read.table(file="/home/nadege/Desktop/panseq/full_no_d_p/binary_table_fixed_labels_accessory.csv", row.names=1, header=T, sep=",")
colnames(acc_genes_t) <- sub("X","", colnames(acc_genes_t)) 
#acc_genes <- t(acc_genes_t) #no need to transpose for jaccard distance
#rownames(acc_genes) <- sub("X","", rownames(acc_genes))

other_acc_genes <- read.table(file="/home/nadege/Desktop/roary/corrected_prokka_db_full_set_exclude_CP017855/gene_presence_absence_accessory.Rtab",row.names=1, header=T, sep=",")
colnames(other_acc_genes) <- sub("X","", colnames(other_acc_genes))

#dist.gene(x, method = "pairwise", pairwise.deletion = FALSE, variance = FALSE)
# x         a matrix or a data frame (will be coerced as a matrix).
# method    a character string specifying the method used to compute the distances; two
#           choices are available: "pairwise" and "percentage", or any unambiguous
#           abbreviation of these.
#pairwise.deletion
#           a logical indicating whether to delete the columns with missing data on a pairwise
#           basis. The default is to delete the columns with at least one missing observation.
#variance   a logical, indicates whether the variance of the distances should be returned
#(default to FALSE).

dist_cg <- dist.gene(core_genes, method = "pairwise", pairwise.deletion = FALSE, variance = FALSE)
dist_acc <- jaccard(as.matrix(acc_genes_t))
dist_other_acc <- jaccard(as.matrix(other_acc_genes))

hc_cg <- hclust(dist_cg, method = "complete")
hc_acc <- hclust(as.dist(dist_acc), method = "complete")
hc_other_acc <- hclust(as.dist(dist_other_acc), method = "complete")

cg_dendro <- as.dendrogram(hc_cg)
acc_dendro <- as.dendrogram(hc_acc)
other_acc_dendro <- as.dendrogram(hc_other_acc)


#One tanglegram
tanglegram(cg_dendro, acc_dendro, common_subtrees_color_branches = TRUE, sort = T, margin_inner = 4,
           main="Core vs. Accessory",main_left="Core",main_right = "Accesory")


#slightly different method and tanglegram

cg_vs_acc_dend <- dendlist(cg_dendro, acc_dendro)
#dend_diff(cg_dendro,acc_dendro)

acc1_vs_acc2_dend <- dendlist(other_acc_dendro, acc_dendro) #this one (other_acc) is Roary currently, watch your labels



#untangle it
#method labels is fast and okay, step1side is better and takes moderate time, step2side is really slow but is the best
untangled_dendro <- untangle(cg_vs_acc_dend,method = "step2side")
#untangled_dendro <- untangle(acc1_vs_acc2_dend,method = "step2side") #use main="Roary acc vs. Panseq acc"


#calc entanglement (lower is better apparently)
tangle <- round(entanglement(untangled_dendro),3)

#plot tanglegram
tanglegram(untangled_dendro,common_subtrees_color_branches = TRUE, margin_inner = 4.5, lab.cex = 0.8, margin_outer = 2,
           main="Accessory vs. Accessory",main_left="Roary",main_right = "Panseq",sub=paste("entanglement = ",tangle))



#computes Robinson-Foulds distance
dist.dendlist(untangled_dendro)
#other
RF.dist(as.phylo(hc_cg), as.phylo(hc_acc), check.labels=TRUE, rooted=FALSE) 
#same result (154)



#cutree and Adjusted Wallace Coefficient (AWC)

max_height = 697
heights = c(0:max_height)
awcs = vector("double",max_height+1)

for (height in c(0:max_height))
{
  cg_groups <- cutree(cg_dendro, h = height )
  acc_groups <- cutree(acc_dendro, h = height ) 
  
  
  try(
    {
      awc <- adj_wallace(cg_groups, acc_groups) #****Error here at certain heights!!! (Not sure why...)
      awcs[height] <- awc[5]
      })  
  
  #if (awc[5] != 1)
  #{
  #str <- sprintf("height: %i", height)
  #print(str)
  #print(awc[5])
  #}
  
}

plot(x=heights, y=awcs, xlabel="height", ylabel="AWC", type="p")

#cg_groups <- cutree(cg_dendro, h = 675 )
#acc_groups <- cutree(cg_dendro, h = 675 )



# Heatmap of Difference Matrix
p_cg_dist <- dist.gene(core_genes, method = "percentage", pairwise.deletion = FALSE, variance = FALSE)
#p_acc_dist <- dist.gene(acc_genes, method = "percentage", pairwise.deletion = FALSE, variance = FALSE)

#scale the acc to range [0-1]
#p_acc_dist_scalled <- linMap(p_acc_dist,0,1)

#heatmap2 applies its own dist function so override it
#run return_same function at bottom of this script
heatmap.2(as.matrix(p_cg_dist), dendrogram="both", trace = "none", col = brewer.pal(9,"PuBu"), main = "Core", distfun=return_same,cexRow = 0.5, cexCol = 0.5)
heatmap.2(as.matrix(dist_acc), dendrogram="both", trace = "none", col = brewer.pal(9,"PuBu"), main = "Panseq Accessory",distfun=return_same,cexRow = 0.5, cexCol = 0.5)
#heatmap.2(as.matrix(p_acc_dist_scalled), dendrogram="both", trace = "none", col = brewer.pal(9,"PuBu"), main = "Accessory",distfun=return_same)

heatmap.2(as.matrix(dist_other_acc), dendrogram="both", trace = "none", col = brewer.pal(9,"PuBu"), main = "Roary Accessory",distfun=return_same,cexRow = 0.5, cexCol = 0.5)

#plot(as.dendrogram(hclust(p_cg_dist,method = "complete")))
#plot(as.dendrogram(hclust(p_acc_dist,method = "complete")))


dif_dist <- as.matrix(p_cg_dist) - dist_acc
#abs_dif_dist <- abs(dif_dist)

#plot(as.dendrogram(hclust(abs_dif_dist,method = "complete")))

heatmap.2(as.matrix(dif_dist), dendrogram="both", trace = "none", col = brewer.pal(9,"RdBu"), 
          main = "Core - Accessory",distfun=return_same, cexRow = 0.5, cexCol = 0.5)
#heatmap.2(as.matrix(abs_dif_dist), dendrogram="both", trace = "none", col = brewer.pal(9,"RdBu"), main = "Core - Accessory",distfun=return_same)







#Other stuff
#plot(cg_dendro, type = "rectangle")
plot(as.phylo(hc_cg), main="Core genes", cex=0.5, label.offset=10)
plot(as.phylo(hc_acc), main="Accessory genes", cex=0.5, label.offset=10)


return_same <- function(x)
{
  as.dist(x)
}


linMap <- function(x, from, to) {
  # Shifting the vector so that min(x) == 0
  x <- x - min(x)
  # Scaling to the range of [0, 1]
  x <- x / max(x)
  # Scaling to the needed amplitude
  x <- x * (to - from)
  # Shifting to the needed level
  x + from
}




