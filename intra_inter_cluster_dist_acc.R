#Intra & inter cluster distances for the accessory

# Nadège Pulgar-Vidal
# March 13, 2017 onwards

library(clv)
library(magrittr)
library(ape)
library(phangorn)
library(ggplot2)
library(RColorBrewer)
library(prabclus)
library(Matrix)


clusts_file <- "/home/nadege/Desktop/acc_clustering/goeburst_clusters.tsv"
#acc_file <- "/home/nadege/Desktop/mist/big_set/acc_json_out.csv"
acc_file <- "/home/dbarker/nadege/acc_clustering/acc_json_out.csv"
binary_acc_file <- "/home/nadege/Desktop/acc_clustering/acc_presence_absence_big_set.Rtab"

clusts <- read.table(file=clusts_file, row.names = 1, header=T, sep="\t")

bin_acc_genes <- read.csv(file=binary_acc_file, row.names=1, header=T, sep=",") 
#colnames(bin_acc_genes) <- sub("X","", colnames(bin_acc_genes)) #Can't really pipe these because of this step
dist_bin_acc <- jaccard(as.matrix(bin_acc_genes))

acc_genes <- read.table(file=acc_file, row.names = 1, header=T, sep=",")
colnames(acc_genes) <- sub("X","", colnames(acc_genes))
acc_genes[acc_genes == 0] <- NA
acc_genes[acc_genes == -1] <- NA
# dist_acc <- dist.gene(acc_genes, method = "pairwise", pairwise.deletion = T, variance = F) %>% as.matrix()

#another acc_genes that includes 0s and -1s
acc_genes_full <- read.table(file=acc_file, row.names = 1, header=T, sep=",")
colnames(acc_genes_full) <- sub("X","", colnames(acc_genes_full))

#redo bin_acc/gene presence absence with the msit data to compare distributions
mist_pres_abs_genes_p0 <- t(acc_genes) #assuming acc_genes are already filtered by wanted genomes
mist_pres_abs_genes_p0[is.na(mist_pres_abs_genes_p0)] <- 0 #partials and null are absence
mist_pres_abs_genes_p0[mist_pres_abs_genes_p0 > 0] <- 1 #alleles are presence
dist_mist_bin_acc_p0 <- jaccard(as.matrix(mist_pres_abs_genes_p0))


mist_pres_abs_genes_p1 <- read.table(file=acc_file, row.names = 1, header=T, sep=",")
colnames(mist_pres_abs_genes_p1) <- sub("X","", colnames(mist_pres_abs_genes_p1))
mist_pres_abs_genes_p1[mist_pres_abs_genes_p1 > 0] <- 1 #alleles and
mist_pres_abs_genes_p1[mist_pres_abs_genes_p1 == -1] <- 1 #partials are presence
#null are already as absence
mist_pres_abs_genes_p1 <- mist_pres_abs_genes_p1[which(rownames(mist_pres_abs_genes_p1) %in% wanted_genomes),]
mist_pres_abs_genes_p1 <- t(mist_pres_abs_genes_p1) #for jaccard
dist_mist_bin_acc_p1 <- jaccard(as.matrix(mist_pres_abs_genes_p1))

#cluster_distances <- mclapply(clusts, cls.scatt.diss.mx, diss.mx= dist_bin_acc, mc.cores = 7)

# cluster_distances <- mclapply(clusts, cls.scatt.diss.mx, diss.mx= as.matrix(dist_bin_acc), mc.cores = 6)

#"/home/dbarker/nadege/acc_clustering/acc_presence_absence_big_set.Rtab"
#"/home/dbarker/nadege/acc_clustering/core_json_out.csv"
#"/home/nadege/Desktop/mist/big_set/core_json_out.csv"
core_genes <- read.table(file="/home/nadege/Desktop/mist/big_set/core_json_out.csv", row.names = 1, header=T, sep=",")
colnames(core_genes) <- sub("X","", colnames(core_genes))



#Merging the two distances to have something that takes into account number of shared loci AND allele similarity for those 
#Assuming both distances are matrices

wanted_genomes <- rownames(clusters)

#filter new_df => if any(curr_genome_name == rownames)
acc_genes <- acc_genes[which(rownames(acc_genes) %in% wanted_genomes),]
core_genes <- core_genes[which(rownames(core_genes) %in% wanted_genomes),]
#bin acc was already filtered


dist_acc <- dist.gene(acc_genes, method = "percent", pairwise.deletion = T, variance = F) %>% as.matrix()
dist_core <- dist.gene(core_genes, method = "percent", pairwise.deletion = T, variance = F) %>% as.matrix()


comp_dist_acc <- sqrt(dist_bin_acc^2 + dist_acc^2) #%>% as.dist()
combined_dist <- sqrt(dist_bin_acc^2 + dist_acc^2 + dist_core^2) #%>% as.dist()
allelic_dist <- sqrt(dist_core^2 + dist_acc^2) #%>% as.dist()

comp_hc_acc <- hclust(comp_dist_acc, method="single")
combined_hc <- hclust(combined_dist, method="single")


#cut_heights <- unique(combined_hc$height) 

comp_acc_clusters <- cutree(comp_hc_acc, h = seq(0, max(round(comp_hc_acc$height + 0.005, digits = 2)),0.001))
combined_clusters <- cutree(combined_hc, h = cut_heights) ###edit these numbers once you can look at the heights
#h = seq(0,1.5,0.001)

colnames(comp_acc_clusters) <- paste0("h_",colnames(comp_acc_clusters))
colnames(combined_clusters) <- paste0("h_",colnames(combined_clusters))



hc_allelic <- hclust(as.dist(allelic_dist), method = "single")
allelic_clusters <- cutree(hc_allelic, h = seq(0, 1.4, 0.0025)) #seq(0, 1.4, 0.0025) generates 561 heights
colnames(allelic_clusters) <- lapply(colnames(allelic_clusters), function(x) sprintf("%1.4f", as.double(x)))
colnames(allelic_clusters) <- paste0("h_",colnames(allelic_clusters))

allelic_df <- stats_table(as.data.frame(allelic_clusters))
m <- melt(combined_df, id.vars = 'Threshold')



#adjust awc function (check dim and number of clusters, try to cut off when you still have three ish clusters)
#adjust as.integer -> as.double to prevent truncation

comp_acc_df <- stats_table(as.data.frame(comp_acc_clusters))
combined_df <- stats_table(as.data.frame(combined_clusters))


#comp_acc_df800 <- comp_acc_df[1:800,]
combined_df_thresh_zoom <- combined_df[1500:2500,]

m <- melt(combined_df, id.vars = 'Threshold')
combined_m_thresh <- melt(combined_df_thresh_zoom, id.vars = 'Threshold')

#scaled_combined_m["Threshold",] <- scaled_combined_m["Threshold",]*460

zoom_df <- combined_df[5075:5100,]
m <- melt(zoom_df, id.vars = 'Threshold')

ggplot(na.omit(m), aes(x = Threshold )) + 
  geom_line(aes(y = value)) +
  #geom_line(aes(y = value, color = Group)) + 
  facet_grid(variable ~ ., scales = 'free_y', switch = 'y') + labs(x = 'heights', y = '') +
  scale_x_continuous(breaks = c(seq(from = 0, to = 1.4, by = 0.1), 115, 225)) #+
  #geom_vline(xintercept = c(0.0968), linetype='dashed')  
  
ggsave(filename = "allelic_dist_stats_table.png",path = "/home/dbarker/nadege/acc_clustering/", height = 8, width = 11)

library(clv)
library(magrittr)
library(ape)
library(prabclus)
library(parallel)
library(pryr)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(Matrix)
library(ggdendro)
library(gplots)
library(RColorBrewer)
library(plyr)

# load(file = "/media/dbarker/cfia2/nadege/acc_clus_Rdat/.RData")

#individual histograms to look at distributiuon of different distances
dist_mat_2_dist_val <- function(dist_matrix)
{
  # dist_values <- unlist(as.list(tril(dist_matrix, k = -1)))
  dist_matrix[lower.tri(dist_matrix, diag = T)] <- NA
  dist_values <- dist_matrix[-which(is.na(dist_matrix))]
}


dist_bin_acc_values <- dist_mat_2_dist_val(dist_bin_acc)
dist_histogram(dist_bin_acc_values, 0, 1, 0.01, 0.1, "bin_acc_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")



dist_acc_values <- dist_mat_2_dist_val(dist_acc)
dist_histogram(dist_acc_values, 0, 1, 0.01, 0.1, "acc_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")



dist_core_values <- dist_mat_2_dist_val(dist_core)
dist_histogram(dist_core_values, 0, 1, 0.01, 0.1, "core_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")



comp_dist_values <- dist_mat_2_dist_val(comp_dist_acc)
dist_histogram(comp_dist_values, 0, 1.42, 0.01, 0.1, "comp_acc_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")



combined_dist_values <- dist_mat_2_dist_val(combined_dist)
dist_histogram(combined_dist_values, 0, 1.74, 0.01, 0.1, "combined_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")



allelic_dist_values <- dist_mat_2_dist_val(allelic_dist)
dist_histogram(allelic_dist_values, 0, 1.42, 0.01, 0.1, "allelic_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")


mist_bin_dist_values_p0 <- dist_mat_2_dist_val(dist_mist_bin_acc_p0)
dist_histogram(mist_bin_dist_values_p0, 0, 1, 0.01, 0.1, "mist_bin_acc_dist_histogram_partials_as_abs.png", "/home/dbarker/nadege/acc_clustering")

mist_bin_dist_values_p1 <- dist_mat_2_dist_val(dist_mist_bin_acc_p1)
dist_histogram(mist_bin_dist_values_p1, 0, 1, 0.01, 0.1, "mist_bin_acc_dist_histogram_partials_as_pres.png", "/home/dbarker/nadege/acc_clustering")



#Combined histograms (faceted)
dist_bin_acc_df <- data.frame("distance" = dist_bin_acc_values, "label" = "binary_acc", stringsAsFactors = F)
dist_acc_df <- data.frame("distance" = dist_acc_values, "label" = "allelic_acc", stringsAsFactors = F)
dist_core_df <- data.frame("distance" = dist_core_values, "label" = "allelic_core", stringsAsFactors = F)
dist_comb_acc_df <- data.frame("distance" = comp_dist_values, "label" = "combined_acc", stringsAsFactors = F)
dist_combined_df <- data.frame("distance" = combined_dist_values, "label" = "combined_all", stringsAsFactors = F)
dist_allelic_df <- data.frame("distance" = allelic_dist_values, "label" = "allelic_core_acc", stringsAsFactors = F)

all_dists_df <- rbind(dist_bin_acc_df,dist_acc_df,dist_core_df,dist_comb_acc_df,dist_combined_df,dist_allelic_df)

ggplot(na.omit(all_dists_df),aes(distance)) + geom_histogram(binwidth = 0.01) +
  scale_x_continuous(limit = c(0,1.74), breaks = c(seq(from = 0, to = 1.75, by = 0.25))) +
  facet_wrap(~ label)
ggsave(file = "all_dist_histograms_0.25.png",path = "/home/dbarker/nadege/acc_clustering", height = 6, width = 10)


#Scatterplot of allelic_core vs. allelic_acc
# core_sim_values <- 1 - unlist(as.list(tril(dist_core, k = -1)))
# acc_sim_values <- 1 - unlist(as.list(tril(dist_acc, k = -1)))

# allelic_core_acc_df <- data.frame("core" = unlist(as.list(tril(dist_core, k = -1))), "acc" = unlist(as.list(tril(dist_acc, k = -1))),
#                                   stringsAsFactors = F)
# 
# core_comb_acc_df <- data.frame("core" = unlist(as.list(tril(dist_core, k = -1))), "combined_acc" = unlist(as.list(tril(comp_dist_acc, k = -1))),
#                                   stringsAsFactors = F)
# 
# core_bin_acc_df <- data.frame("core" = unlist(as.list(tril(dist_core, k = -1))), "bin_acc" = unlist(as.list(tril(dist_bin_acc, k = -1))),
#                                stringsAsFactors = F)
# 
# ggplot(core_bin_acc_df, aes(x = core, y = bin_acc)) + 
#   geom_point(aes(x=core,y = bin_acc)) +
#   geom_smooth(method="lm") +
#   ggtitle("Corelation Between Allelic Core and Binary Accessory Distance")+
#   xlab("Allelic Core Distance") +
#   ylab("Binary Accessory Distance")
# ggsave(file = "allelic_core_bin_acc_dist_scatterplot_line.png",path = "/home/dbarker/nadege/acc_clustering")


draw_dist_vs_dist_scatt(dist_bin_acc, dist_acc, "Binary Accessory", "Allelic Accessory",
                        "bin_acc_vs_acc_dist_scatterplot_line.png", "/home/dbarker/nadege/acc_clustering/", percent = F)


draw_dist_vs_dist_scatt <- function(dist1, dist2, title1, title2, outfile, outpath, percent = T)
{
  list1 <- dist_mat_2_dist_val(dist1)
  list2 <- dist_mat_2_dist_val(dist2)
  
  if (percent)
  {
    list1 <- list1/max(list1)
    list2 <- list2/max(list2)
  }
  
  dist_df <- data.frame(label1 = list1, label2 = list2, stringsAsFactors = F)
  
  g <- ggplot(dist_df, aes(x = label1, y = label2)) + 
    geom_point() +
    geom_smooth(method="lm") +
    ggtitle(paste0("Corelation Between ", title1 ," and ", title2, " Distance")) +
    xlab(paste0(title1," Distance")) +
    ylab(paste0(title2, " Distance"))
  
   if (percent)
   {
     g <- g + xlim(0,1) + ylim(0,1)
   }

  ggsave(plot = g, filename = outfile, path = outpath)
}


#Heatmaps (all dists onto core clustering)
return_same <- function(x)
{
  as.dist(x)
}

cg_dendro <- hclust(as.dist(dist_core), method = "single") %>% as.dendrogram()

allelic_dendro <-hclust(as.dist(allelic_dist), method = "single") %>% as.dendrogram()


# Coloured heatmaps and histograms

#combined_dist
# draw_heatmap(combined_dist, 
#              "/home/dbarker/nadege/acc_clustering/coloured_combined_dist_heatmap.png",
#              "Combined", 0, colours = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#B3DE69", "#80B1D3", "#FDB462", "#FB8072"), 
#              col_breaks = c(0, 0.35, 0.8, 0.95, 1.2, 1.4, 1.55, 1.7))
#dist_combined_df
# colour_dist_histogram(dist_values_df = dist_combined_df, min = 0, max = max(dist_combined_df$distance),
#                       bin_width = 0.01, break_width = 0.1,
#                       colours = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#B3DE69", "#80B1D3", "#FDB462", "#FB8072"),
#                       col_breaks = c(0, 0.35, 0.8, 0.95, 1.2, 1.4, 1.55, round(max(dist_combined_df$distance)+0.005, digits = 2)),
#                       out_file = "combined_dist_coloured_histogram.png",
#                       out_path = "/home/dbarker/nadege/acc_clustering/")

#dist_core
# draw_heatmap(dist_core,
#              "/home/dbarker/nadege/acc_clustering/coloured_core_dist_heatmap.png",
#              "Allelic Core", 0, colours = brewer.pal(6,"Set3"),
#              col_breaks = c(0, 0.2, 0.45, 0.57, 0.77, 0.92, 1),
#              dendrogram = allelic_dendro,
#              clstrd_by = "Allelic Core and Accessory")
#dist_core_df
# colour_dist_histogram(dist_values_df = dist_core_df, min = 0, max = max(dist_core_df$distance),
#                       bin_width = 0.01, break_width = 0.1,
#                       colours = brewer.pal(6,"Set3"),
#                       col_breaks = c(0, 0.2, 0.45, 0.57, 0.77, 0.92, round(max(dist_core_df$distance)+0.005, digits = 2)),
#                       out_file = "core_dist_coloured_histogram.png",
#                       out_path = "/home/dbarker/nadege/acc_clustering/")

#dist_acc
# draw_heatmap(dist_acc,
#              "/home/dbarker/nadege/acc_clustering/coloured_acc_dist_heatmap.png",
#              "Allelic Accessory", 0, colours = brewer.pal(6,"Set3"),
#              col_breaks = c(0, 0.25, 0.53, 0.65, 0.8, 0.93, 1),
             # dendrogram = allelic_dendro,
             # clstrd_by = "Allelic Core and Accessory")
#dist_acc_df
# colour_dist_histogram(dist_values_df = dist_acc_df, min = 0, max = max(dist_acc_df$distance),
#                       bin_width = 0.01, break_width = 0.1,
#                       colours = brewer.pal(6,"Set3"),
#                       col_breaks = c(0, 0.25, 0.53, 0.65, 0.8, 0.93, round(max(dist_acc_df$distance)+0.005, digits = 2)),
#                       out_file = "acc_dist_coloured_histogram.png",
#                       out_path = "/home/dbarker/nadege/acc_clustering/")

#dist_bin_acc
#brewer.pal(5,"Set3")
#c(0, 0.2, 0.35, 0.5, 0.67, 0.8)
# draw_heatmap(dist_bin_acc,
#              "/home/dbarker/nadege/acc_clustering/coloured_bin_acc_dist_heatmap.png",
#              "Binary Accessory", 0, colours = brewer.pal(5,"Set3"),
#              col_breaks = c(0, 0.2, 0.35, 0.5, 0.67, 0.8),
#              dendrogram = allelic_dendro,
#              clstrd_by = "Allelic Core and Accessory")
#dist_bin_acc_df
# colour_dist_histogram(dist_values_df = dist_bin_acc_df, min = 0, max = max(dist_bin_acc_df$distance),
#                       bin_width = 0.01, break_width = 0.1,
#                       colours = brewer.pal(5,"Set3"),
#                       col_breaks = c(0, 0.2, 0.35, 0.5, 0.67, round(max(dist_bin_acc_df$distance)+0.005, digits = 2)),
#                       out_file = "bin_acc_dist_coloured_histogram.png",
#                       out_path = "/home/dbarker/nadege/acc_clustering/")

#comp_dist_acc
# draw_heatmap(comp_dist_acc,
#              "/home/dbarker/nadege/acc_clustering/coloured_combined_acc_dist_heatmap.png",
#              "Combined Allelic", 0, colours = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#B3DE69", "#80B1D3", "#FDB462", "#FB8072"),
#              col_breaks = c(0, 0.3, 0.61, 0.75, 0.9, 1.06, 1.2, 1.3))
# colour_dist_histogram(dist_values_df = dist_comb_acc_df, min = 0, max = max(dist_comb_acc_df$distance),
#                       bin_width = 0.01, break_width = 0.1,
#                       colours = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#B3DE69", "#80B1D3", "#FDB462", "#FB8072"),
#                       col_breaks = c(0, 0.3, 0.61, 0.75, 0.9, 1.06, 1.2, round(max(dist_comb_acc_df$distance)+0.005, digits = 2)),
#                       out_file = "comp_acc_dist_coloured_histogram.png",
#                       out_path = "/home/dbarker/nadege/acc_clustering/")


#allelic_dist
# draw_heatmap(allelic_dist,
#              "/home/dbarker/nadege/acc_clustering/coloured_allelic_dist_heatmap.png",
#              "Allelic", 0, colours = brewer.pal(6,"Set3"),
#              col_breaks = c(0, 0.3, 0.74, 0.87, 1.14, 1.34, 1.42),
#              dendrogram = allelic_dendro,
#              clstrd_by = "Allelic Core and Accessory")
#dist_allelic_df
# colour_dist_histogram(dist_values_df = dist_allelic_df, min = 0, max = max(dist_allelic_df$distance),
#                       bin_width = 0.01, break_width = 0.1,
#                       colours = brewer.pal(6,"Set3"),
#                       col_breaks = c(0, 0.3, 0.74, 0.87, 1.14, 1.34, round(max(dist_allelic_df$distance)+0.005, digits = 2)),
#                       out_file = "allelic_dist_coloured_histogram.png",
#                       out_path = "/home/dbarker/nadege/acc_clustering/")


# allelic_dendro
#ggdendrogram(cg_dendro)
draw_heatmap <- function(dist_matrix, filename, title, max, colours, col_breaks, dendrogram, clstrd_by)
{
  png(filename = filename, width=12, height=12, units="in", res = 450)
  dist_matrix[1,1] <- max
  heatmap.2(as.matrix(dist_matrix),Rowv = dendrogram, Colv = "Rowv", distfun=return_same, dendrogram = "both", symm = T,
            trace = "none", col = colours, breaks = col_breaks, main = paste0(title," Distance Matrix Clustered by ",clstrd_by),
            key.title = "Distance Distribution", key.xlab = paste0(title, " Distance"))
  dist_bin_acc[1,1] <- 0
  dev.off()
}

colour_dist_histogram(dist_values_df = dist_core_df, min = 0, max = max(dist_core_df$distance),
                      bin_width = 0.01, break_width = 0.1,
                      colours = brewer.pal(6,"Set3"),
                      col_breaks = c(0, 0.2, 0.45, 0.57, 0.77, 0.92, round(max(dist_core_df$distance)+0.005, digits = 2)),
                      out_file = "core_dist_coloured_histogram.png",
                      out_path = "/home/dbarker/nadege/acc_clustering/")

#coloured_histograms
#requires a dist_value_df of this format:
# distance        label
# 1 1.3679782 combined_all
# 2 1.3793407 combined_all
# 3 1.0883384 combined_all
colour_dist_histogram <- function(dist_values_df, min, max, bin_width, break_width, 
                                  colours, col_breaks, out_file, out_path)
{
  # dist_values_df$group <- cut(dist_values_df$distance, breaks = col_breaks, right = FALSE)

  g <- ggplot(dist_values_df, aes(distance)) + 
    geom_histogram(binwidth = bin_width, fill = "white", colour = "black") + 
    scale_x_continuous(limit = c(min,max), breaks = c(seq(from = min, to = max, by = break_width))) 
  
  for (i in c(1:length(col_breaks)-1))
  {
    g <- g + geom_histogram(data = subset(dist_values_df, distance > col_breaks[i] & distance < col_breaks[i+1]), 
                            binwidth = bin_width, fill = colours[i])
  }

  ggsave(plot = g, file = out_file, path = out_path)
}


#subset heatmaps
# subset_dist <- combined_dist

subset_dist <- ifelse(dist_bin_acc < 0.67 | dist_bin_acc > 1, 0, dist_bin_acc)

# subset_dist[subset_dist < 0.2 | subset_dist > 0.3] <- 0
# subset_dist[subset_dist < 0.8 | subset_dist > 0.95] <- 0
# subset_dist[subset_dist < 0.95 | subset_dist > 1.2] <- 0
# subset_dist[subset_dist < 1.25 | subset_dist > 1.4] <- 0
# subset_dist[subset_dist < 1.4 | subset_dist > 1.55] <- 0
# subset_dist[subset_dist < 1.55 | subset_dist > 1.6] <- 0


draw_heatmap(subset_dist, "/home/dbarker/nadege/acc_clustering/select_bin_acc_dist_0.67-1.00_heatmap.png", 
             "Allelic Accessory", 0)

# png(filename = "/home/dbarker/nadege/acc_clustering/bin_acc_heatmap.png", width=12, height=12, units="in", res = 450)
# dist_bin_acc[1,1] <- 1.7
# heatmap.2(as.matrix(dist_bin_acc),Rowv = cg_dendro, Colv = "Rowv", distfun=return_same, dendrogram = "both", symm = T,
#           trace = "none", col = brewer.pal(9,"PuBu"), main = "Accessory Distance Matrix Clustered by Core",
#           key.title = "Accessory Distance Distribution", key.xlab = "Accessory Distance", scale ="row")
# dist_bin_acc[1,1] <- 0 
# dev.off()


#Intra/Inter cluster distances at height 45 
# h_45_clusters <- core_clusters[,"h_45"] %>% matrix(ncol = 1, nrow = 5257)
# colnames(h_45_clusters)  <- colnames(core_clusters)
#h_45_cluster_distances <- cls.scatt.diss.mx(diss.mx= as.matrix(___), h_45_clusters)

clusters_combined_distances <- mclapply(core_clusters, cls.scatt.diss.mx, diss.mx= as.matrix(combined_dist), mc.cores = 12)
clusters_bin_acc_distances <- mclapply(core_clusters, cls.scatt.diss.mx, diss.mx= as.matrix(dist_bin_acc), mc.cores = 12)

#vvv returns errors on 4+ threads (2 and 3 were not tested)
clusters_allelic_distances <- mclapply(as.data.frame(allelic_clusters), cls.scatt.diss.mx, 
                                       diss.mx= allelic_dist, mc.cores = 1)




#clusters_combined_distances[[height_index]]$int__cls.method[cluster_number]

dist_histogram(clusters_combined_distances[[45]]$intercls.single[1000,], 0, 1.6, 0.01, 0.1, 
               "h_45_cls_1000_intercls_single_combined_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")




#histograms of intercluster distribution at given height
dist_histogram_at_h <- function(dist_values, min, max, bin_width, break_width, out_file, out_path, h)
{
  ggplot(as.data.frame(dist_values), aes(dist_values)) + geom_histogram(binwidth = bin_width)+
    scale_x_continuous(limit = c(min,max), breaks = c(seq(from = min, to = max, by = break_width))) +
    labs(title = paste0("Intercluster Distances at ", h),
         subtitle = paste0("Standard Deviation = ", round(sd(dist_values), digits = 4), 
                           "     Min = ", round(min(dist_values), digits = 4),
                           "     Mean = ", round(mean(dist_values), digits = 4),
                           "     Median = ", round(median(dist_values), digits = 4),
                           "     Max = ", round(max(dist_values), digits = 4)),
         x = "Intercluster Single Linkage Distances",
         y = "Count")
  ggsave(file = out_file, path = out_path)
}

#All intercluster distances at given height (all heights)
for (i in c(1:length(clusters_combined_distances)))
{
  intercls_single_combined_dists <- dist_mat_2_dist_val(clusters_combined_distances[[i]]$intercls.single)
  
  height <- colnames(core_clusters)[i]

  #histograms of all intercls distance at a given height plus standard dev.
  dist_histogram_at_h(dist_values = intercls_single_combined_dists, 
                 min = 0, max = max(intercls_single_combined_dists), 
                 bin_width = 0.01, 
                 break_width = 0.1,
                 out_file = paste0(height,"_intercls_single_combined_dist_histogram.png"), 
                 out_path = "/home/dbarker/nadege/acc_clustering",
                 h = height)
}


#All intercluster distances at given height and given cluster
# for (h in c(1:length(clusters_combined_distances)))
for (h in c(1:5))
{
  height <- colnames(core_clusters)[h]
  
  clusts_at_percent <- 0.2
  
  for (cls in seq(1, length(unique(core_clusters[,h])), length(unique(core_clusters[,h]))*clusts_at_percent))
  {
    cls <- as.integer(cls)
    intercls_single_combined_dists <- clusters_combined_distances[[h]]$intercls.single[cls,]
    
    #histograms of all intercls distance at a given height plus standard dev.
    dist_histogram_at_h(dist_values = intercls_single_combined_dists, 
                        min = 0, max = max(intercls_single_combined_dists), 
                        bin_width = 0.01, 
                        break_width = 0.1,
                        out_file = paste0(height,"_cls_",cls,"_intercls_single_combined_dist_histogram.png"), 
                        out_path = "/home/dbarker/nadege/acc_clustering",
                        h = paste0(height,"_cls_",cls))
  }
}



#Scatterplot of height & cls # vs. standard deviations of intercls distances 
for (h in c(1:length(clusters_combined_distances)))
# for (h in c(1, 45, 230, 400))
{
  height <- colnames(core_clusters)[h]
   
  df <- data.frame(x = sapply(sort(unique(core_clusters[,h])), function(x) sum(core_clusters[,h] == x)), 
                   y = apply(clusters_combined_distances[[h]]$intercls.single, 1, sd))
  draw_scatt(df, 
             xlb = "Cluster Size",
             ylb = "Std Dev. of Intercluster Distances",
             xmax = 5300,
             ymax = 1.2,
             heading = paste0("Correlation Btw Cluster Size and Intercls Distance Std. Dev at ", height),
             out_file = paste0(height,"cls_size_vs_intercls_dist_scatterplot.png"),
             out_path = "/home/dbarker/nadege/acc_clustering")
}


draw_scatt <- function(df, xlb, ylb, xmax, ymax, heading, out_file, out_path)
{
  ggplot(df, aes(df$x, df$y)) +
    geom_point(shape = 21, colour = "black", fill = "white", alpha = 0.6) +
    labs(title = heading,
         x = xlb,
         y = ylb)+
    expand_limits(x = c(0, xmax), y = c(0, ymax))
  
  ggsave(file = out_file, path = out_path)
}


# 
#   ggplot(df, aes((intracls_average), unlist(intercls_single))) +
#   geom_point(aes(size = cluster_sizes), shape = 21, colour = "black", fill = "white", alpha = 1/2) +
#   ggtitle(paste0("Intra vs Inter Distances (Combined Core and Accessory) at Height ", h))+
#   xlab("Average Intracluster Distance") +
#   ylab("Minimum Single Linkage Interclsuter Distance")+
#   scale_y_continuous(limits = c(0,1.6)) +
#   scale_x_continuous(limits = c(0,0.7)) +
#   scale_size_area()
# ggsave(file = paste0(savefile,"_avg_intra_vs_sing_inter_combined_dist_scatterplot.png") ,path = savepath)
#   


dist_histogram <- function(dist_values, min, max, bin_width, break_width, out_file, out_path)
{
  ggplot(na.omit(as.data.frame(dist_values)),aes(dist_values)) + geom_histogram(binwidth = bin_width)+
    scale_x_continuous(limit = c(min,max), breaks = c(seq(from = min, to = max, by = break_width)))
  ggsave(file = out_file, path = out_path)
}



# draw_plot(plot_func = intra_inter_scatter_at_h,
#           clusters = core_clusters, 
#           cluster_distances = clusters_combined_distances, 
#           height = i, 
#           savefile = filename,
#           savepath = "/home/dbarker/nadege/acc_clustering" )

#min(intercls[intercls[,cluster] > 0, cluster])
#make a list of the smallest non-zero values by column
non_zero_col_min <- function(matrix)
{
  tmp <- list()
  for (i in 1:ncol(matrix))
  {
    tmp[[i]] <- min(matrix[matrix[,i] > 0, i])
  }
  as.matrix(tmp)
}


#Scatterplots of intra vs min(inter) at given height
#Density plots of intra and inter cluster distances at a given height

#Plot function wrapper
draw_plot <- function(plot_func, clusters, cluster_distances, height, savepath, dist_type)
{
  intracls <- t(cluster_distances[[height]]$intracls.average)
  intercls <- non_zero_col_min(cluster_distances[[height]]$intercls.single)

  #make a data frame to make a scatterplot of intracls average vs. intercls single (at height 45)
  intra_vs_inter_df <- data.frame("x" = c(1:length(intracls)),
                                  "intracls_average" = intracls,
                                  "intercls_single" = as.numeric(intercls),
                                  stringsAsFactors = F)
  
  #Add information about cluster size
  intra_vs_inter_df <- join(intra_vs_inter_df, count(clusters[,height]), by = "x")
  colnames(intra_vs_inter_df)[1] <- "cluster_number"
  
  plot_func(intra_vs_inter_df, height, savepath, dist_type)
}


#Interchangeable plot functions
#Scatterplot
intra_inter_scatter_at_h <- function(df, h, savepath, dist_type)
{
  seg_coords <- data.frame(x1 = 0, x2 = 0.5, x3 = 0.1, y1 = 0, y2 = 0.5)
  
  ggplot(df, aes((intracls_average), unlist(intercls_single))) +
    geom_point(aes(size = freq), shape = 21, colour = "black", fill = "white", alpha = 1/2) +
    labs(title = paste0("Intra vs Inter Distances at ", h),
         subtitle = paste0("(",dist_type," Distances)"),
         x = "Average Intracluster Distance",
         y = "Minimum Single Linkage Interclsuter Distance") +
    scale_y_continuous(limits = c(0,1.6)) +
    scale_x_continuous(limits = c(0,0.7)) +
    scale_size_area() +
    geom_segment(aes(x = 0, xend = 0.6, y = 0, yend = 0.6), color = "red") + 
    geom_segment(aes(x = 0, xend = 0.3, y = 0, yend = 1.5), color = "blue")
  ggsave(file = paste0(h,"_avg_intra_vs_sing_inter_combined_dist_scatterplot.png") ,path = savepath)
}

#Density plot
intra_inter_density_at_h <- function(intra_vs_inter_df, h, savepath, dist_type)
{
  df <- melt(intra_vs_inter_df, id = c("cluster_number","cluster_sizes"))
  ggplot(df, aes(value)) +
    geom_density(aes(fill = variable), alpha = 0.8) +
    labs(title = paste0("Intra & Inter Distance Distribution at ", h),
         subtitle = paste0("(",dist_type," Distances)"),
         fill = "Cluster Distances") +
    #scale_y_continuous(limits = c(0,100)) +
    scale_x_continuous(limits = c(0,1.65)) 
    ggsave(file = paste0(h,"_avg_intra_sing_inter_combined_dist_density_plot.png") ,path = savepath, height = 8, width= 10)
}


#Histogram of intra/min(inter) (minimization)
intra_inter_ratio_histogram_at_h <- function(intra_vs_inter_df, h, savepath, dist_type)
{
  intra_vs_inter_df$ratio <- intra_vs_inter_df$intracls_average/intra_vs_inter_df$intercls_single
  # intra_vs_inter_df$ratio <- intra_vs_inter_df$intercls_single/intra_vs_inter_df$intracls_average
  ggplot(intra_vs_inter_df, aes(ratio)) +
    geom_histogram(binwidth = 0.005) +
    labs(title = paste0("Intra/Inter Ratio Distribution at ", h),
         subtitle = paste0("(",dist_type," Distances)")) +
    scale_x_continuous(limit = c(0,0.5), breaks = c(seq(from = 0, to = 0.5, by = 0.05))) +
    scale_y_continuous(limit = c(0, 20), breaks = c(seq(0, 20, 2)))
    theme_bw()
  ggsave(file = paste0(h,"_sing_inter_avg_intra_allelic_dist_histogram.png") ,path = savepath)
}


#Multiple heights wrapper
#c(1,45,seq(25,430,25),430)
#c(1:430)
# heights <- c(1,45,100,250,430)
# heights <- c(1,45,seq(25,430,25),430)

#c(1:560) for allelic clusters

# heights <- c(1,50,100,150,200,250,300,350,400,450,500,550)
heights <- c(1:560)
for (i in heights)
{
  draw_plot(plot_func = intra_inter_ratio_histogram_at_h,
            clusters = allelic_clusters, 
            cluster_distances = clusters_allelic_distances,  #clusters_combined_distances, clusters_bin_acc_distances
            height = paste0(colnames(allelic_clusters)[i]), 
            savepath = "/home/dbarker/nadege/acc_clustering" ,
            dist_type = "Allelic Core and Accessory")  #"Combined Core and Accessory", "Binary Accessory"
}


#Investigating the genomes at the edge of the heatmap/dendro and what's causing the little 
#spike on the right side of the dist_bin_acc histogram
target_genome_indices <- order.dendrogram(cg_dendro)[1:14]
target_genomes <- rownames(dist_core)[target_genome_indices]

# counts <- lapply(target_genomes, function(genome) 
#   { 
#     apply(bin_acc_genes, 1, function(gene_row) 
#       { 
#         if (gene_row[genome] == 1)
#         {
#           sum(gene_row)
#         } #if
#       }) #apply (gene rows)
#   }) #lapply (target_genomes)

row_counts <- apply(bin_acc_genes, 1, sum) %>% as.list()


# genes <- row_count[which(bin_acc_genes[,target_genomes] == 1)]


genome_genes <- lapply(target_genomes, function(x) { row_counts[which(bin_acc_genes[,x] == 1)]})
# genome_genes_control <- lapply(c(3000:3005), function(x) { row_counts[which(bin_acc_genes[,x] == 1)]})

gene_count_histogram <- function(dist_values, xmin, xmax, ymin, ymax, bin_width, break_width, out_file, out_path)
{
  ggplot(na.omit(as.data.frame(dist_values)),aes(as.data.frame(dist_values))) + 
    geom_histogram(binwidth = bin_width, center = 0)+
    scale_x_continuous(limit = c(xmin,xmax), breaks = c(seq(from = xmin, to = xmax, by = break_width))) +
    scale_y_continuous(limit = c(ymin, ymax))
  ggsave(file = out_file, path = out_path)
}

for (i in c(1:14))
{
  gene_count_histogram(dist_values = t(data.frame(genome_genes[[i]])/max(as.data.frame(genome_genes))),
                 xmin = 0,
                 xmax = 1,
                 ymin = 0,
                 ymax = 100,
                 bin_width = 0.005,
                 break_width = 0.1,
                 out_file = paste0(target_genomes[i],"_genes_occurrences_histogram.png"),
                 out_path = "/home/dbarker/nadege/acc_clustering")
}





#Get the values to build the dataframe used to plot ratios (and stability)
ratios_weighted_sums_at_h <- function(height, cluster_distances, clusters)
{
  intracls <- t(cluster_distances[[height]]$intracls.average)
  intercls <- non_zero_col_min(cluster_distances[[height]]$intercls.single)
  
  df <- data.frame("x" = c(1:length(intracls)),
                                  "intracls_average" = intracls,
                                  "intercls_single" = as.numeric(intercls),
                                  stringsAsFactors = F)
  
  #Add information about cluster size
  df <- join(df, count(clusters[,height]), by = "x")
  colnames(df)[1] <- "cluster_number"
  
  #intra/inter ratio
  df$ratio <- df$intracls_average/df$intercls_single
  
  #weighted
  df$cls_ratio <- df$freq * df$ratio
  
  abs_height <- colnames(allelic_clusters)[height]
  ratio_w_sum <- sum(df$cls_ratio)
  intracls_w_sum <- sum(df$intracls_average * df$freq)
  intercls_w_sum <- sum(df$intercls_single * df$freq)
  sum_ratio <- intracls_w_sum/intercls_w_sum
  stability <- intercls_w_sum - (height-1)*0.0025*5257 #0.0025 is the step size and 5257 is the number os strains/genomes
  #the math (rearrangement) behind this is in my lab journal Apr 19, 2017
    
  data.frame(abs_height, height, ratio_w_sum, intracls_w_sum, intercls_w_sum, sum_ratio, stability)
}

ratios_lists <- lapply(c(1:557), function(x) ratios_weighted_sums_at_h(x, clusters_allelic_distances, allelic_clusters))

ratios_df <- do.call("rbind", ratios_lists)

colnames(ratios_df) <- c("abs_height", "rel_height", "ratio_w_sum", "intracls_w_sum", "intercls_w_sum", "sum_ratio", "stability")


m <- melt(ratios_df, id.vars = c("abs_height", "rel_height"), measure.vars = c("ratio_w_sum", "intracls_w_sum", "intercls_w_sum"))


#tri-line: weighted sum of ratios of intra/inter and individual weighted sums of inter and intra at height
ggplot(m, aes(x = rel_height)) +
  geom_line(aes(y = value, color = variable)) +
  labs(title = "Weighted Sums by Height",
       x = "i-th height (constant 0.0025 steps)",
       y = "Weighted sums") +
  scale_x_continuous(breaks = seq(0, 560, 25))
ggsave(file = "intra_inter_w_sums_at_h_allelic_dist.png", path = "/home/dbarker/nadege/acc_clustering", width = 9)

#ratio of weigthed sums of intra/inter at height
ggplot(ratios_df, aes(x = rel_height)) +
  geom_line(aes(y = sum_ratio)) +
  scale_x_continuous(breaks = seq(0, 560, 25)) +
  labs(title = "Ratio of Weighted Sums by Height",
       x = "i-th height (constant 0.0025 steps)",
       y = "Weighted sum")
ggsave(file = "intra_inter_w_sums_ratio_at_h_allelic_dist.png", path = "/home/dbarker/nadege/acc_clustering")



n.unique <- length %.% unique
#Allele distribution in the accessory

#histograms of distribution of number of alleles per gene
num_allele_df <- data.frame(num_alleles = apply(acc_genes, 2, n.unique))

ggplot(num_allele_df, aes(num_alleles)) +
  geom_histogram(binwidth = 5) +
  labs(title = "Distribution of Number of Alleles per Gene",
       subtitle = "(Accessory Genes)",
       x = "Alleles per gene",
       y = "Frequency") 
ggsave(filename = "num_alleles_in_acc_histogram.png", path = "/home/dbarker/nadege/acc_clustering")



#bar/col plot of allele distribution for each gene

allele_distribution_histogram <- function(genes_df, gene_col_index)
{
  alleles_df <- count(genes_df[,gene_col_index])
  
  ggplot(alleles_df, aes(x, freq)) +
    geom_col() +
    labs(title = colnames(genes_df)[gene_col_index],
         subtitle = "Distribution of alleles",
         x = "Alleles",
         y = "Frequency")
  ggsave(filename = paste0("allele_distribution_in_", colnames(genes_df)[gene_col_index], ".png"), path = "/home/dbarker/nadege/acc_clustering")
}

lapply(c(1:ncol(acc_genes_full)), function(x) allele_distribution_histogram(acc_genes_full, x))

# lapply(c(1:3), function(x) allele_distribution_histogram(acc_genes_full, x))



#Allele stability line plot (reuse ratios_df)
ggplot(ratios_df, aes(x = rel_height))+
  geom_line(aes(y = stability)) +
  labs(title = "Stability", 
       subtitle = "(weighted sum of min intercls distance single linkage) - (distance from threshold)",
       x = "i-th height (constant 0.0025 steps)",
       y = "Stability") +
  scale_x_continuous(breaks = seq(0,575, 25))
ggsave(filename = "stability_line_allelic.png", path = "/home/dbarker/nadege/acc_clustering", width = 12)

