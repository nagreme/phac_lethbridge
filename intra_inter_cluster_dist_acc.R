#Intra & inter cluster distances for the accessory

# Nadège Pulgar-Vidal
# March 13, 2017

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
binary_acc_file <- "/home/nadege/Desktop/acc_clustering/acc_presence_absence_big_set.Rtab"

clusts <- read.table(file=clusts_file, row.names = 1, header=T, sep="\t")

bin_acc_genes <- read.csv(file=binary_acc_file, row.names=1, header=T, sep=",") 
#colnames(bin_acc_genes) <- sub("X","", colnames(bin_acc_genes)) #Can't really pipe these because of this step
dist_bin_acc <- jaccard(as.matrix(bin_acc_genes))

#acc_genes <- read.table(file=acc_file, row.names = 1, header=T, sep=",")
#colnames(acc_genes) <- sub("X","", colnames(acc_genes))
#acc_genes[acc_genes == 0] <- NA 
#dist_acc <- dist.gene(acc_genes, method = "pairwise", pairwise.deletion = T, variance = F) %>% as.matrix()


#cluster_distances <- mclapply(clusts, cls.scatt.diss.mx, diss.mx= dist_bin_acc, mc.cores = 7)

cluster_distances <- mclapply(clusts, cls.scatt.diss.mx, diss.mx= as.matrix(dist_bin_acc), mc.cores = 6)

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


dist_acc <- dist.gene(acc_genes, method = "percent", pairwise.deletion = T, variance = F) %>% as.matrix()
dist_core <- dist.gene(core_genes, method = "percent", pairwise.deletion = T, variance = F) %>% as.matrix()


comp_dist_acc <- sqrt(dist_bin_acc^2 + dist_acc^2) %>% as.dist()
combined_dist <- sqrt(dist_bin_acc^2 + dist_acc^2 + dist_core^2) %>% as.dist()
allelic_dist <- sqrt(dist_core^2 + dist_acc^2) %>% as.dist()

comp_hc_acc <- hclust(comp_dist_acc, method="single")
combined_hc <- hclust(combined_dist, method="single")

#cut_heights <- unique(combined_hc$height) 

comp_acc_clusters <- cutree(comp_hc_acc, h = seq(0,1.129,0.001))
combined_clusters <- cutree(combined_hc, h = cut_heights) ###edit these numbers once you can look at the heights
#h = seq(0,1.5,0.001)

colnames(comp_acc_clusters) <- paste0("h_",colnames(comp_acc_clusters))
colnames(combined_clusters) <- paste0("h_",colnames(combined_clusters))


#adjust awc function
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
  scale_x_continuous(breaks = c(seq(from = 0, to = max(combined_df$Threshold), by = 0.01), 115, 225)) #+
  #geom_vline(xintercept = c(0.0968), linetype='dashed')  
  
ggsave(filename = "zoom_0_95_stats_table.png",path = "/home/dbarker/nadege/acc_clustering/", height = 8, width = 11)

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

# load(file = "/media/dbarker/cfia2/nadege/acc_clus_Rdat/.RData")

#individual histograms to look at distributiuon of different distances
dist_mat_2_dist_val <- function(dist_matrix)
{
  # dist_values <- unlist(as.list(tril(dist_matrix, k = -1)))
  dist_matrix[lower.tri(dist_matrix, diag = T)] <- NA
  dist_values <- dist_matrix[-which(is.na(dist_matrix))]
}

dist_histogram <- function(dist_values, min, max, bin_width, break_width, out_file, out_path)
{
  ggplot(na.omit(as.data.frame(dist_values)),aes(dist_values)) + geom_histogram(binwidth = bin_width)+
    scale_x_continuous(limit = c(min,max), breaks = c(seq(from = min, to = max, by = break_width))) +
    annotate("text", x = max*0.2, y = max(dist_values)*0.85, label = paste0("st.dev = ",sd(dist_values)))
  ggsave(file = out_file, path = out_path)
}


dist_bin_acc_values <- dist_mat_2_dist_val(dist_bin_acc)
dist_histogram(dist_bin_acc_values, 0, 1, 0.01, 0.1, "bin_acc_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")

# dist_bin_acc_values <- unlist(as.list(tril(dist_bin_acc, k = -1)))
# dist_bin_acc_values[dist_bin_acc_values == 0] <- NA
# ggplot(na.omit(as.data.frame(dist_bin_acc_values)),aes(dist_bin_acc_values)) + geom_histogram(binwidth = 0.01)+
#   scale_x_continuous(limit = c(0,1), breaks = c(seq(from = 0, to = 1, by = 0.1)))
# ggsave(file = "bin_acc_dist_histogram.png",path = "/home/dbarker/nadege/acc_clustering")

dist_acc_values <- dist_mat_2_dist_val(dist_acc)
dist_histogram(dist_acc_values, 0, 1, 0.01, 0.1, "acc_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")

# dist_acc_values <- unlist(as.list(tril(dist_acc, k = -1)))
# dist_acc_values[dist_acc_values == 0] <- NA
# ggplot(na.omit(as.data.frame(dist_acc_values)),aes(dist_acc_values)) + geom_histogram(binwidth = 0.01) +
#   scale_x_continuous(limit = c(0,1), breaks = c(seq(from = 0, to = 1, by = 0.1)))
# ggsave(file = "acc_dist_histogram.png",path = "/home/dbarker/nadege/acc_clustering")

dist_core_values <- dist_mat_2_dist_val(dist_core)
dist_histogram(dist_core_values, 0, 1, 0.01, 0.1, "core_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")

# dist_core_values <- unlist(as.list(tril(dist_core, k = -1)))
# dist_core_values[dist_core_values == 0] <- NA
# ggplot(na.omit(as.data.frame(dist_core_values)),aes(dist_core_values)) + geom_histogram(binwidth = 0.01) +
#   scale_x_continuous(limit = c(0,1), breaks = c(seq(from = 0, to = 1, by = 0.1)))
# ggsave(file = "core_dist_histogram.png",path = "/home/dbarker/nadege/acc_clustering")

comp_dist_values <- dist_mat_2_dist_val(comp_dist_acc)
dist_histogram(comp_dist_values, 0, 1.42, 0.01, 0.1, "comp_acc_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")

# comp_dist_values <- unlist(as.list(tril(as.matrix(comp_dist_acc), k = -1)))
# comp_dist_values[comp_dist_values == 0] <- NA
# ggplot(na.omit(as.data.frame(comp_dist_values)),aes(comp_dist_values)) + geom_histogram(binwidth = 0.01) +
#   scale_x_continuous(limit = c(0,1.42), breaks = c(seq(from = 0, to = 1.42, by = 0.1)))
# ggsave(file = "comp_acc_dist_histogram.png",path = "/home/dbarker/nadege/acc_clustering")

combined_dist_values <- dist_mat_2_dist_val(combined_dist)
dist_histogram(combined_dist_values, 0, 1.74, 0.01, 0.1, "combined_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")

# combined_dist_values <- unlist(as.list(tril(as.matrix(combined_dist), k = -1)))
# combined_dist_values[combined_dist_values == 0] <- NA
# ggplot(na.omit(as.data.frame(combined_dist_values)),aes(combined_dist_values)) + geom_histogram(binwidth = 0.01) +
#   scale_x_continuous(limit = c(0,1.74), breaks = c(seq(from = 0, to = 1.74, by = 0.1)))
# ggsave(file = "combined_dist_histogram.png",path = "/home/dbarker/nadege/acc_clustering")

allelic_dist_values <- dist_mat_2_dist_val(allelic_dist)
dist_histogram(allelic_dist_values, 0, 1.42, 0.01, 0.1, "allelic_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")

# allelic_dist_values <- unlist(as.list(tril(as.matrix(allelic_dist), k = -1)))
# allelic_dist_values[allelic_dist_values == 0] <- NA
# ggplot(na.omit(as.data.frame(allelic_dist_values)),aes(allelic_dist_values)) + geom_histogram(binwidth = 0.01) +
#   scale_x_continuous(limit = c(0,1.42), breaks = c(seq(from = 0, to = 1.42, by = 0.1)))
# ggsave(file = "allelic_dist_histogram.png",path = "/home/dbarker/nadege/acc_clustering")


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

allelic_core_acc_df <- data.frame("core" = unlist(as.list(tril(dist_core, k = -1))), "acc" = unlist(as.list(tril(dist_acc, k = -1))),
                                  stringsAsFactors = F)

core_comb_acc_df <- data.frame("core" = unlist(as.list(tril(dist_core, k = -1))), "combined_acc" = unlist(as.list(tril(comp_dist_acc, k = -1))),
                                  stringsAsFactors = F)

core_bin_acc_df <- data.frame("core" = unlist(as.list(tril(dist_core, k = -1))), "bin_acc" = unlist(as.list(tril(dist_bin_acc, k = -1))),
                               stringsAsFactors = F)

ggplot(core_bin_acc_df, aes(x = core, y = bin_acc)) + 
  geom_point(aes(x=core,y = bin_acc)) +
  geom_smooth(method="lm") +
  ggtitle("Corelation Between Allelic Core and Binary Accessory Distance")+
  xlab("Allelic Core Distance") +
  ylab("Binary Accessory Distance")
ggsave(file = "allelic_core_bin_acc_dist_scatterplot_line.png",path = "/home/dbarker/nadege/acc_clustering")


draw_dist_vs_dist_scatt(dist_core, dist_bin_acc, "Allelic Core", "Binary Accessory",
                        "allelic_core_bin_acc_dist_scatterplot_line.png", "/home/dbarker/nadege/acc_clustering/")


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

#ggdendrogram(cg_dendro)
draw_heatmap <- function(dist_matrix, filename, title, max)
{
  png(filename = filename, width=12, height=12, units="in", res = 450)
  dist_matrix[1,1] <- max
  heatmap.2(as.matrix(dist_matrix),Rowv = cg_dendro, Colv = "Rowv", distfun=return_same, dendrogram = "both", symm = T,
            trace = "none", col = brewer.pal(9,"PuBu"), main = paste0(title," Distance Matrix Clustered by Core"),
            key.title = "Distance Distribution", key.xlab = paste0(title, " Distance"))
  dist_bin_acc[1,1] <- 0
  dev.off()
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


draw_heatmap(subset_dist, "/home/dbarker/nadege/acc_clustering/select_bin_acc_dist_0.67-1.00_heatmap.png", "Allelic Accessory", 0)

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

#clusters_combined_distances[[height_index]]$int__cls.method[cluster_number]

dist_histogram(clusters_combined_distances[[45]]$intercls.single[1000,], 0, 1.6, 0.01, 0.1, 
               "h_45_cls_1000_intercls_single_combined_dist_histogram.png", "/home/dbarker/nadege/acc_clustering")





#histograms of intercluster distribution for a given cluster at a given height
dist_histogram_at_h <- function(dist_values, min, max, bin_width, break_width, out_file, out_path, h)
{
  ggplot(as.data.frame(dist_values), aes(dist_values)) + geom_histogram(binwidth = bin_width)+
    scale_x_continuous(limit = c(min,max), breaks = c(seq(from = min, to = max, by = break_width))) +
    labs(title = paste0("Intercluster Distances at Height ", h),
         subtitle = paste0("Standard Deviation = ", round(sd(dist_values), digits = 4), 
                           "     Min = ", round(min(dist_values), digits = 4),
                           "     Mean = ", round(mean(dist_values), digits = 4),
                           "     Median = ", round(median(dist_values), digits = 4),
                           "     Max = ", round(max(dist_values), digits = 4)),
         x = "Intercluster Single Linkage Distances",
         y = "Count")
  ggsave(file = out_file, path = out_path)
}

intercls_single_combined_dists <- dist_mat_2_dist_val(clusters_combined_distances[[45]]$intercls.single)

#histograms of all intercls distance at a given height plus standard dev.
dist_histogram_at_h(dist_values = intercls_single_combined_dists, 
               min = 0, max = max(intercls_single_combined_dists), 
               bin_width = 0.01, 
               break_width = 0.1,
               out_file = "h_45_intercls_single_combined_dist_histogram.png", 
               out_path = "/home/dbarker/nadege/acc_clustering",
               h = 45)





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
draw_plot <- function(plot_func, clusters, cluster_distances, height, savefile, savepath)
{
  intracls <- t(cluster_distances[[height]]$intracls.average)
  intercls <- non_zero_min(cluster_distances[[height]]$intercls.single)
  cluster_sizes <-  sapply(sort(unique(core_clusters[,height]),decreasing=F), function(x) {sum(core_clusters[,height] == x)})
  
  #make a data frame to make a scatterplot of intracls average vs. intercls single (at height 45)
  intra_vs_inter_df <- data.frame("cluster_number" = c(1:length(intracls)),
                                  "intracls_average" = intracls,
                                  "intercls_single" = as.numeric(intercls),
                                  "cluster_sizes" = cluster_sizes,
                                  stringsAsFactors = F)
  
  plot_func(intra_vs_inter_df, height, savefile, savepath)
}


#Interchangeable plot functions
#Scatterplot
intra_inter_scatter_at_h <- function(df, h, savefile, savepath)
{
  ggplot(df, aes((intracls_average), unlist(intercls_single))) +
    geom_point(aes(size = cluster_sizes), shape = 21, colour = "black", fill = "white", alpha = 1/2) +
    ggtitle(paste0("Intra vs Inter Distances (Combined Core and Accessory) at Height ", h))+
    xlab("Average Intracluster Distance") +
    ylab("Minimum Single Linkage Interclsuter Distance")+
    scale_y_continuous(limits = c(0,1.6)) +
    scale_x_continuous(limits = c(0,0.7)) +
    scale_size_area()
  ggsave(file = paste0(savefile,"_avg_intra_vs_sing_inter_combined_dist_scatterplot.png") ,path = savepath)
}

#Density plot
intra_inter_density_at_h <- function(intra_vs_inter_df, h, savefile, savepath)
{
  df <- melt(intra_vs_inter_df, id = c("cluster_number","cluster_sizes"))
  ggplot(df, aes(value)) +
    geom_density(aes(fill = variable), alpha = 0.8) +
    labs(title = paste0("Intra & Inter Distance Distribution at Height ", h),
         subtitle = "(Combined Core and Accessory Distances)",
         fill = "Cluster Distances") +
    #scale_y_continuous(limits = c(0,100)) +
    scale_x_continuous(limits = c(0,1.65)) 
    ggsave(file = paste0(savefile,"_avg_intra_sing_inter_combined_dist_density_plot.png") ,path = savepath, height = 8, width= 10)
}


#Histogram of min(inter)/intra
inter_intra_ratio_histogram_at_h <- function(intra_vs_inter_df, h, savefile, savepath)
{
  intra_vs_inter_df$ratio <- intra_vs_inter_df$intercls_single/intra_vs_inter_df$intracls_average
  ggplot(intra_vs_inter_df, aes(ratio)) +
    geom_histogram(binwidth = 1) +
    labs(title = paste0("Inter/Intra Ratio Distribution at Height ", h),
         subtitle = "(Combined Core and Accessory Distances)") +
    scale_x_continuous(limit = c(0,200), breaks = c(seq(from = 0, to = 200, by = 10)))
  ggsave(file = paste0(savefile,"_avg_intra_sing_inter_combined_dist_histogram.png") ,path = savepath)
}


#Multiple heights wrapper
#c(1,45,seq(25,430,25),430)
#c(1:430)
# heights <- c(1,45,100,250,430)
heights <- c(1,45,seq(25,430,25),430)

#heights <- c(1:430)
for (i in heights)
{
  filename = paste0("h_",i)
  draw_plot(plot_func = intra_inter_scatter_at_h,
            clusters = core_clusters, 
            cluster_distances = clusters_combined_distances, 
            height = i, 
            savefile = filename,
            savepath = "/home/dbarker/nadege/acc_clustering" )
}




