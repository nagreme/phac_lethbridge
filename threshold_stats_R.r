library(pryr)

#this needs the adj_wallace function

stats_table <- function(clusters) {

  #as int for core, as double for accessory
    threshold <- as.double(gsub(pattern = 'h_', replacement = '', colnames(clusters)))
    entropy <- sapply(clusters, shannon)
    nawc <- c(NA, sapply(2:(ncol(clusters) - 5), neighbour_awc, calls = clusters), rep(NA, 5))

    df <- data.frame('Threshold' = threshold,
                     'Neighbour AWC' = nawc,
                     'Shannon' = entropy,
                     'Number of Clusters' = sapply(clusters, function(x) { length(unique(x)) }),
                     'Number of Singletons' = sapply(clusters, number_singletons),
                     check.names = FALSE)

    df

}

shannon <- function(x) {

    p_i <- function(i) {
        sum(x == i) / length(x)
    }

    -sum(sapply(unique(x), function(z) p_i(z) * log2(p_i(z))))  
}

mean_cluster_size <- function(x) {

    N <- unique(x)

    mean(sapply(N, function(i) sum(x == i)))

}

max_cluster_size <- function(x) {

    max(sapply(unique(x), function(i) sum(x == i)))

}

median_cluster_size <- function(x) {

    N <- unique(x)

    median(sapply(N, function(i) sum(x == i)))

}

#adjust hard-coded upper bound in if for different clusters
neighbour_awc <- function(i, calls) {

    if (i - 1 == 0 || i > 5240) {
        return(NA)
    }

    a <- calls[,i]
    b <- calls[,i-1]

    adj_wallace(a, b)$Adjusted_Wallace_A_vs_B
}

st_awc <- function(i, st) {

    adj_wallace(st, i)$Adjusted_Wallace_A_vs_B

}

number_singletons <- function(x) {

    N <- table(x)

    sum(N == 1)

}

counter <- function(x) {

    N <- unique(x)

    lapply(N, function(i) sum(x == i))
}

top5_and_remainder <- function(x) {
     counts <- sort(table(x), decreasing = TRUE)
     l <- c('first' = counts[1], 'second' = counts[2], 'third' = counts[3],
               'fourth' = counts[4], 'fifth' = counts[5],
               'balance' = sum(counts[6:length(counts)]))
     # l <- c(counts[1:5], sum(counts[6:length(counts)]))
     l
}
top_and_remainder <- function(x) {
    counts <- sort(table(x), decreasing = TRUE)
    l <- c('first' = counts[1], 'balance' = sum(counts[2:length(counts)]))
    # l <- c(counts[1:5], sum(counts[6:length(counts)]))
    l
}

sub_dists <- function(sub_mat) {
    lt <- sub_mat[lower.tri(sub_mat)]
}

st_dists <- function(dm, mlst_calls, st) {

    subset_ <- subset(mlst_calls, subset = mlst_calls$ST == st)
    strains <- rownames(subset_)

    dsub <- dm[rownames(dm) %in% strains, colnames(dm) %in% strains]

    distances <- sub_dists(dsub)

    df <- data.frame('st' = rep(st, length(distances)), 'distance' = distances)

    df
}

# df <- read.csv('~/Dropbox/cgMLST for Frontiers/data/threshold_statistics.csv',
#                check.names = FALSE,
#                row.names = 1)

df <- read.csv("/home/nadege/Desktop/acc_clustering/stats_table.csv", check.names = F, row.names = 1)

#df800 <- df[1:800,]

m <- melt(df, id.vars = 'Threshold')
ggplot(na.omit(m), aes(x = Threshold )) + geom_line(aes(y = value)) +
    facet_grid(variable ~ ., scales = 'free_y', switch = 'y') + labs(x = 'goeBURST Threshold', y = '') +
    scale_x_continuous(breaks = c(seq(from = 0, to = max(df$Threshold), by = 100), 115, 225)) +
    #geom_vline(xintercept = c(144), linetype="dashed")
#+
    geom_vline(xintercept = c(115, 225), linetype='dotted') + # ours
    geom_vline(xintercept = c(222, 539), linetype = 'dashed') # ST & CC
