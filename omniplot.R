library(parallel)
library(magrittr)
library(pryr)
library(ggplot2)
library(clv)
library(reshape2)
library(ggthemes)

###### Data setup ##############################################################
dm <-
    'distance_matrices/pristine_distance_matrix.csv' %>%
    read.csv(row.names = 1)

dm[] <- lapply(dm, as.numeric)
dm <- as.matrix(dm)

clusters <-
    'goeburst_clusters.tsv' %>%
    read.table(row.names = 1, header = TRUE, sep = '\t')

mlst_clusters <-
    'mlst_results.txt' %>%
    read.table(row.names = 1, header = TRUE,
               sep = '\t', stringsAsFactors = FALSE) %>%
    extract(c('ST', 'CC')) %>%
    na.omit

mlst_clusters[] <- lapply(mlst_clusters, as.integer)
mlst_clusters <- subset(mlst_clusters,
                        rownames(mlst_clusters) %in% rownames(dm))


cluster_distances <-
    mclapply(clusters, cls.scatt.diss.mx, diss.mx = dm, mc.cores = 24)


threshold_stats <-
    '/home/nadege/Desktop/acc_clustering/stats_table.csv' %>%
    read.csv(row.names = 1)

##### Shared functions #########################################################

n.unique <- length %.% unique

int.height <- as.integer %.% partial(gsub, pattern = 'h_', replacement = '')

list.to.df <- function(l, cols) {
    l %>% unlist %>% matrix(byrow = TRUE, ncol = cols) %>% data.frame
}

subclusts <- function(cur, h, clusts) {
    # returns the clusts in h-1 that fall into cluster cur
    left <- colnames(clusts)[which(colnames(clusts) == h) - 1]

    sub <- subset(clusts,
                  clusts[[h]] == cur,
                  select = left,
                  drop = TRUE)
    sub

}

##### Proportion of pairwise strain distances below threshold ##################

dmlt <- dm %>% extract(lower.tri(.))

prop_lower_than_threshold <-
    threshold_stats$Threshold %>%
    sapply(function(x) {sum(dmlt < x) / length(dmlt)})

data.frame(x = threshold_stats$Threshold, y = prop_lower_than_threshold) %>%
    ggplot(aes(x, y)) +
    geom_line() +
    labs(x = 'Threshold',
         y = 'Proportion Pairwise Strain\nDistances Below Threshold')

##### Minimum intercluster distance including gaps #############################

# TODO: Fix redundant minima

min_inters <-
    cluster_distances %>%
    names %>%
    sapply(function(h) {

        mat <- cluster_distances[[c(h, 'intercls.single')]]
        height <- int.height(h)

        if (length(mat) == 1) {
            list('height' = height, 'minimum' = NA, 'count' = 1)
        } else {
            lt <- mat %>% extract(lower.tri(.))
            minimum <- min(lt)
            count <- sum(lt == minimum)

            list('height' = height, 'minimum' = minimum, 'count' = count)
        }
    })

min_inter_df <-
    min_inters %>%
    list.to.df(cols = 3) %>%
    na.omit

min_inter_df %>%
ggplot(aes(X1, X2)) +
    geom_point(aes(size = X3), alpha = 0.5) +
    labs(x = 'Threshold',
         y = 'Minimum Intercluster\nSingle-Linkage Distance') +
    scale_size(name = 'Count\n(Redundant)')

##### Strain cluster membership change count ###################################

count_cluster_changes <- function(h, clusts) {
    # Returns a list
    #
    # Keys: the clusters at height h which were formed
    # from collapsing clusters at h-1
    #
    # Values: the h-1 clusters that merged


    count_changes <- function(cur, h, clusts) {
        # For a particular cluster `cur` at height `h`,
        # how many strains experienced a cluster membership change

        subclusters <- subclusts(cur, h, clusts)

        n_subclusts <- n.unique(subclusters)

        if (n_subclusts > 1) {
            sum(clusts[, h] == cur)
        } else {
            0
        }
    }

    if (h %in% colnames(clusts)) {
        cs <- unique(clusts[, h])

        changes <- sapply(cs, count_changes,
                          h = h,
                          clusts = clusts)

        out <- list('height' = int.height(h), 'changes' = sum(changes))
    } else {
        out <- list('height' = int.height(h), 'changes' = 0)
    }

    out
}

changes <-
    paste0('h_', 1:max(threshold_stats$Threshold)) %>%
    mclapply(count_cluster_changes, clusts = clusters, mc.cores = 24)

changes.df <-
    changes %>%
    list.to.df(cols = 2) %>%
    set_colnames(c('height', 'count'))

changes.df %>%
    subset(height %in% 110:120) %>%
    ggplot(aes(height, count)) +
    geom_line() +
    labs(x = 'Threshold',
         y = 'Number of Strains Experiencing\nAltered Cluster Membership') #+
    # geom_vline(xintercept = c(115, 225), linetype = 'dotted',
               # alpha = 0.5, colour = 'red')

##### Distribution of inter:intra ##############################################
intra.inter.ratio <- function(h, clusts) {

    extract_multiples <- function(cur, h, clusts) {

        subclusters <- subclusts(cur, h, clusts)

        n_subclusts <- n.unique(subclusters)

        if (n_subclusts > 1) {
            unique(subclusters)
        } else {
            NULL
        }
    }

    calc.ratio <- function(cur, h, merges, clusts) {

        left <- colnames(clusts)[which(colnames(clusts) == h) - 1]
        left.inters <- cluster_distances[[c(left, 'intercls.single')]]

        merge.intra <-
            cluster_distances %>%
            extract2(c(h, 'intracls.complete')) %>%
            extract(as.integer(cur))

        subclusters <- merges[[cur]]

        sub.inters <-
            left.inters %>%
            extract(subclusters, subclusters) %>%
            extract(lower.tri(.))

        merge.intra / sub.inters
    }

    cs <- clusts %>% extract(, h) %>% unique

    l <-
        cs %>%
        lapply(extract_multiples, h = h, clusts = clusts) %>%
        set_names(cs)


    # names are merged cluster, values are subclusters
    merges <- Filter(`!` %.% is.null, l)

    iir <-
        merges %>%
        names %>%
        lapply(calc.ratio, h = h, merges = merges, clusts = clusts) %>%
        set_names(names(merges))

    iir
}

iir.df <-
    cluster_distances %>%
    names %>%
    mclapply(intra.inter.ratio, clusts = clusters, mc.cores = 24) %>%
    set_names(names(cluster_distances)) %>%
    Filter(f = `!` %.% (partial(equals, 0) %.% length)) %>%
    melt

iir.df$L1 <- int.height(iir.df$L1)

iir.df %>%
    # subset(L1 == 74) %>%
    ggplot(aes(L1, value, group = L1)) +
    geom_point() +
    labs(x = 'Threshold',
         y = paste('Ratios of Merged Intracluster Complete-linkage Distance',
                   'to Subcluster Intercluster Single-linkage Distance',
                   sep = '\n'))

##### Inter:Intra Scatterplot ##################################################

hs <- c(31, 74, 221, 226, 284, 287)

intra.inter.ratio2 <- function(h, clusts) {

    extract_multiples <- function(cur, h, clusts) {

        subclusters <- subclusts(cur, h, clusts)

        n_subclusts <- n.unique(subclusters)

        if (n_subclusts > 1) {
            unique(subclusters)
        } else {
            NULL
        }
    }

    calc.ratio <- function(cur, h, merges, clusts) {

        left <- colnames(clusts)[which(colnames(clusts) == h) - 1]
        left.inters <- cluster_distances[[c(left, 'intercls.single')]]

        merge.intra <-
            cluster_distances %>%
            extract2(c(h, 'intracls.complete')) %>%
            extract(as.integer(cur))

        subclusters <- merges[[cur]]

        sub.inters <-
            left.inters %>%
            extract(subclusters, subclusters) %>%
            extract(lower.tri(.))

        data.frame('height' = int.height(h),
                   'intra'  = merge.intra,
                   'inters' = mean(sub.inters))

    }

    cs <- clusts %>% extract(, h) %>% unique

    l <-
        cs %>%
        lapply(extract_multiples, h = h, clusts = clusts) %>%
        set_names(cs)


    # names are merged cluster, values are subclusters
    merges <- Filter(`!` %.% is.null, l)

    iir <-
        merges %>%
        names %>%
        lapply(calc.ratio, h = h, merges = merges, clusts = clusts) %>%
        set_names(names(merges))

    iir
}

iir.df2 <-
    cluster_distances %>%
    names %>%
    mclapply(intra.inter.ratio2, clusts = clusters, mc.cores = 24) %>%
    Filter(f = `!` %.% (partial(equals, 0) %.% length)) %>%
    lapply(rbind) %>%
    unlist %>%
    matrix(byrow = TRUE, ncol = 3) %>%
    as.data.frame %>%
    set_colnames(c('height', 'intra', 'inter'))

iir.df2 %>%
    subset(height %in% hs) %>%
    ggplot(aes(intra, inter)) +
    # facet_grid(height ~ .) +
    geom_jitter(aes(fill = factor(height)), shape = 21,
               colour = 'black', stroke = 0.5, size = 5,
               alpha = 0.7, height = 1.5, width = 1.5) +
    scale_fill_brewer('Threshold', type = 'qual',
                        palette = 'Paired', direction = -1) +
    labs(x = paste('Complete-linkage Intracluster Distance',
                   'of Composite Cluster', sep = '\n'),
         y = paste('Single-linkage Intercluster Distance',
                   'of Merging Clusters', sep = '\n'))
##### Lower quartile of intercluster distances #################################

ps <- seq(0, 1, 0.005)

prob_name <- ps %>% extract(2) %>% multiply_by(100) %>%
    as.character %>% paste0('%')

bottom_quartile_distances <-
    cluster_distances %>%
    names %>%
    sapply(function(h) {
        mat <- cluster_distances[[c(h, 'intercls.single')]]
        lt <- mat %>% extract(lower.tri(.))

        bottom_quart <- lt[lt <= quantile(lt, probs = ps)[2]]

        if (length(bottom_quart) == 0) {
            bottom_quart <- NA
        } else {
            bottom_quart <- bottom_quart - int.height(h)
        }

        bottom_quart
    }, simplify = FALSE)

bottom_quartile_melt <- melt(bottom_quartile_distances) %>% na.omit
bottom_quartile_melt$L1 <- int.height(bottom_quartile_melt$L1)

bottom_quartile_melt %>%
    subset(L1 %in% c(45, 74, 221, 226, 284, 287) ) %>%
    ggplot(aes(factor(L1), value)) +

    geom_jitter(aes(fill = factor(L1)), shape = 21, colour = 'black',
                stroke = 0.5, width = 0.4, height = 0.1) +
    geom_violin(scale = 'area', adjust = 0.1, aes(fill = factor(L1)), alpha = 0.5) +
    labs(x = 'Threshold',
         y = paste('Lowest', prob_name,
                   'of Adjusted Intercluster Single-linkage Distances')) +
    # scale_y_continuous(breaks = seq(from = 0, to = 25, by = 2)) +
    scale_fill_brewer('Threshold', type = 'qual',
                      palette = 'Paired', direction = -1) +
    geom_hline(yintercept = 0)#+
    # scale_colour_brewer('', type = 'qual', palette = 'Paired', direction = -1)

##### Complete:average intra ###################################################

c.a.ratio <- sapply(cluster_distances, function(h) {
    (h$intracls.complete / h$intracls.average)
})

c.a.melt <- na.omit(melt(c.a.ratio))
c.a.melt$L1 <- int.height(c.a.melt$L1)

ggplot(c.a.melt, aes(L1, value, group = L1)) +
    geom_boxplot(outlier.size = 0.2) +
    labs(x = 'Threshold',
         y = paste('Ratio of Complete Intercluster Distance',
                   'to Average Intercluster Distance', sep = '\n'))

##### Stability Vector Magnitude ###############################################

# Uses changes.df defined above

hypotenuse <- function(a, b) sqrt(a ^ 2 + b ^ 2)

calculate_vectors <- function(h1, dat) {


    b1 <- dat[h1, 'count']

    heights <- (h1 + 1):max(dat$height)

    max_h <-
        dat %>%
        extract(heights, 'count') %>%
        is_weakly_greater_than(b1) %>%
        which %>%
        extract(1) %>%
        add(h1)
    max_h

    cs <-
        heights %>%
        sapply(function(h2) {


            b2 <- dat[h2, 'count']

            if (h2 >= max_h) {
                NA
            } else {

                a <- (h2 - h1)
                b <- (b1 - b2)

                hypotenuse(a, b)
            }

        }) %>%
        set_names(as.character(heights)) %>%
        na.omit
        attr(cs, 'na.action') <- NULL

    shadows <-

        cs %>%
        names %>%
        extract(2:n.unique(cs)) %>%
        sapply(function(h2) {

        # trig to get the maximum permitted height of h2
            preceeding_names <- names(cs)[2:which(names(cs) == h2)]

            sapply(preceeding_names, function(h3) {

                # triangle formed by the far and intermediate peaks
                a. <- as.integer(h3) - as.integer(h2)
                b. <- dat[h3, 'count'] - dat[h2, 'count']
                print(paste(h2, h3, b.))
                c. <- hypotenuse(a., b.)

                if (b. < 0 & h2 != h3) {


                    c. < (cs[h3] - cs[h2])
                } else {
                    TRUE
                }


                # print(paste(h3, cs[h3]))


       }) #%>% all

    })

    cs
    shadows
}

clusts_at_h <- sapply(names(cluster_distances), function(h) {
    ds <-
        cluster_distances[[h]]$intercls.single %>%
        extract(lower.tri(.))

    sum(ds == int.height(h))
})
n_clusters <- sapply(names(cluster_distances), function(h) {
    nrow(cluster_distances[[h]]$intercls.single)
})

df <- data.frame(x = int.height(names(clusts_at_h)),
                 y = clusts_at_h / n_clusters)

ggplot(df, aes(x, y)) + geom_line() +
    labs(x = 'Threshold', y = 'Clusters at Threshold / Number of Clusters') +
    geom_vline(xintercept = c(115, 225), colour = 'red', linetype = 'dotted', alpha = 0.5)

df <- data.frame(x = hs[2:length(hs)],
                 y = sapply(2:length(hs), function(x) hs[x] - hs[x-1]))

ggplot(df, aes(x, y)) +
    geom_line() +
    labs(x = 'Threshold', y = 'Thresholds since last change') +
    geom_vline(xintercept = 225)

y <- cluster_distances$h_74$intercls.single %>% extract(lower.tri(.))
x <- cluster_distances$h_74$intracls.complete %>% factor

dat <- cluster_distances$h_74$intercls.single
ys <- sapply(1:ncol(dat), function(i) k.nearest.neighbors(i, dat, k = 1))
y <- sapply(1:ncol(dat), function(i) dat[ys[i], i])

plot(x/698, y/698, type = 'p', xlab = 'complete linkage intracluster',
     ylab = 'closest single linkage distance', xlim = c(0,1), ylim = c(0,1))
