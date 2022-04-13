library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(traj)
library(ggdendro)

# https://uc-r.github.io/hc_clustering
# https://uc-r.github.io/kmeans_clustering
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# https://cran.r-project.org/web/packages/traj/vignettes/trajVignette.pdf
# https://bradleyboehmke.github.io/HOML/hierarchical.html
# wards method: https://www.statisticshowto.com/wards-method/
# wards method: https://online.stat.psu.edu/stat505/lesson/14/14.7
# how to get centroid avg: https://sciencing.com/centroid-clustering-analysis-10070345.html
# https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#a-dendrogram-is-a-nested-list-of-lists-with-attributes
# https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
# https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html

df <- USArrests %>% drop_na() %>% mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1]))
df
df %>% glimpse()

# check scaling
scale(USArrests$Murder, center = TRUE, scale = TRUE)[ , 1]


#//////////////////


# cluster with hclust ####

# get dissimilarity matrix
d <- dist(df, method = "euclidean")
d
d %>% attributes()

# inspect distance_matrix_tbl
distance_matrix_tbl <- as.matrix(d) %>% as_tibble() %>% 
        mutate(state = row.names(as.matrix(d))) %>%
        relocate(state, .before = everything())

# visualize distance
fviz_dist(dist.obj = d, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# check euclidena distance calculation for alaska & alabama (2.7037)
df %>% rownames_to_column() %>% filter(rowname %in% c("Alaska", "Alabama"))
sqrt( (1.2425641 - 0.5078625)^2 + (0.7828393 - 1.1068225)^2 + (-0.5209066 - -1.2117642)^2 + (-0.003416473 - 2.484202941)^2 )

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d = d, method = "complete" )
hc1
hc1 %>% attributes()
hc1$merge
hc1$height

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)


#////////////////////


# cluster with agnes ####

hc2 <- agnes(df, method = "complete")

# Agglomerative coefficient 
# "Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 suggest a 
# more balanced clustering structure such as the complete linkage and Wards method dendrograms in Figure 21.3. 
# Values closer to 0 suggest less well-formed clusters such as the single linkage dendrogram in Figure 21.3. 
# However, the AC tends to become larger as n increases, so it should not be used to compare across data sets of very different sizes."
hc2$ac


#//////////////////



# use agnes to look for clustering method with best agglomerative coefficient ####
method_list <- c( "average", "single", "complete", "ward")
names(method_list) <- c( "average", "single", "complete", "ward")
method_list

# function to compute coefficient
cluster_methods <- function(current_method) {
        print(current_method)
        agnes(df, method = current_method)$ac
}

# run cluster_methods
# note that wards method has the highest/best agglomerative coefficient
map(.x = method_list, .f = ~ cluster_methods(current_method = .x))


#///////////////////


# visualize clustering ####
hc3 <- agnes(df, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 


# Ward's method
# note that for some reason agnes(method = "ward") is slightly different dendrogram than hclust(method = "ward.D2")
#  ward.D2 is an updated criteria
hc5 <- hclust(d, method = "ward.D2" )
plot(hc5, cex = 0.6, hang = -1)


#///////////////////


# Cut tree ####
sub_grp <- cutree(hc5, k = 4)
sub_grp

# Number of members in each cluster
table(sub_grp)

# can add cutree output to df
df %>% mutate(cluster = sub_grp)

# plot dendrogram with border around 4 clusters
plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 4, border = 2:5)

# visualize clusters along 2 primary PCA dimensions
fviz_cluster(list(data = df, cluster = sub_grp))


#/////////////////////


# can also use cutree on agnes

# Cut agnes() tree into 4 groups
hc_a <- agnes(df, method = "ward")
cutree(as.hclust(hc_a), k = 4)


#//////////////////////////


# can use tanglegram to compared two dendrograms ####

# Compute distance matrix
res.dist <- dist(df, method = "euclidean")

# Compute 2 hierarchical clusterings
hc1 <- hclust(res.dist, method = "complete")
hc2 <- hclust(res.dist, method = "ward.D2")

# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)

# get tanglegram
# note that solid lines are cluster relationships in both dendrograms (common subtrees)
# dotted lines show cluster relationships that are not present in the other dendrogram; 
tanglegram(dend1, dend2)

# can also use tanglegram to get the entanglement similarity score between two dendrograms
dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_lines = FALSE, # Turn-off line colors
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2))
)


#//////////////////


# get optimal number of clusters ####

# note wss is within-cluster sum of squares, and just sums the euclidean distance btw each point in cluster and its centroid
# the centroid is just a vector with the average value for each variable over all observations in the cluster

# note the gap statistic is recently developed by tibshirani et al
fviz_nbclust(df, FUN = hcut, method = "wss")
fviz_nbclust(df, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(df, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)


#////////////////////


# use factoextra convenience functions ####

# get cluster tendency to see if data is likely to have clusters that are more meaningful/tight than uniform random sample
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# note that the hopkins stat uses .5 as a threshold, with higher values indicating there is a clustering tendency; 
# a hopkins stat over .75 is 90% confidence that there is clustering tendency that is real, and not random
clust_tendency <- get_clust_tendency(df, n = nrow(df)-1, graph = FALSE)
clust_tendency
clust_tendency$hopkins_stat

# use factoextra to create and cut clusters, so the factoextra object can be used in enhanced plots
f_clust <- hcut(USArrests, k = 4, stand = TRUE)
f_clust
f_clust %>% attributes()

# visualize dendrogram
fviz_dend(f_clust, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

# view PCA for clusters
fviz_cluster(object = f_clust)

# view silhouette plot
# closer to 1 is tighter cluster
fviz_silhouette(sil.obj = f_clust)


#/////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////


# modified ggdendro ####

# https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/


# load custom functions

dendro_data_k <- function(hc, k) {
        
        hcdata    <-  ggdendro::dendro_data(hc, type = "rectangle")
        seg       <-  hcdata$segments
        labclust  <-  cutree(hc, k)[hc$order]
        segclust  <-  rep(0L, nrow(seg))
        heights   <-  sort(hc$height, decreasing = TRUE)
        height    <-  mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
        
        for (i in 1:k) {
                xi      <-  hcdata$labels$x[labclust == i]
                idx1    <-  seg$x    >= min(xi) & seg$x    <= max(xi)
                idx2    <-  seg$xend >= min(xi) & seg$xend <= max(xi)
                idx3    <-  seg$yend < height
                idx     <-  idx1 & idx2 & idx3
                segclust[idx] <- i
        }
        
        idx                    <-  which(segclust == 0L)
        segclust[idx]          <-  segclust[idx + 1L]
        hcdata$segments$clust  <-  segclust
        hcdata$segments$line   <-  as.integer(segclust < 1L)
        hcdata$labels$clust    <-  labclust
        
        hcdata
}

set_labels_params <- function(nbLabels,
                              direction = c("tb", "bt", "lr", "rl"),
                              fan       = FALSE) {
        if (fan) {
                angle       <-  360 / nbLabels * 1:nbLabels + 90
                idx         <-  angle >= 90 & angle <= 270
                angle[idx]  <-  angle[idx] + 180
                hjust       <-  rep(0, nbLabels)
                hjust[idx]  <-  1
        } else {
                angle       <-  rep(0, nbLabels)
                hjust       <-  0
                if (direction %in% c("tb", "bt")) { angle <- angle + 45 }
                if (direction %in% c("tb", "rl")) { hjust <- 1 }
        }
        list(angle = angle, hjust = hjust, vjust = 0.5)
}

plot_ggdendro <- function(hcdata,
                          direction   = c("lr", "rl", "tb", "bt"),
                          fan         = FALSE,
                          scale.color = NULL,
                          branch.size = 1,
                          label.size  = 3,
                          nudge.label = 0.01,
                          expand.y    = 0.1) {
        
        direction <- match.arg(direction) # if fan = FALSE
        ybreaks   <- pretty(segment(hcdata)$y, n = 5)
        ymax      <- max(segment(hcdata)$y)
        
        ## branches
        p <- ggplot() +
                geom_segment(data         =  segment(hcdata),
                             aes(x        =  x,
                                 y        =  y,
                                 xend     =  xend,
                                 yend     =  yend,
                                 linetype =  factor(line),
                                 colour   =  factor(clust)),
                             lineend      =  "round",
                             show.legend  =  FALSE,
                             size         =  branch.size)
        
        ## orientation
        if (fan) {
                p <- p +
                        coord_polar(direction = -1) +
                        scale_x_continuous(breaks = NULL,
                                           limits = c(0, nrow(label(hcdata)))) +
                        scale_y_reverse(breaks = ybreaks)
        } else {
                p <- p + scale_x_continuous(breaks = NULL)
                if (direction %in% c("rl", "lr")) {
                        p <- p + coord_flip()
                }
                if (direction %in% c("bt", "lr")) {
                        p <- p + scale_y_reverse(breaks = ybreaks)
                } else {
                        p <- p + scale_y_continuous(breaks = ybreaks)
                        nudge.label <- -(nudge.label)
                }
        }
        
        # labels
        labelParams <- set_labels_params(nrow(hcdata$labels), direction, fan)
        hcdata$labels$angle <- labelParams$angle
        
        p <- p +
                geom_text(data        =  label(hcdata),
                          aes(x       =  x,
                              y       =  y,
                              label   =  label,
                              colour  =  factor(clust),
                              angle   =  angle),
                          fontface = "bold",
                          vjust       =  labelParams$vjust,
                          hjust       =  labelParams$hjust,
                          nudge_y     =  ymax * nudge.label,
                          size        =  label.size,
                          show.legend =  FALSE)
        
        # colors and limits
        if (!is.null(scale.color)) {
                p <- p + scale_color_manual(values = scale.color)
        }
        
        ylim <- -round(ymax * expand.y, 1)
        p    <- p + expand_limits(y = ylim)
        
        p
}


#/////////////////////


# example
mtc <- scale(mtcars)
D   <- dist(mtc)
hc  <- hclust(D)

hcdata <- dendro_data_k(hc, 3)

# inspect
hcdata
hcdata %>% attributes()
hcdata$segments %>% head()
hcdata$labels %>% head()
hcdata$leaf_labels %>% head()
hcdata$class %>% head()


#///////////////////////


# basic dendrogram
p <- plot_ggdendro(hcdata,
                   direction   = "lr",
                   expand.y    = 0.2)
p


#///////////////////////


# customized

# get tree colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
hcdata$labels
cols <- c("#a9a9a9", "#ff7f0e", "#1f77b4", "#2ca02c")

p <- plot_ggdendro(hcdata,
                   direction   = "tb",
                   scale.color = cols,
                   label.size  = 2.5,
                   branch.size = 0.5,
                   expand.y    = 0.2)

p <- p + theme_void() + expand_limits(x = c(-1, 32))
p


#///////////////////


# alternate customization

# get tree colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
hcdata$labels
cols <- c("#1f77b4", "#ff00ff", "#ff7f0e", "#2ca02c")


p <- plot_ggdendro(hcdata,
                   direction   = "tb",
                   scale.color = cols,
                   label.size  = 2.5,
                   branch.size = 0.5,
                   expand.y    = 0.2) +
        theme(axis.text = element_text(color = "#50505030"),
                 panel.grid.major.y = element_line(color = "#50505030",
                                                   size  = 0.25))
p


#/////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////////////////////


# cluster trajectories ####

# https://cran.r-project.org/web/packages/traj/vignettes/trajVignette.pdf

# each row is an individual trajectory, with the columns being the time observation
# note sure why the second dataframe $time is needed???
head(example.data$data)
head(example.data$time)

# step 1: get 24 summary measures for the trajectories
s1 = step1measures(example.data$data, example.data$time, ID = TRUE)
s1

# step 2: use factor analysis to get subset of 24 measures that best describes the main features of the trajectories
s2 = step2factors(s1)
head(s2$factors)

# step 3: cluster based on the selected subset of trajectory measures
# note you can pre-specify the number of clusters, or leave blank to get optimal number
s3 = step3clusters(s2, nclusters = 4)
s3 = step3clusters(s2)
s3

# inspect clusters
head(s3$clusters)
s3$clust.distr

# plot trajectory clusters
?plot.traj
plot(s3, num.samples = 10)
plot(s3, num.samples = 2, color.vect = c("#ff0000", "#00ff00"))
plot(s3, num.samples = 10, clust.num = 2)

# plot trajectory mean
?plotMeanTraj
?plotCom
plotMeanTraj(s3, clust.num = NULL)
plotCombTraj(s3)

# plot trajectory distribution
?plotMedTraj
plotMedTraj(s3)
plotBoxplotTraj(s3)





