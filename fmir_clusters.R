extrafont::loadfonts(device="win") 
library(tidyverse)
library(lubridate)
library(readxl)
library(skimr)
library(rlang)
library(haven)
library(fs)
library(officer)
library(devEMF)
library(ggrepel)
library(rvest)
library(scales)
library(testthat)
library(patchwork)
library(vdemdata)
library(naniar)
library(dichromat)
library(RColorBrewer)
library(fmsb)
library(ggcorrplot)
library(bazar)
library(cowplot)
library(ggplotify)
library(elementalist)
library(ggforce)
library(magick)
library(ggpubr)
library(grid)
library(rworldxtra)
library(sf)        
library(lwgeom)    
library(rworldmap)
library(pals)
library(openxlsx)
library(viridis)
library(rpart)
library(rattle)				
library(rpart.plot)
library(party)					
library(partykit)				
library(caret)
library(cluster)
library(factoextra)
library(dendextend)
library(traj)
library(ggdendro)


# https://uc-r.github.io/hc_clustering
# https://uc-r.github.io/kmeans_clustering
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# https://cran.r-project.org/web/packages/traj/vignettes/trajVignette.pdf
# includes definition of agglomerative coefficient: http://strata.uga.edu/8370/lecturenotes/clusterAnalysis.html
# ac: agglomerative coefficient, a measure of the clustering structure. The agglomerative coefficient measures the dissimilarity 
# of an object to the first cluster it joins, divided by the dissimilarity of the final merger in the cluster analysis, 
# averaged across all samples. Low values reflect tight clustering of objects, larger values indicate less well-formed clusters. 
# The agglomerative coefficient increases with sample sizes, making comparisons among data sets difficult.
# https://bradleyboehmke.github.io/HOML/hierarchical.html
# wards method: https://www.statisticshowto.com/wards-method/
# wards method: https://online.stat.psu.edu/stat505/lesson/14/14.7
# how to get centroid avg: https://sciencing.com/centroid-clustering-analysis-10070345.html
# https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#a-dendrogram-is-a-nested-list-of-lists-with-attributes
# https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
# https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html


# setwd
setwd("C:/Users/sdevine/Desktop/usaid/mcp/malign_influence")
options(scipen = 999)


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# get align_legend() from claus wilke ####

# https://stackoverflow.com/questions/48000292/center-align-legend-title-and-legend-keys-in-ggplot2-for-long-legend-titles

align_legend <- function(p, hjust = 0.5) {
        # extract legend
        g <- cowplot::plot_to_gtable(p)
        grobs <- g$grobs
        legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
        legend <- grobs[[legend_index]]
        
        # extract guides table
        guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
        
        # there can be multiple guides within one legend box  
        for (gi in guides_index) {
                guides <- legend$grobs[[gi]]
                
                # add extra column for spacing
                # guides$width[5] is the extra spacing from the end of the legend text
                # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
                # both sides, we get an aligned legend
                spacing <- guides$width[5]
                guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
                guides$widths[6] <- (1-hjust)*spacing
                title_index <- guides$layout$name == "title"
                guides$layout$l[title_index] <- 2
                
                # reconstruct guides and write back
                legend$grobs[[gi]] <- guides
        }
        
        # reconstruct legend and write back
        g$grobs[[legend_index]] <- legend
        g
}


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# load modified ggdendro functions
current_wd <- getwd()
setwd("C:/Users/sdevine/Desktop/R/clustering")
source("modified_ggdendro_functions.R")
setwd(current_wd)


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# load add_group_index() ####
current_wd <- getwd()
setwd("C:/Users/sdevine/Desktop/R/assorted_helper_scripts")
source("add_group_index.R")
setwd(current_wd)


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# create custom color_palette ####
color_palette <- tibble(hex = c("#083D7F", "#2474B6", "#8BBFD0",
                                "#CBCBCB", "#7D7D7D",
                                "#99ba78", "#35B779FF", "#006629", 
                                "#E4DC68", "#FDA159", "#EF6712", "#CE1B1E",
                                "#8B008B", "#DA70D6"))
color_palette
color_palette %>% pull(hex) %>% show_col()

# color_palette supports 11 colors, plus possible extensions via fill/line type
show_col(color_palette %>% slice(1, 3) %>% pull(hex)) # 2 colors
show_col(color_palette %>% slice(1, 2, 3) %>% pull(hex)) # 3 colors
show_col(color_palette %>% slice(1, 2, 3, 4) %>% pull(hex)) # 4 colors
show_col(color_palette %>% slice(1, 2, 3, 4, 5) %>% pull(hex)) # 5 colors
show_col(color_palette %>% slice(1, 2, 3, 4, 5, 6) %>% pull(hex)) # 6 colors
show_col(color_palette %>% slice(1, 2, 3, 4, 5, 6, 7) %>% pull(hex)) # 7 colors
show_col(color_palette %>% slice(1, 2, 3, 4, 5, 6, 7, 8) %>% pull(hex)) # 8 colors
show_col(color_palette %>% slice(1, 2, 3, 4, 5, 6, 7, 8, 9) %>% pull(hex)) # 9 colors
show_col(color_palette %>% slice(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) %>% pull(hex)) # 10 colors
show_col(color_palette %>% slice(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) %>% pull(hex)) # 11 colors


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# read country_crosswalk ####
current_wd <- getwd()
setwd("C:/Users/sdevine/Desktop/usaid/mcp/useful_info/country_crosswalk")
country_crosswalk <- read_csv("country_crosswalk.csv", lazy = FALSE)
setwd(current_wd)


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# read final fmir data ####
fmir <- read.csv(file = "data/fmir/fmir_20220307.csv") %>% as_tibble()

# inspect
fmir
fmir %>% glimpse()
fmir %>% nrow() # 30745
fmir %>% ncol() # 63


fmir %>% count(obj)
fmir %>% count(obj_num)
fmir %>% count(obj_short_name)
fmir %>% count(sub_obj)
fmir %>% count(sub_obj_num)
fmir %>% count(sub_obj_short_name)
fmir %>% count(concept)
fmir %>% count(indicator_name) %>% print(n = nrow(.))


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# create get_label_position() for use in facet charts ####
get_label_position <- function(chart_data) {
        
        # get initial label_position
        chart_data <- chart_data %>% arrange(year, values) %>% 
                group_by(year) %>%
                mutate(label_position = round(values, digits = 3),
                       dist_to_lower = round(label_position - lag(label_position, n = 1), digits = 3),
                       dist_to_higher = round(lead(label_position, n = 1) - label_position, digits = 3),
                       overlap_flag = case_when(dist_to_lower < .04 ~ 1,
                                                TRUE ~ 0)) %>%
                ungroup()
        
        # set counter
        counter <- 0
        
        # run loop to update label_position
        while((chart_data %>% mutate(sum_overlap_flag = sum(overlap_flag)) %>% slice(1) %>% pull(sum_overlap_flag)) != 0) {
                
                counter <- counter + 1
                print(counter + 1)
                
                chart_data <- chart_data %>% arrange(year, values) %>% 
                        group_by(year) %>%
                        mutate(label_position = case_when(is.na(dist_to_lower) | 
                                                                  dist_to_lower >= .04 ~ round(label_position, digits = 3),
                                                          TRUE ~ round(label_position + (.04 - dist_to_lower), digits = 3)),
                               dist_to_higher = round(lead(label_position, n = 1) - label_position, digits = 3),
                               dist_to_lower = round(label_position - lag(label_position, n = 1), digits = 3),
                               overlap_flag = case_when(dist_to_lower < .04 ~ 1,
                                                        TRUE ~ 0)) %>%
                        ungroup() %>%
                        arrange(desc(year), values)
        }
        
        # return
        return(chart_data)
}



#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# create round_to_digits() #### 
round_to_digits <- function(x, digits, big_mark = ",") {
        sign <- sign(x)
        z <- abs(x) * 10^digits
        z <- z + 0.5
        z <- trunc(z)
        z <- z / 10^digits
        z <- z * sign
        return(str_trim(string = format(x = z, nsmall = digits, big.mark = big_mark), side = "both"))
}


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_clusters ####

# get fmir_obj
fmir_obj <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, obj_num, obj_avg) %>%
        pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
        column_to_rownames(var = "country") %>%
        mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1]))


#//////////////////////


# inspect
fmir_obj
fmir_obj %>% glimpse()
fmir_obj %>% nrow()
fmir_obj %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get fmir_obj_distance_matrix 
fmir_obj_distance_matrix <- dist(fmir_obj, method = "euclidean")
# fmir_obj_distance_matrix <- get_dist(fmir_obj, method = "euclidean")


#//////////////


# inspect
fmir_obj_distance_matrix
fmir_obj_distance_matrix %>% str()
fmir_obj_distance_matrix %>% attributes()

# inspect fmir_obj_distance_matrix_tbl
fmir_obj_distance_matrix_tbl <- as.matrix(fmir_obj_distance_matrix) %>% as_tibble() %>% 
        mutate(country = row.names(as.matrix(fmir_obj_distance_matrix))) %>%
        relocate(country, .before = everything())

fmir_obj_distance_matrix_tbl
fmir_obj_distance_matrix_tbl %>% skim() %>% as_tibble() %>% 
        select(skim_variable, numeric.mean, starts_with(match = "numeric.p")) %>%
        # arrange(numeric.mean) %>% 
        arrange(numeric.p25) %>%
        print(n = nrow(.))
fmir_obj_distance_matrix_tbl %>% 
        # select(country, Serbia) %>% arrange(Serbia) %>% print(n = nrow(.))
        # select(country, BiH) %>% arrange(BiH) %>% print(n = nrow(.))
        # select(country, Kosovo) %>% arrange(Kosovo) %>% print(n = nrow(.))
        # select(country, Armenia) %>% arrange(Armenia) %>% print(n = nrow(.))
        select(country, Ukraine) %>% arrange(Ukraine) %>% print(n = nrow(.))
# filter(country == "Serbia") %>% select(ends_with("stan"))


# compare indicator_standardized_values
fmir %>% filter(obj_num == "obj_3", country %in% c("Poland", "Kazakhstan"), year == 2020) %>% 
        select(country, year, indicator_name, values, values_z_std, indicator_standardized_values) %>%
        arrange(country) %>% 
        # unite(col = var, indicator_name, country) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = indicator_name, names_from = country, values_from = indicator_standardized_values)


# good example of value of clustering instead of just looking at obj_avg
# based on obj_avg for obj_1: armenia (.42), georgia (.41), moldova (.41), ukraine (.36)
# so georgia looks very similar to armenia and moldova; but based on indicator-level distance it's 
# most similar to ukraine (dist = 3.13), then moldova (dist = .4.25), then armenia (5.57)
# note that indicator-level distance does not reflect implicit weighting based on indicator count per concept/sub-obj
# but it nonetheless captures similarity better, which could be useful info for regional programming/planning
fmir_obj_distance_matrix_tbl %>% 
        filter(country == "Ukraine") %>% select(Georgia, Armenia, Moldova)
fmir_obj_distance_matrix_tbl %>% 
        filter(country == "Armenia") %>% select(Georgia, Ukraine, Moldova)
fmir_obj_distance_matrix_tbl %>% 
        filter(country == "Georgia") %>% select(Ukraine, Armenia, Moldova)
fmir_obj_distance_matrix_tbl %>% 
        filter(country == "Moldova") %>% select(Ukraine, Armenia, Georgia)


# visualize distance
viridis_pal()(3) %>% show_col()
fviz_dist(dist.obj = fmir_obj_distance_matrix, 
          gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_distance_matrix_tile_chart ####


# get country_order from fviz_dist_plot
fmir_obj_fviz_dist_plot <- fviz_dist(dist.obj = fmir_obj_distance_matrix, 
                                         gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))
fmir_obj_fviz_dist_plot
fmir_obj_fviz_dist_plot %>% attributes()
fmir_obj_fviz_dist_plot$data %>% as_tibble()

fmir_obj_country_order_tbl <- fmir_obj_fviz_dist_plot$data %>% as_tibble() %>% slice(1:(fmir_obj %>% nrow())) %>% 
        mutate(country = str_sub(string = Var1, start = 1, end = -2),
               country_order = row_number()) %>% 
        select(country, country_order)
fmir_obj_country_order_tbl


#/////////////////////////////


# create chart
fmir_obj_distance_matrix_tile_chart <- fmir_obj_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "var", values_to = "values") %>%
        left_join(., fmir_obj_country_order_tbl, by = c("var" = "country")) %>%
        rename(var_order = country_order) %>%
        left_join(., fmir_obj_country_order_tbl, by = "country") %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = country, .x = country_order, .fun = min), 
                                       y = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Dissimiliarity\n(higher is more dissimilar)", y = NULL, x = NULL) +
        coord_fixed(ratio = 1/1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 13, r = 0, b = 5, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_distance_matrix_tile_chart <- ggdraw(align_legend(fmir_obj_distance_matrix_tile_chart))


# inspect
fmir_obj_distance_matrix_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_distance_matrix_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_distance_matrix_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get cluster tendency to see if data is likely to have clusters that are more meaningful/tight than uniform random sample
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# note that the hopkins stat uses .5 as a threshold, with higher values indicating there is a clustering tendency; 
# a hopkins stat over .75 is 90% confidence that there is clustering tendency that is real, and not random
cluster_tendency <- get_clust_tendency(fmir_obj, n = nrow(fmir_obj) - 1, graph = FALSE)


#///////////////////


# inspect
cluster_tendency
cluster_tendency$hopkins_stat


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters
fmir_obj_hclust <- hclust(d = fmir_obj_distance_matrix, method = "complete" )


#/////////////////////


# inspect
fmir_obj_hclust
fmir_obj_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_hclust$merge
fmir_obj_hclust$height

# plot dendrogram
plot(fmir_obj_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters with agnes, which allows enhanced fviz visualizations
fmir_obj_agnes <- agnes(x = fmir_obj, method = "complete")


#//////////////////


# inspect
fmir_obj_agnes
fmir_obj_agnes %>% attributes()

# get agglomerative coefficient from agnes object 
# "Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 suggest a 
# more balanced clustering structure such as the complete linkage and Wards method dendrograms in Figure 21.3. 
# Values closer to 0 suggest less well-formed clusters such as the single linkage dendrogram in Figure 21.3. 
# However, the AC tends to become larger as n increases, so it should not be used to compare across data sets of very different sizes."
fmir_obj_agnes$ac

# plot
pltree(fmir_obj_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 



#//////////////////////////////////////////////////////////////////////////////////////////////////


# use agnes to look for clustering method with best agglomerative coefficient 
method_list <- c( "average", "single", "complete", "ward")
names(method_list) <- c( "average", "single", "complete", "ward")
method_list

# function to compute coefficient
cluster_methods <- function(current_method) {
        print(current_method)
        agnes(fmir_obj_distance_matrix, method = current_method)$ac
}

# run cluster_methods
# note that wards method has the highest/best agglomerative coefficient
map(.x = method_list, .f = ~ cluster_methods(current_method = .x))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with agnes, which had highest (best) agglomerative coefficient
fmir_obj_agnes <- agnes(fmir_obj, method = "ward")
fmir_obj_agnes


#///////////////////


# visualize clustering 
pltree(fmir_obj_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
fviz_dend(fmir_obj_agnes, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

# inspect optimal number of clusters
# note the gap statistic is recently developed by tibshirani et al
fviz_nbclust(fmir_obj, FUN = hcut, method = "wss")
fviz_nbclust(fmir_obj, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(fmir_obj, FUN = hcut, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# set optimal k
fmir_obj_k <- 4


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with hclust
fmir_obj_hclust <- hclust(d = fmir_obj_distance_matrix, method = "ward.D2")


#/////////////////////


# inspect
fmir_obj_hclust
fmir_obj_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_hclust$merge
fmir_obj_hclust$height

# plot dendrogram
plot(fmir_obj_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut tree 

# build and cut agnes tree with hcut for enhanced fviz
fmir_obj_agnes_cut <- hcut(fmir_obj, k = fmir_obj_k, hc_method = "ward.D2", stand = TRUE)
fmir_obj_agnes_cut 

# view PCA for clusters, which requires an hcut object
fviz_cluster(object = fmir_obj_agnes_cut)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut hclust tree with cutree
fmir_obj_hclust_cut <- cutree(fmir_obj_hclust, k = fmir_obj_k)
fmir_obj_hclust_cut

# Number of members in each cluster
table(fmir_obj_hclust_cut)

# can add cutree output to df
fmir_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value) %>%
        left_join(fmir_obj %>% rownames_to_column(var = "country"), ., by = "country")

# plot dendrogram with border around 4 clusters
plot(fmir_obj_hclust, cex = 0.6, hang = -1)
rect.hclust(fmir_obj_hclust, k = fmir_obj_k, border = 2:5)

# plot dendrogram w labels
dendrogram_obj <- as.dendrogram(fmir_obj_hclust)
dendrogram_obj %>% unclass() %>% str()
dendrogram_obj %>% unclass() %>% class()

subdendrogram_list <- get_subdendrograms(dendrogram_obj, k = fmir_obj_k)
subdendrogram_list
dendrogram_labels <- map(.x = subdendrogram_list, .f = ~ labels(.x)) %>% 
        enframe() %>% unnest(cols = value) %>% 
        rename(subdendrogram_plot_order = name, country = value) %>%
        left_join(fmir_obj %>% mutate(cluster = fmir_obj_hclust_cut) %>%
                          rownames_to_column(var = "country"),
                  ., by = "country") %>% 
        distinct(cluster, subdendrogram_plot_order) %>% 
        arrange(subdendrogram_plot_order) %>%
        mutate(cluster = str_c("\n\n\n", cluster)) %>%
        pull(cluster) %>% unname()
dendrogram_labels        

plot(fmir_obj_hclust, hang = -1)
rect.dendrogram(as.dendrogram(fmir_obj_hclust), k = fmir_obj_k, text = dendrogram_labels)

# visualize clusters along 2 primary PCA dimensions
fviz_cluster(list(data = fmir_obj, cluster = fmir_obj_hclust_cut))


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_dendrogram ####

# get fmir_obj_dendrogram_data
# note that dendro_data_k() takes the uncut hclust as first argument
fmir_obj_dendrogram_data <- dendro_data_k(hc = fmir_obj_hclust, 
                                              k = fmir_obj_k)

# inspect
fmir_obj_dendrogram_data
fmir_obj_dendrogram_data %>% attributes()
fmir_obj_dendrogram_data$segments %>% head()
fmir_obj_dendrogram_data$labels %>% head()
fmir_obj_dendrogram_data$leaf_labels %>% head()
fmir_obj_dendrogram_data$class %>% head()


#///////////////////


# get cluster_colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
cluster_colors <- c("#333333", 
                    color_palette %>% slice(2) %>% pull(hex), 
                    color_palette %>% slice(7) %>% pull(hex),
                    color_palette %>% slice(10) %>% pull(hex),
                    color_palette %>% slice(13) %>% pull(hex))
cluster_colors


# plot
fmir_obj_dendrogram <- plot_ggdendro(fmir_obj_dendrogram_data,
                                         direction = "tb",
                                         scale.color = cluster_colors,
                                         label.size = 2.5,
                                         branch.size = 0.5,
                                         expand.y = 0) +
        labs(y = "Dissimilarity\n(higher groupings are more dissimilar)", x = "") +
        coord_fixed(ratio = 1 / 1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 5, 0, 0), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(size = .5),
                axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "bold", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 0, hjust = .5),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                # axis.title.x = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 5, r = 0, b = 10, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 0, r = 15, b = 0, l = 5)),
                # axis.title.y = element_blank(),
                # axis.text.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 9, family = "Calibri", face = "plain", color = "#333333"),
                legend.text = element_text(size = 9, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# guides(color = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333", keywidth = 4),
#        linetype = guide_legend(keywidth = 4))

# inspect
fmir_obj_dendrogram


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_dendrogram)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_dendrogram.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_country_indicator_tile_chart ####

# get chart_data
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, obj_num, obj_avg) %>%
        left_join(., fmir_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        left_join(., fmir_obj_country_order_tbl, by = "country") %>%
        relocate(country_order, .after = cluster) %>% arrange(country_order) %>%
        pivot_wider(id_cols = c(country, cluster, country_order), 
                    names_from = obj_num, values_from = obj_avg) %>%
        pivot_longer(cols = -c(country, cluster, country_order), names_to = "var", values_to = "values") %>%
        mutate(color = case_when(cluster == 1 ~ color_palette %>% slice(2) %>% pull(hex),
                                 cluster == 2 ~ color_palette %>% slice(7) %>% pull(hex),
                                 cluster == 3 ~ color_palette %>% slice(10) %>% pull(hex),
                                 cluster == 4 ~ color_palette %>% slice(13) %>% pull(hex)),
               country_label = str_c("Group ", cluster, ": ", country),
               var_order = var)

# inspect
chart_data


#/////////////////////////////


# create chart
fmir_obj_country_indicator_tile_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       y = fct_reorder(.f = country_label, .x = country_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Score\n(higher is better)\n", y = NULL, x = "Indicator") +
        coord_fixed(ratio = 1 / 3, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 5, 
                                           color = chart_data %>% distinct(country, color) %>% pull(color), 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 5, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 5, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_country_indicator_tile_chart <- ggdraw(align_legend(fmir_obj_country_indicator_tile_chart))


# inspect
fmir_obj_country_indicator_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_country_indicator_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_country_indicator_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_cluster_boxplot_chart ####

fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, sub_obj_num, sub_obj_avg) %>%
        left_join(., fmir_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country")

# get chart_data
# note that this uses the original ISV values from the index, not the scaled ISV used in fmir_obj for clustering; 
# original/unscaled ISV avg is better for plot because values are then interpretable in the framework of the index (e.g. ISV btw 0-1)
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, sub_obj_num, sub_obj_avg) %>%
        pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
        
        
        # add clusters
        left_join(., fmir_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15") %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          distinct(country, year, obj_num, obj_avg) %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg),
                  ., by = "country") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -c(country, cluster), names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster),
               color_bin = as.character(cluster_name),
               color = case_when(color_bin == "Group 1" ~ color_palette %>% slice(2) %>% pull(hex),
                                 color_bin == "Group 2" ~ color_palette %>% slice(7) %>% pull(hex),
                                 color_bin == "Group 3" ~ color_palette %>% slice(10) %>% pull(hex)))

# inspect
chart_data
chart_data %>% glimpse()

# create color_list for to pass to scale_color_manual
chart_data_color_list <- chart_data %>% count(color_bin, color) %>% pull(color)
names(chart_data_color_list) <- chart_data %>% count(color_bin, color) %>% pull(color_bin)
chart_data_color_list


#//////////////////////////


# create chart
fmir_obj_cluster_boxplot_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = values, fill = color_bin)) + 
        geom_boxplot(width = .5, lwd = .01, outlier.size = .5) + facet_wrap(facets = vars(cluster_name)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = .2), limits = c(0, 1), expand = c(0, 0)) +
        scale_fill_manual(values = chart_data_color_list) + 
        # scale_color_manual(values = chart_data_color_list) + 
        coord_fixed(ratio = 10/.6, clip = "off") +
        labs(x = "Indicator", y = "Score (higher is better)", fill = NULL, color = NULL) +
        theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                # strip.background = element_blank(),
                strip.background = element_rect(fill = "#DDDDDD", color = NA),
                strip.text = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333"),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = -15, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 5, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = 1,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 


# inspect
fmir_obj_cluster_boxplot_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_cluster_boxplot_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_cluster_boxplot_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_cluster_tile_chart ####

# get chart_data 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, sub_obj_num, sub_obj_avg) %>%
        pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
        
        
        # add clusters
        left_join(., fmir_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        group_by(cluster) %>%
        mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
        ungroup() %>%
        select(-country) %>% distinct() %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15") %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -cluster, names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster))

# inspect
chart_data


#///////////////


# get obj/sub_obj segment coordinates

# note that in the entire fmir dataset, the maximum sub_obj per obj is 3, so will have 3 segment_tbls,
# and for those obj with less than 3 sub_obj, the segment_tbl will just have NA coordinates and no output in ggplot
fmir %>% distinct(obj_num, sub_obj_num) %>% 
        add_count(obj_num, name = "sub_obj_count") %>% 
        distinct(obj_num, sub_obj_count) %>%
        arrange(desc(sub_obj_count)) 

# segment_1_tbl
# segment_1_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(obj_num, sub_obj_num) %>%
#         mutate(nrow = nrow(.),
#                x = nrow + 1.5,
#                xend = nrow + 1.5,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5) %>%
#         select(x, xend, y, yend) %>%
#         slice(1)
# segment_1_tbl

# # segment_2_tbl
# segment_2_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                 cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(1)
# segment_2_tbl
#                
# # segment_3_tbl
# segment_3_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(2)
# segment_3_tbl


#///////////////////////


# plot

# get order_obj_and_sub_obj_avg_first()
order_obj_and_sub_obj_avg_first <- function(var){
        
        tibble(var = var) %>% 
                left_join(., fmir %>% distinct(indicator_name, concept) %>% 
                                  add_group_index(group_vars = concept, group_name = "concept_index"), 
                          by = c("var" = "indicator_name")) %>%
                mutate(order = case_when(str_detect(string = var, pattern = "^obj_") ~ 1,
                                         str_detect(string = var, pattern = "^sub_obj_[0-9]_[0-9]$") ~ 2,
                                         TRUE ~ 3)) %>%
                arrange(order, concept_index, var) %>% pull(var)
}

# create chart
fmir_obj_cluster_tile_chart <- chart_data %>% 
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = factor(cluster_name), fill = values, label = round_to_digits(values, digits = 2))) + 
        geom_tile(color = "#DDDDDD") + 
        geom_text(color = "#ffffff", size = 1.75, family = "Calibri") +
        scale_y_discrete(expand = c(0, 0)) +
        # geom_segment(x = segment_1_tbl$x, xend = segment_1_tbl$xend, 
        #              y = segment_1_tbl$y, yend = segment_1_tbl$yend, color = "#ffffff", size = 3) +
        # geom_segment(x = segment_2_tbl$x, xend = segment_2_tbl$xend,
        #              y = segment_2_tbl$y, yend = segment_2_tbl$yend, color = "#CBCBCB", size = 3) +
        # geom_segment(x = segment_3_tbl$x, xend = segment_3_tbl$xend,
        #              y = segment_3_tbl$y, yend = segment_3_tbl$yend, color = "#CBCBCB", size = 3) +
        labs(x = "Indicator", y = NULL, fill = "Avg. score\n(higher is better)\n") +
        coord_fixed(ratio = 2 / 1, clip = "off") +
        theme_bw() +
        # theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD"),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "right",
                legend.key.size = unit(4, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 5, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# fmir_obj_cluster_tile_chart <- ggdraw(align_legend(fmir_obj_cluster_tile_chart))

# inspect
fmir_obj_cluster_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_cluster_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_cluster_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_sub_obj_clusters ####

# get fmir_sub_obj
fmir_sub_obj <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, sub_obj_num, sub_obj_avg) %>%
        pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
        column_to_rownames(var = "country") %>%
        mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1]))


#//////////////////////


# inspect
fmir_sub_obj
fmir_sub_obj %>% glimpse()
fmir_sub_obj %>% nrow()
fmir_sub_obj %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get fmir_sub_obj_distance_matrix 
fmir_sub_obj_distance_matrix <- dist(fmir_sub_obj, method = "euclidean")
# fmir_sub_obj_distance_matrix <- get_dist(fmir_sub_obj, method = "euclidean")


#//////////////


# inspect
fmir_sub_obj_distance_matrix
fmir_sub_obj_distance_matrix %>% str()
fmir_sub_obj_distance_matrix %>% attributes()

# inspect fmir_sub_obj_distance_matrix_tbl
fmir_sub_obj_distance_matrix_tbl <- as.matrix(fmir_sub_obj_distance_matrix) %>% as_tibble() %>% 
        mutate(country = row.names(as.matrix(fmir_sub_obj_distance_matrix))) %>%
        relocate(country, .before = everything())

fmir_sub_obj_distance_matrix_tbl
fmir_sub_obj_distance_matrix_tbl %>% skim() %>% as_tibble() %>% 
        select(skim_variable, numeric.mean, starts_with(match = "numeric.p")) %>%
        # arrange(numeric.mean) %>% 
        arrange(numeric.p25) %>%
        print(n = nrow(.))
fmir_sub_obj_distance_matrix_tbl %>% 
        # select(country, Serbia) %>% arrange(Serbia) %>% print(n = nrow(.))
        # select(country, BiH) %>% arrange(BiH) %>% print(n = nrow(.))
        # select(country, Kosovo) %>% arrange(Kosovo) %>% print(n = nrow(.))
        # select(country, Armenia) %>% arrange(Armenia) %>% print(n = nrow(.))
        select(country, Ukraine) %>% arrange(Ukraine) %>% print(n = nrow(.))
# filter(country == "Serbia") %>% select(ends_with("stan"))


# compare indicator_standardized_values
fmir %>% filter(obj_num == "obj_3", country %in% c("Poland", "Kazakhstan"), year == 2020) %>% 
        select(country, year, indicator_name, values, values_z_std, indicator_standardized_values) %>%
        arrange(country) %>% 
        # unite(col = var, indicator_name, country) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = indicator_name, names_from = country, values_from = indicator_standardized_values)


# good example of value of clustering instead of just looking at obj_avg
# based on obj_avg for obj_1: armenia (.42), georgia (.41), moldova (.41), ukraine (.36)
# so georgia looks very similar to armenia and moldova; but based on indicator-level distance it's 
# most similar to ukraine (dist = 3.13), then moldova (dist = .4.25), then armenia (5.57)
# note that indicator-level distance does not reflect implicit weighting based on indicator count per concept/sub-obj
# but it nonetheless captures similarity better, which could be useful info for regional programming/planning
fmir_sub_obj_distance_matrix_tbl %>% 
        filter(country == "Ukraine") %>% select(Georgia, Armenia, Moldova)
fmir_sub_obj_distance_matrix_tbl %>% 
        filter(country == "Armenia") %>% select(Georgia, Ukraine, Moldova)
fmir_sub_obj_distance_matrix_tbl %>% 
        filter(country == "Georgia") %>% select(Ukraine, Armenia, Moldova)
fmir_sub_obj_distance_matrix_tbl %>% 
        filter(country == "Moldova") %>% select(Ukraine, Armenia, Georgia)


# visualize distance
viridis_pal()(3) %>% show_col()
fviz_dist(dist.obj = fmir_sub_obj_distance_matrix, 
          gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_sub_obj_distance_matrix_tile_chart ####


# get country_order from fviz_dist_plot
fmir_sub_obj_fviz_dist_plot <- fviz_dist(dist.obj = fmir_sub_obj_distance_matrix, 
                                       gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))
fmir_sub_obj_fviz_dist_plot
fmir_sub_obj_fviz_dist_plot %>% attributes()
fmir_sub_obj_fviz_dist_plot$data %>% as_tibble()

fmir_sub_obj_country_order_tbl <- fmir_sub_obj_fviz_dist_plot$data %>% as_tibble() %>% slice(1:(fmir_sub_obj %>% nrow())) %>% 
        mutate(country = str_sub(string = Var1, start = 1, end = -2),
               country_order = row_number()) %>% 
        select(country, country_order)
fmir_sub_obj_country_order_tbl


#/////////////////////////////


# create chart
fmir_sub_obj_distance_matrix_tile_chart <- fmir_sub_obj_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "var", values_to = "values") %>%
        left_join(., fmir_sub_obj_country_order_tbl, by = c("var" = "country")) %>%
        rename(var_order = country_order) %>%
        left_join(., fmir_sub_obj_country_order_tbl, by = "country") %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = country, .x = country_order, .fun = min), 
                                       y = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Dissimiliarity\n(higher is more dissimilar)", y = NULL, x = NULL) +
        coord_fixed(ratio = 1/1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 13, r = 0, b = 5, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_sub_obj_distance_matrix_tile_chart <- ggdraw(align_legend(fmir_sub_obj_distance_matrix_tile_chart))


# inspect
fmir_sub_obj_distance_matrix_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_sub_obj_distance_matrix_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_sub_obj_distance_matrix_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get cluster tendency to see if data is likely to have clusters that are more meaningful/tight than uniform random sample
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# note that the hopkins stat uses .5 as a threshold, with higher values indicating there is a clustering tendency; 
# a hopkins stat over .75 is 90% confidence that there is clustering tendency that is real, and not random
cluster_tendency <- get_clust_tendency(fmir_sub_obj, n = nrow(fmir_sub_obj) - 1, graph = FALSE)


#///////////////////


# inspect
cluster_tendency
cluster_tendency$hopkins_stat


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters
fmir_sub_obj_hclust <- hclust(d = fmir_sub_obj_distance_matrix, method = "complete" )


#/////////////////////


# inspect
fmir_sub_obj_hclust
fmir_sub_obj_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_sub_obj_hclust$merge
fmir_sub_obj_hclust$height

# plot dendrogram
plot(fmir_sub_obj_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters with agnes, which allows enhanced fviz visualizations
fmir_sub_obj_agnes <- agnes(x = fmir_sub_obj, method = "complete")


#//////////////////


# inspect
fmir_sub_obj_agnes
fmir_sub_obj_agnes %>% attributes()

# get agglomerative coefficient from agnes object 
# "Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 suggest a 
# more balanced clustering structure such as the complete linkage and Wards method dendrograms in Figure 21.3. 
# Values closer to 0 suggest less well-formed clusters such as the single linkage dendrogram in Figure 21.3. 
# However, the AC tends to become larger as n increases, so it should not be used to compare across data sets of very different sizes."
fmir_sub_obj_agnes$ac

# plot
pltree(fmir_sub_obj_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 



#//////////////////////////////////////////////////////////////////////////////////////////////////


# use agnes to look for clustering method with best agglomerative coefficient 
method_list <- c( "average", "single", "complete", "ward")
names(method_list) <- c( "average", "single", "complete", "ward")
method_list

# function to compute coefficient
cluster_methods <- function(current_method) {
        print(current_method)
        agnes(fmir_sub_obj_distance_matrix, method = current_method)$ac
}

# run cluster_methods
# note that wards method has the highest/best agglomerative coefficient
map(.x = method_list, .f = ~ cluster_methods(current_method = .x))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with agnes, which had highest (best) agglomerative coefficient
fmir_sub_obj_agnes <- agnes(fmir_sub_obj, method = "ward")
fmir_sub_obj_agnes


#///////////////////


# visualize clustering 
pltree(fmir_sub_obj_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
fviz_dend(fmir_sub_obj_agnes, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

# inspect optimal number of clusters
# note the gap statistic is recently developed by tibshirani et al
fviz_nbclust(fmir_sub_obj, FUN = hcut, method = "wss")
fviz_nbclust(fmir_sub_obj, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(fmir_sub_obj, FUN = hcut, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# set optimal k
fmir_sub_obj_k <- 3


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with hclust
fmir_sub_obj_hclust <- hclust(d = fmir_sub_obj_distance_matrix, method = "ward.D2")


#/////////////////////


# inspect
fmir_sub_obj_hclust
fmir_sub_obj_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_sub_obj_hclust$merge
fmir_sub_obj_hclust$height

# plot dendrogram
plot(fmir_sub_obj_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut tree 

# build and cut agnes tree with hcut for enhanced fviz
fmir_sub_obj_agnes_cut <- hcut(fmir_sub_obj, k = fmir_sub_obj_k, hc_method = "ward.D2", stand = TRUE)
fmir_sub_obj_agnes_cut 

# view PCA for clusters, which requires an hcut object
fviz_cluster(object = fmir_sub_obj_agnes_cut)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut hclust tree with cutree
fmir_sub_obj_hclust_cut <- cutree(fmir_sub_obj_hclust, k = fmir_sub_obj_k)
fmir_sub_obj_hclust_cut

# Number of members in each cluster
table(fmir_sub_obj_hclust_cut)

# can add cutree output to df
fmir_sub_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value) %>%
        left_join(fmir_sub_obj %>% rownames_to_column(var = "country"), ., by = "country")

# plot dendrogram with border around 4 clusters
plot(fmir_sub_obj_hclust, cex = 0.6, hang = -1)
rect.hclust(fmir_sub_obj_hclust, k = fmir_sub_obj_k, border = 2:5)

# plot dendrogram w labels
dendrogram_obj <- as.dendrogram(fmir_sub_obj_hclust)
dendrogram_obj %>% unclass() %>% str()
dendrogram_obj %>% unclass() %>% class()

subdendrogram_list <- get_subdendrograms(dendrogram_obj, k = fmir_sub_obj_k)
subdendrogram_list
dendrogram_labels <- map(.x = subdendrogram_list, .f = ~ labels(.x)) %>% 
        enframe() %>% unnest(cols = value) %>% 
        rename(subdendrogram_plot_order = name, country = value) %>%
        left_join(fmir_sub_obj %>% mutate(cluster = fmir_sub_obj_hclust_cut) %>%
                          rownames_to_column(var = "country"),
                  ., by = "country") %>% 
        distinct(cluster, subdendrogram_plot_order) %>% 
        arrange(subdendrogram_plot_order) %>%
        mutate(cluster = str_c("\n\n\n", cluster)) %>%
        pull(cluster) %>% unname()
dendrogram_labels        

plot(fmir_sub_obj_hclust, hang = -1)
rect.dendrogram(as.dendrogram(fmir_sub_obj_hclust), k = fmir_sub_obj_k, text = dendrogram_labels)

# visualize clusters along 2 primary PCA dimensions
fviz_cluster(list(data = fmir_sub_obj, cluster = fmir_sub_obj_hclust_cut))


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_sub_obj_dendrogram ####

# get fmir_sub_obj_dendrogram_data
# note that dendro_data_k() takes the uncut hclust as first argument
fmir_sub_obj_dendrogram_data <- dendro_data_k(hc = fmir_sub_obj_hclust, 
                                                      k = fmir_sub_obj_k)

# inspect
fmir_sub_obj_dendrogram_data
fmir_sub_obj_dendrogram_data %>% attributes()
fmir_sub_obj_dendrogram_data$segments %>% head()
fmir_sub_obj_dendrogram_data$labels %>% head()
fmir_sub_obj_dendrogram_data$leaf_labels %>% head()
fmir_sub_obj_dendrogram_data$class %>% head()


#///////////////////


# get cluster_colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
cluster_colors <- c("#333333", 
                    color_palette %>% slice(2) %>% pull(hex), 
                    color_palette %>% slice(7) %>% pull(hex),
                    color_palette %>% slice(10) %>% pull(hex))
cluster_colors


# plot
fmir_sub_obj_dendrogram <- plot_ggdendro(fmir_sub_obj_dendrogram_data,
                                                 direction = "tb",
                                                 scale.color = cluster_colors,
                                                 label.size = 2.5,
                                                 branch.size = 0.5,
                                                 expand.y = 0) +
        labs(y = "Dissimilarity\n(higher groupings are more dissimilar)", x = "") +
        coord_fixed(ratio = 1 / 1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 5, 0, 0), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(size = .5),
                axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "bold", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 0, hjust = .5),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                # axis.title.x = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 5, r = 0, b = 10, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 0, r = 15, b = 0, l = 5)),
                # axis.title.y = element_blank(),
                # axis.text.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 9, family = "Calibri", face = "plain", color = "#333333"),
                legend.text = element_text(size = 9, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# guides(color = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333", keywidth = 4),
#        linetype = guide_legend(keywidth = 4))

# inspect
fmir_sub_obj_dendrogram


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_sub_obj_dendrogram)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_sub_obj_dendrogram.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_sub_obj_country_indicator_tile_chart ####

# get chart_data
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, sub_obj_num, sub_obj_avg) %>%
        left_join(., fmir_sub_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        left_join(., fmir_sub_obj_country_order_tbl, by = "country") %>%
        relocate(country_order, .after = cluster) %>% arrange(country_order) %>%
        pivot_wider(id_cols = c(country, cluster, country_order), 
                    names_from = sub_obj_num, values_from = sub_obj_avg) %>%
        left_join(., fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15") %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          distinct(country, obj_num, obj_avg) %>%
                          pivot_wider(id_cols = c(country), names_from = obj_num, values_from = obj_avg), 
                  by = "country") %>%
        pivot_longer(cols = -c(country, cluster, country_order), names_to = "var", values_to = "values") %>%
        mutate(color = case_when(cluster == 1 ~ color_palette %>% slice(2) %>% pull(hex),
                                 cluster == 2 ~ color_palette %>% slice(7) %>% pull(hex),
                                 cluster == 3 ~ color_palette %>% slice(10) %>% pull(hex)),
               country_label = str_c("Group ", cluster, ": ", country),
               var_order = case_when(str_detect(string = var, pattern = "^obj_[0-9]|[a-z]$") ~ 1,
                                     TRUE ~ 2))

# inspect
chart_data


#/////////////////////////////


# create chart
fmir_sub_obj_country_indicator_tile_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       y = fct_reorder(.f = country_label, .x = country_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Score\n(higher is better)\n", y = NULL, x = "Indicator") +
        coord_fixed(ratio = 1 / 3, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 5, 
                                           color = chart_data %>% distinct(country, color) %>% pull(color), 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 5, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 5, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_sub_obj_country_indicator_tile_chart <- ggdraw(align_legend(fmir_sub_obj_country_indicator_tile_chart))


# inspect
fmir_sub_obj_country_indicator_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_sub_obj_country_indicator_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_sub_obj_country_indicator_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_sub_obj_cluster_boxplot_chart ####

fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, sub_obj_num, sub_obj_avg) %>%
        left_join(., fmir_sub_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country")

# get chart_data
# note that this uses the original ISV values from the index, not the scaled ISV used in fmir_sub_obj for clustering; 
# original/unscaled ISV avg is better for plot because values are then interpretable in the framework of the index (e.g. ISV btw 0-1)
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, sub_obj_num, sub_obj_avg) %>%
        pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
        
        
        # add clusters
        left_join(., fmir_sub_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15") %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          distinct(country, year, obj_num, obj_avg) %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg),
                  ., by = "country") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -c(country, cluster), names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster),
               color_bin = as.character(cluster_name),
               color = case_when(color_bin == "Group 1" ~ color_palette %>% slice(2) %>% pull(hex),
                                 color_bin == "Group 2" ~ color_palette %>% slice(7) %>% pull(hex),
                                 color_bin == "Group 3" ~ color_palette %>% slice(10) %>% pull(hex)))

# inspect
chart_data
chart_data %>% glimpse()

# create color_list for to pass to scale_color_manual
chart_data_color_list <- chart_data %>% count(color_bin, color) %>% pull(color)
names(chart_data_color_list) <- chart_data %>% count(color_bin, color) %>% pull(color_bin)
chart_data_color_list


#//////////////////////////


# create chart
fmir_sub_obj_cluster_boxplot_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = values, fill = color_bin)) + 
        geom_boxplot(width = .5, lwd = .01, outlier.size = .5) + facet_wrap(facets = vars(cluster_name)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = .2), limits = c(0, 1), expand = c(0, 0)) +
        scale_fill_manual(values = chart_data_color_list) + 
        # scale_color_manual(values = chart_data_color_list) + 
        coord_fixed(ratio = 10/.6, clip = "off") +
        labs(x = "Indicator", y = "Score (higher is better)", fill = NULL, color = NULL) +
        theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                # strip.background = element_blank(),
                strip.background = element_rect(fill = "#DDDDDD", color = NA),
                strip.text = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333"),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = -15, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 5, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = 1,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 


# inspect
fmir_sub_obj_cluster_boxplot_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_sub_obj_cluster_boxplot_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_sub_obj_cluster_boxplot_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_sub_obj_cluster_tile_chart ####

# get chart_data 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        distinct(country, sub_obj_num, sub_obj_avg) %>%
        pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
        
        
        # add clusters
        left_join(., fmir_sub_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        group_by(cluster) %>%
        mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
        ungroup() %>%
        select(-country) %>% distinct() %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15") %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_sub_obj_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -cluster, names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster))

# inspect
chart_data


#///////////////


# get obj/sub_obj segment coordinates

# note that in the entire fmir dataset, the maximum sub_obj per obj is 3, so will have 3 segment_tbls,
# and for those obj with less than 3 sub_obj, the segment_tbl will just have NA coordinates and no output in ggplot
fmir %>% distinct(obj_num, sub_obj_num) %>% 
        add_count(obj_num, name = "sub_obj_count") %>% 
        distinct(obj_num, sub_obj_count) %>%
        arrange(desc(sub_obj_count)) 

# segment_1_tbl
# segment_1_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(obj_num, sub_obj_num) %>%
#         mutate(nrow = nrow(.),
#                x = nrow + 1.5,
#                xend = nrow + 1.5,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5) %>%
#         select(x, xend, y, yend) %>%
#         slice(1)
# segment_1_tbl

# # segment_2_tbl
# segment_2_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                 cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(1)
# segment_2_tbl
#                
# # segment_3_tbl
# segment_3_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(2)
# segment_3_tbl


#///////////////////////


# plot

# get order_obj_and_sub_obj_avg_first()
order_obj_and_sub_obj_avg_first <- function(var){
        
        tibble(var = var) %>% 
                left_join(., fmir %>% distinct(indicator_name, concept) %>% 
                                  add_group_index(group_vars = concept, group_name = "concept_index"), 
                          by = c("var" = "indicator_name")) %>%
                mutate(order = case_when(str_detect(string = var, pattern = "^obj_") ~ 1,
                                         str_detect(string = var, pattern = "^sub_obj_[0-9]_[0-9]$") ~ 2,
                                         TRUE ~ 3)) %>%
                arrange(order, concept_index, var) %>% pull(var)
}

# create chart
fmir_sub_obj_cluster_tile_chart <- chart_data %>% 
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = factor(cluster_name), fill = values, label = round_to_digits(values, digits = 2))) + 
        geom_tile(color = "#DDDDDD") + 
        geom_text(color = "#ffffff", size = 1.75, family = "Calibri") +
        scale_y_discrete(expand = c(0, 0)) +
        # geom_segment(x = segment_1_tbl$x, xend = segment_1_tbl$xend, 
        #              y = segment_1_tbl$y, yend = segment_1_tbl$yend, color = "#ffffff", size = 3) +
        # geom_segment(x = segment_2_tbl$x, xend = segment_2_tbl$xend,
        #              y = segment_2_tbl$y, yend = segment_2_tbl$yend, color = "#CBCBCB", size = 3) +
        # geom_segment(x = segment_3_tbl$x, xend = segment_3_tbl$xend,
        #              y = segment_3_tbl$y, yend = segment_3_tbl$yend, color = "#CBCBCB", size = 3) +
        labs(x = "Indicator", y = NULL, fill = "Avg. score\n(higher is better)\n") +
        coord_fixed(ratio = 2 / 1, clip = "off") +
        theme_bw() +
        # theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD"),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "right",
                legend.key.size = unit(4, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 5, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# fmir_sub_obj_cluster_tile_chart <- ggdraw(align_legend(fmir_sub_obj_cluster_tile_chart))

# inspect
fmir_sub_obj_cluster_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_sub_obj_cluster_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_sub_obj_cluster_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_1_indicator_clusters ####

# set current_obj_num
current_obj_num <- "obj_1"

# get fmir_obj_1_indicator with scaled version of indicator_standardized_values
fmir_obj_1_indicator <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        distinct(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        column_to_rownames(var = "country") %>%
        mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1]))


#//////////////////////


# inspect
fmir_obj_1_indicator
fmir_obj_1_indicator %>% glimpse()
fmir_obj_1_indicator %>% nrow()
fmir_obj_1_indicator %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get fmir_obj_1_indicator_distance_matrix 
fmir_obj_1_indicator_distance_matrix <- dist(fmir_obj_1_indicator, method = "euclidean")
# fmir_obj_1_indicator_distance_matrix <- get_dist(fmir_obj_1_indicator, method = "euclidean")


#//////////////


# inspect
fmir_obj_1_indicator_distance_matrix
fmir_obj_1_indicator_distance_matrix %>% str()
fmir_obj_1_indicator_distance_matrix %>% attributes()

# inspect fmir_obj_1_indicator_distance_matrix_tbl
fmir_obj_1_indicator_distance_matrix_tbl <- as.matrix(fmir_obj_1_indicator_distance_matrix) %>% as_tibble() %>% 
        mutate(country = row.names(as.matrix(fmir_obj_1_indicator_distance_matrix))) %>%
        relocate(country, .before = everything())

fmir_obj_1_indicator_distance_matrix_tbl
fmir_obj_1_indicator_distance_matrix_tbl %>% skim() %>% as_tibble() %>% 
        select(skim_variable, numeric.mean, starts_with(match = "numeric.p")) %>%
        # arrange(numeric.mean) %>% 
        arrange(numeric.p25) %>%
        print(n = nrow(.))
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        # select(country, Serbia) %>% arrange(Serbia) %>% print(n = nrow(.))
        # select(country, BiH) %>% arrange(BiH) %>% print(n = nrow(.))
        # select(country, Kosovo) %>% arrange(Kosovo) %>% print(n = nrow(.))
        # select(country, Armenia) %>% arrange(Armenia) %>% print(n = nrow(.))
        select(country, Ukraine) %>% arrange(Ukraine) %>% print(n = nrow(.))
# filter(country == "Serbia") %>% select(ends_with("stan"))


# compare indicator_standardized_values
fmir %>% filter(obj_num == "obj_3", country %in% c("Poland", "Kazakhstan"), year == 2020) %>% 
        select(country, year, indicator_name, values, values_z_std, indicator_standardized_values) %>%
        arrange(country) %>% 
        # unite(col = var, indicator_name, country) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = indicator_name, names_from = country, values_from = indicator_standardized_values)


# good example of value of clustering instead of just looking at obj_avg
# based on obj_avg for obj_1: armenia (.42), georgia (.41), moldova (.41), ukraine (.36)
# so georgia looks very similar to armenia and moldova; but based on indicator-level distance it's 
# most similar to ukraine (dist = 3.13), then moldova (dist = .4.25), then armenia (5.57)
# note that indicator-level distance does not reflect implicit weighting based on indicator count per concept/sub-obj
# but it nonetheless captures similarity better, which could be useful info for regional programming/planning
fmir %>% filter(country %in% c("Ukraine", "Armenia", "Moldova", "Georgia"), 
                year == 2020,
                obj_num == current_obj_num) %>%
        distinct(country, year, obj_num, obj_avg)
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        filter(country == "Ukraine") %>% select(Georgia, Armenia, Moldova)
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        filter(country == "Armenia") %>% select(Georgia, Ukraine, Moldova)
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        filter(country == "Georgia") %>% select(Ukraine, Armenia, Moldova)
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        filter(country == "Moldova") %>% select(Ukraine, Armenia, Georgia)


fmir %>% filter(country %in% c("Ukraine", "Armenia", "Moldova", "Georgia"), 
                year == 2020,
                obj_num == current_obj_num) %>%
        distinct(country, year, obj_num, obj_avg)
fmir_obj_1_indicator_distance_matrix_tbl %>% select(country, Poland) %>% arrange(Poland)
fmir_obj_1_indicator_distance_matrix_tbl %>% select(country, Albania) %>% arrange(Albania)
fmir_obj_1_indicator_distance_matrix_tbl %>% select(country, Hungary) %>% arrange(Hungary)
fmir_obj_1_indicator_distance_matrix_tbl %>% select(country, `N. Macedonia`) %>% arrange(`N. Macedonia`)

# check distance manually
fmir %>% filter(year == 2020, obj_num == current_obj_num, mcp_grouping != "EU-15") %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        column_to_rownames(var = "country") %>%
        mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1])) %>%
        rownames_to_column(var = "country") %>% 
        filter(country %in% c("Hungary", "N. Macedonia")) %>%
        pivot_longer(cols = -country, names_to = "var", values_to = "values") %>%
        arrange(var, country) %>%
        group_by(var) %>%
        mutate(diff = values - lag(values, n = 1)) %>%
        ungroup() %>% 
        select(var, diff) %>% filter(!is.na(diff)) %>%
        mutate(diff_squared = diff^2,
               sum_diff_squared = sum(diff_squared),
               sqrt_sum_diff_squared = sqrt(sum_diff_squared),
               distance = sqrt_sum_diff_squared)

fmir_obj_1_indicator_distance_matrix_tbl %>% select(country, Poland) %>% arrange(Poland)
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        filter(country == "Poland") %>% select(Georgia, Armenia, Moldova)
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        filter(country == "Hungary") %>% select(Georgia, Ukraine, Moldova)
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        filter(country == "Georgia") %>% select(Ukraine, Armenia, Moldova)
fmir_obj_1_indicator_distance_matrix_tbl %>% 
        filter(country == "Moldova") %>% select(Ukraine, Armenia, Georgia)


# visualize distance
viridis_pal()(3) %>% show_col()
fviz_dist(dist.obj = fmir_obj_1_indicator_distance_matrix, 
          gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_1_indicator_distance_matrix_tile_chart ####


# get country_order from fviz_dist_plot
fmir_obj_1_fviz_dist_plot <- fviz_dist(dist.obj = fmir_obj_1_indicator_distance_matrix, 
                                       gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))
fmir_obj_1_fviz_dist_plot
fmir_obj_1_fviz_dist_plot %>% attributes()
fmir_obj_1_fviz_dist_plot$data %>% as_tibble()

fmir_obj_1_country_order_tbl <- fmir_obj_1_fviz_dist_plot$data %>% as_tibble() %>% slice(1:(fmir_obj_1_indicator %>% nrow())) %>% 
        mutate(country = str_sub(string = Var1, start = 1, end = -2),
               country_order = row_number()) %>% 
        select(country, country_order)
fmir_obj_1_country_order_tbl


#/////////////////////////////


# create chart
fmir_obj_1_indicator_distance_matrix_tile_chart <- fmir_obj_1_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "var", values_to = "values") %>%
        left_join(., fmir_obj_1_country_order_tbl, by = c("var" = "country")) %>%
        rename(var_order = country_order) %>%
        left_join(., fmir_obj_1_country_order_tbl, by = "country") %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = country, .x = country_order, .fun = min), 
                                       y = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Similiarity\n(lower is more similar)", y = NULL, x = NULL) +
        coord_fixed(ratio = 1/1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 13, r = 0, b = 5, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_1_indicator_distance_matrix_tile_chart <- ggdraw(align_legend(fmir_obj_1_indicator_distance_matrix_tile_chart))


# inspect
fmir_obj_1_indicator_distance_matrix_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_1_indicator_distance_matrix_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_1_indicator_distance_matrix_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_1_indicator_similarity_table ####

fmir_obj_1_indicator_similarity_table <- fmir_obj_1_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "counterpart", values_to = "distance") %>%
        filter(country != counterpart) %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping), by = "country") %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping) %>% rename(counterpart_mcp_grouping = mcp_grouping), 
                  by = c("counterpart" = "country")) %>%
        group_by(country) %>%
        mutate(distance_in_ee = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_flag = case_when(distance == min(distance) ~ 1,
                                               TRUE ~ 0),
               max_similarity = case_when(max_similarity_flag == 1 ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_counterpart = case_when(max_similarity_flag == 1 ~ counterpart,
                                          TRUE ~ NA_character_),
               max_similarity_in_ee_flag = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") &
                                                        distance_in_ee == min(distance_in_ee, na.rm = TRUE) ~ 1,
                                                TRUE ~ NA_real_),
               max_similarity_in_ee = case_when(max_similarity_in_ee_flag == 1 ~ distance_in_ee,
                                                TRUE ~ NA_real_),
               max_similarity_in_ee_counterpart = case_when(max_similarity_in_ee_flag == 1 ~ counterpart,
                                                      TRUE ~ NA_character_)) %>%
        fill(max_similarity, .direction = "updown") %>%
        fill(max_similarity_counterpart, .direction = "updown") %>%
        fill(max_similarity_in_ee, .direction = "updown") %>%
        fill(max_similarity_in_ee_counterpart, .direction = "updown") %>%
        ungroup()


#//////////////////


# inspect
fmir_obj_1_indicator_similarity_table
fmir_obj_1_indicator_similarity_table %>% glimpse()
fmir_obj_1_indicator_similarity_table %>% nrow() # 756
fmir_obj_1_indicator_similarity_table %>% ncol() # 12


fmir_obj_1_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping)) %>% 
        print(n = 30)
fmir_obj_1_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping, max_similarity, max_similarity_counterpart)) %>% 
        print(n = 30)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get distinct max_similarity/_in_ee values/counterparts and text_output for table
fmir_obj_1_indicator_similarity_table <- fmir_obj_1_indicator_similarity_table %>% 
        distinct(country, mcp_grouping, max_similarity, max_similarity_counterpart,
                 max_similarity_in_ee, max_similarity_in_ee_counterpart) %>%
        mutate(max_similarity_text_output = str_c("   ", max_similarity_counterpart, " (", round_to_digits(max_similarity, digits = 2), ")"),
               max_similarity_in_ee_text_output = str_c("   ", max_similarity_in_ee_counterpart, " (", 
                                                        round_to_digits(max_similarity_in_ee, digits = 2), ")"),
               mcp_grouping_order = case_when(mcp_grouping == "E&E Balkans" ~ 1,
                                         mcp_grouping == "E&E Eurasia" ~ 2,
                                         mcp_grouping == "CARs" ~ 3,
                                         mcp_grouping == "E&E graduates" ~ 4)) %>%
        arrange(mcp_grouping_order, country) %>%
        select(country, mcp_grouping, 
               max_similarity_text_output, max_similarity, 
               max_similarity_in_ee_text_output, max_similarity_in_ee)


#///////////////////


# inspect
fmir_obj_1_indicator_similarity_table
fmir_obj_1_indicator_similarity_table %>% glimpse()
fmir_obj_1_indicator_similarity_table %>% nrow() # 28
fmir_obj_1_indicator_similarity_table %>% ncol() # 5

fmir_obj_1_indicator_similarity_table %>% print(n = nrow(.))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get table_palette


# use viridis palette
viridis_pal()(10) %>% show_col()


# historic
# note viridis is not used to avoid confusing color scales when compared to dissimilarity matrix, which has much higher
# range of values than the similarity_table
# table_palette <- colorRampPalette(c(brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(4) %>% pull(value), 
#                                      brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(9) %>% pull(value)))


#////////////////////////


# inspect
viridis_pal()(10) %>% show_col()
brewer.pal(n = 9, "Blues") %>% show_col()
table_palette(5) %>% show_col()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get values_normalized_quantile_cutoffs, creating one set of color thresholds to be applied uniformly to all obj in the table
# this just makes it easier to automatically assign color-coded quantiles for the full range of values for E&E countries
# whereas using the actual values values themselves would require checks after updates to make sure the cutoffs are still balanced
# across the full range of values values
# note that this doesn't include miri_avg, but the table does, but this is fine because miri_avg by definition must fall within
# range of max/min values, and since this only gets the color quantiles based on max/min values, it's fine to exclude miri_avg
values_normalized_quantile_cutoffs <- fmir_obj_1_indicator_similarity_table %>%
        pivot_longer(cols = c(max_similarity, max_similarity_in_ee), names_to = "var", values_to = "values") %>%
        mutate(values_min = min(values),
               values_max = max(values),
               values_normalized = (values - values_min) / (values_max - values_min),
               q0_to_q10 = .10 * (values_max - values_min) + values_min,
               q11_to_q20 = .20 * (values_max - values_min) + values_min,
               q21_to_q30 = .30 * (values_max - values_min) + values_min,
               q31_to_q40 = .40 * (values_max - values_min) + values_min,
               q41_to_q50 = .50 * (values_max - values_min) + values_min,
               q51_to_q60 = .60 * (values_max - values_min) + values_min,
               q61_to_q70 = .70 * (values_max - values_min) + values_min,
               q71_to_q80 = .80 * (values_max - values_min) + values_min,
               q81_to_q90 = .90 * (values_max - values_min) + values_min,
               q91_to_q100 = 1 * (values_max - values_min) + values_min) %>%
        select(-c(max_similarity_text_output, max_similarity_in_ee_text_output))


#/////////////////////////


# inspect
values_normalized_quantile_cutoffs %>% data.frame()
# values_normalized = (values - values_min) / (values_max - values_min)

# check q0-q10 threshold
(2.68 - 2.31) / (6.04 - 2.31)

# check q11-q20 threshold
(3.06 - 2.31) / (6.04 - 2.31)

# check q31-q40 threshold
(3.80 - 2.31) / (6.04 - 2.31)

# check q41-q50 threshold
(4.18 - 2.31) / (6.04 - 2.31)

# check q51-q60 threshold
(4.55 - 2.31) / (6.04 - 2.31)

# check q81-q90 threshold
(5.66 - 2.31) / (6.04 - 2.31)

# check q91-q100 threshold
(6.04 - 2.31) / (6.04 - 2.31)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# create table_palette_tbl
# manually apply conditional formatting in excel to highlight and change font color
# =IF(AND(B4>5, B4<10),TRUE,FALSE)
table_palette_tbl <- values_normalized_quantile_cutoffs %>% select(starts_with("q")) %>% distinct() %>%
        pivot_longer(cols = everything(), names_to = "quantile", values_to = "values_threshold") %>%
        left_join(., tibble(hex = viridis_pal()(10), 
                            quantile = c("q0_to_q10", "q11_to_q20", "q21_to_q30", "q31_to_q40", "q41_to_q50",
                                         "q51_to_q60", "q61_to_q70", "q71_to_q80", "q81_to_q90", "q91_to_q100")) %>%
                          mutate(rgb = map(.x = hex, .f = ~ col2rgb(col = .x) %>% str_c(., collapse = "."))) %>%
                          unnest(rgb),
                  by = "quantile")
table_palette_tbl 


#//////////////////////////////////////////////////////////////////////////////////////////////////


# write to file

# note csv will be copied into xlsx shell and formatted
# old: blue header is rgb 31/73/125
# current: gray header is rgb 128/128/128
# grey grid lines are rgb 217/217/217, medium weight
# country rows on left side gets a rgb 217/217/217 heavy border
# grey top/bottom row above/under table is rgb 217/217/217, height = 15
# font is Calibri, size 15, bolded, centered horizontally/vertically in cells; 
# grey font color = rgb 51/51/51, white for column headers
# flag column is width = 15
# country name column is width = 32
# data columns are width = 28
# column header row height = 95
# data row height = 35
# data cells formatted as number with two decimal places, to avoid truncation (e.g. 0.70 converted to 0.7)

# remember to fully spell out macedonia and bih
# output is copied into powerpoint using paste special -> Picture (Enhanced Metafile), aligned center
# flags are copy/pasted in from state dept website, w locked aspect ratio, and resized to 70%

# save
# fmir_obj_1_indicator_similarity_table %>% write_csv(file = "output/charts/fmir_obj_1_indicator_similarity_table.csv")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# get cluster tendency to see if data is likely to have clusters that are more meaningful/tight than uniform random sample
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# note that the hopkins stat uses .5 as a threshold, with higher values indicating there is a clustering tendency; 
# a hopkins stat over .75 is 90% confidence that there is clustering tendency that is real, and not random
cluster_tendency <- get_clust_tendency(fmir_obj_1_indicator, n = nrow(fmir_obj_1_indicator) - 1, graph = FALSE)


#///////////////////


# inspect
cluster_tendency
cluster_tendency$hopkins_stat


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters
fmir_obj_1_indicator_hclust <- hclust(d = fmir_obj_1_indicator_distance_matrix, method = "complete" )


#/////////////////////


# inspect
fmir_obj_1_indicator_hclust
fmir_obj_1_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_1_indicator_hclust$merge
fmir_obj_1_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_1_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters with agnes, which allows enhanced fviz visualizations
fmir_obj_1_indicator_agnes <- agnes(x = fmir_obj_1_indicator, method = "complete")


#//////////////////


# inspect
fmir_obj_1_indicator_agnes
fmir_obj_1_indicator_agnes %>% attributes()

# get agglomerative coefficient from agnes object 
# "Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 suggest a 
# more balanced clustering structure such as the complete linkage and Wards method dendrograms in Figure 21.3. 
# Values closer to 0 suggest less well-formed clusters such as the single linkage dendrogram in Figure 21.3. 
# However, the AC tends to become larger as n increases, so it should not be used to compare across data sets of very different sizes."
fmir_obj_1_indicator_agnes$ac

# plot
pltree(fmir_obj_1_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 



#//////////////////////////////////////////////////////////////////////////////////////////////////


# use agnes to look for clustering method with best agglomerative coefficient 
method_list <- c( "average", "single", "complete", "ward")
names(method_list) <- c( "average", "single", "complete", "ward")
method_list

# function to compute coefficient
cluster_methods <- function(current_method) {
        print(current_method)
        agnes(fmir_obj_1_indicator_distance_matrix, method = current_method)$ac
}

# run cluster_methods
# note that wards method has the highest/best agglomerative coefficient
map(.x = method_list, .f = ~ cluster_methods(current_method = .x))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with agnes, which had highest (best) agglomerative coefficient
fmir_obj_1_indicator_agnes <- agnes(fmir_obj_1_indicator, method = "ward")
fmir_obj_1_indicator_agnes


#///////////////////


# visualize clustering 
pltree(fmir_obj_1_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
fviz_dend(fmir_obj_1_indicator_agnes, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

# inspect optimal number of clusters
# note the gap statistic is recently developed by tibshirani et al
fviz_nbclust(fmir_obj_1_indicator, FUN = hcut, method = "wss")
fviz_nbclust(fmir_obj_1_indicator, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(fmir_obj_1_indicator, FUN = hcut, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# set optimal k
fmir_obj_1_indicator_k <- 3


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with hclust
fmir_obj_1_indicator_hclust <- hclust(d = fmir_obj_1_indicator_distance_matrix, method = "ward.D2")


#/////////////////////


# inspect
fmir_obj_1_indicator_hclust
fmir_obj_1_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_1_indicator_hclust$merge
fmir_obj_1_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_1_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut tree 

# build and cut agnes tree with hcut for enhanced fviz
fmir_obj_1_indicator_agnes_cut <- hcut(fmir_obj_1_indicator, k = fmir_obj_1_indicator_k, hc_method = "ward.D2", stand = TRUE)
fmir_obj_1_indicator_agnes_cut 

# view PCA for clusters, which requires an hcut object
fviz_cluster(object = fmir_obj_1_indicator_agnes_cut)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut hclust tree with cutree
fmir_obj_1_indicator_hclust_cut <- cutree(fmir_obj_1_indicator_hclust, k = fmir_obj_1_indicator_k)
fmir_obj_1_indicator_hclust_cut

# Number of members in each cluster
table(fmir_obj_1_indicator_hclust_cut)

# can add cutree output to df
fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value) %>%
        left_join(fmir_obj_1_indicator %>% rownames_to_column(var = "country"), ., by = "country")

# plot dendrogram with border around 4 clusters
plot(fmir_obj_1_indicator_hclust, cex = 0.6, hang = -1)
rect.hclust(fmir_obj_1_indicator_hclust, k = fmir_obj_1_indicator_k, border = 2:5)

# plot dendrogram w labels
dendrogram_obj <- as.dendrogram(fmir_obj_1_indicator_hclust)
dendrogram_obj %>% unclass() %>% str()
dendrogram_obj %>% unclass() %>% class()

subdendrogram_list <- get_subdendrograms(dendrogram_obj, k = fmir_obj_1_indicator_k)
subdendrogram_list
dendrogram_labels <- map(.x = subdendrogram_list, .f = ~ labels(.x)) %>% 
        enframe() %>% unnest(cols = value) %>% 
        rename(subdendrogram_plot_order = name, country = value) %>%
        left_join(fmir_obj_1_indicator %>% mutate(cluster = fmir_obj_1_indicator_hclust_cut) %>%
                          rownames_to_column(var = "country"),
                  ., by = "country") %>% 
        distinct(cluster, subdendrogram_plot_order) %>% 
        arrange(subdendrogram_plot_order) %>%
        mutate(cluster = str_c("\n\n\n", cluster)) %>%
        pull(cluster) %>% unname()
dendrogram_labels        

plot(fmir_obj_1_indicator_hclust, hang = -1)
rect.dendrogram(as.dendrogram(fmir_obj_1_indicator_hclust), k = fmir_obj_1_indicator_k, text = dendrogram_labels)

# visualize clusters along 2 primary PCA dimensions
fviz_cluster(list(data = fmir_obj_1_indicator, cluster = fmir_obj_1_indicator_hclust_cut))


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_1_indicator_dendrogram ####

# get fmir_obj_1_indicator_dendrogram_data
# note that dendro_data_k() takes the uncut hclust as first argument
fmir_obj_1_indicator_dendrogram_data <- dendro_data_k(hc = fmir_obj_1_indicator_hclust, 
                                                      k = fmir_obj_1_indicator_k)

# inspect
fmir_obj_1_indicator_dendrogram_data
fmir_obj_1_indicator_dendrogram_data %>% attributes()
fmir_obj_1_indicator_dendrogram_data$segments %>% head()
fmir_obj_1_indicator_dendrogram_data$labels %>% head()
fmir_obj_1_indicator_dendrogram_data$leaf_labels %>% head()
fmir_obj_1_indicator_dendrogram_data$class %>% head()


#///////////////////


# get cluster_colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
cluster_colors <- c("#333333", 
                    color_palette %>% slice(2) %>% pull(hex), 
                    color_palette %>% slice(7) %>% pull(hex),
                    color_palette %>% slice(10) %>% pull(hex))
cluster_colors


# plot
fmir_obj_1_indicator_dendrogram <- plot_ggdendro(fmir_obj_1_indicator_dendrogram_data,
                                                 direction = "tb",
                                                 scale.color = cluster_colors,
                                                 label.size = 2.5,
                                                 branch.size = 0.5,
                                                 expand.y = 0) +
        labs(y = "Similarity\n(lower groupings are more similar)", x = "") +
        coord_fixed(ratio = 1 / 1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 5, 0, 0), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(size = .5),
                axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "bold", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 0, hjust = .5),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                # axis.title.x = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 5, r = 0, b = 10, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 0, r = 15, b = 0, l = 5)),
                # axis.title.y = element_blank(),
                # axis.text.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 9, family = "Calibri", face = "plain", color = "#333333"),
                legend.text = element_text(size = 9, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# guides(color = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333", keywidth = 4),
#        linetype = guide_legend(keywidth = 4))

# inspect
fmir_obj_1_indicator_dendrogram


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_1_indicator_dendrogram)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_1_indicator_dendrogram.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_1_country_indicator_tile_chart ####

# get chart_data
# note this is 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        left_join(., fmir_obj_1_country_order_tbl, by = "country") %>%
        relocate(country_order, .after = cluster) %>% arrange(country_order) %>%
        pivot_wider(id_cols = c(country, cluster, country_order), 
                    names_from = indicator_name, values_from = indicator_standardized_values) %>%
        left_join(., fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          distinct(country, obj_num, sub_obj_num, sub_obj_avg, obj_avg) %>%
                          mutate(sub_obj_num_label = case_when(sub_obj_num == "sub_obj_1_1" ~ 
                                                                       "Sub-obj. 1.1: Checks & balances and rule of law",
                                                               sub_obj_num == "sub_obj_1_2" ~ "Sub-obj. 1.2: Civil society",
                                                               sub_obj_num == "sub_obj_1_3" ~ 
                                                                     "Sub-obj. 1.3: Reslience to electoral/political interference")) %>%
                          select(country, obj_num, obj_avg, sub_obj_num_label, sub_obj_avg) %>%
                          pivot_wider(id_cols = c(country, obj_num, obj_avg), 
                                      names_from = sub_obj_num_label, values_from = sub_obj_avg) %>%
                          rename(!!sym("Obj. 1: Democratic") := obj_avg) %>% select(-obj_num), 
                  by = "country") %>%
        pivot_longer(cols = -c(country, cluster, country_order), names_to = "var", values_to = "values") %>%
        mutate(color = case_when(cluster == 1 ~ color_palette %>% slice(2) %>% pull(hex),
                                 cluster == 2 ~ color_palette %>% slice(7) %>% pull(hex),
                                 cluster == 3 ~ color_palette %>% slice(10) %>% pull(hex)),
               country_label = str_c("Group ", cluster, ": ", country),
               var_order = case_when(var == "Obj. 1: Democratic" ~ 1,
                                     str_detect(string = var, pattern = "^Sub-obj.") ~ 2,
                                     TRUE ~ 3))

# inspect
chart_data


#/////////////////////////////


# create chart
fmir_obj_1_country_indicator_tile_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       y = fct_reorder(.f = country_label, .x = country_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Score\n(higher is better)\n", y = NULL, x = "Indicator") +
        coord_fixed(ratio = 1 / 3, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 5, 
                                           color = chart_data %>% distinct(country, color) %>% pull(color), 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 5, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 5, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_1_country_indicator_tile_chart <- ggdraw(align_legend(fmir_obj_1_country_indicator_tile_chart))


# inspect
fmir_obj_1_country_indicator_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_1_country_indicator_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_1_country_indicator_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_1_indicator_cluster_boxplot_chart ####

# get chart_data
# note that this uses the original ISV values from the index, not the scaled ISV used in fmir_obj_1_indicator for clustering; 
# original/unscaled ISV avg is better for plot because values are then interpretable in the framework of the index (e.g. ISV btw 0-1)
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        
        
        # add cluster avg of sub_obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -c(country, cluster), names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster),
               color_bin = as.character(cluster_name),
               color = case_when(color_bin == "Group 1" ~ color_palette %>% slice(2) %>% pull(hex),
                                 color_bin == "Group 2" ~ color_palette %>% slice(7) %>% pull(hex),
                                 color_bin == "Group 3" ~ color_palette %>% slice(10) %>% pull(hex)))

# inspect
chart_data
chart_data %>% glimpse()

# create color_list for to pass to scale_color_manual
chart_data_color_list <- chart_data %>% count(color_bin, color) %>% pull(color)
names(chart_data_color_list) <- chart_data %>% count(color_bin, color) %>% pull(color_bin)
chart_data_color_list


#//////////////////////////


# create chart
fmir_obj_1_indicator_cluster_boxplot_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = values, fill = color_bin)) + 
        geom_boxplot(width = .5, lwd = .01, outlier.size = .5) + facet_wrap(facets = vars(cluster_name)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = .2), limits = c(0, 1), expand = c(0, 0)) +
        scale_fill_manual(values = chart_data_color_list) + 
        # scale_color_manual(values = chart_data_color_list) + 
        coord_fixed(ratio = 10/.6, clip = "off") +
        labs(x = "Indicator", y = "Score (higher is better)", fill = NULL, color = NULL) +
        theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                # strip.background = element_blank(),
                strip.background = element_rect(fill = "#DDDDDD", color = NA),
                strip.text = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333"),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = -15, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 5, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = 1,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 


# inspect
fmir_obj_1_indicator_cluster_boxplot_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_1_indicator_cluster_boxplot_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_1_indicator_cluster_boxplot_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_1_indicator_cluster_tile_chart ####

# get chart_data 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        group_by(cluster) %>%
        mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
        ungroup() %>%
        select(-country) %>% distinct() %>%
        
        
        # add cluster avg of sub_obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -cluster, names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster))

# inspect
chart_data


#///////////////


# get obj/sub_obj segment coordinates

# note that in the entire fmir dataset, the maximum sub_obj per obj is 3, so will have 3 segment_tbls,
# and for those obj with less than 3 sub_obj, the segment_tbl will just have NA coordinates and no output in ggplot
fmir %>% distinct(obj_num, sub_obj_num) %>% 
        add_count(obj_num, name = "sub_obj_count") %>% 
        distinct(obj_num, sub_obj_count) %>%
        arrange(desc(sub_obj_count)) 

# segment_1_tbl
# segment_1_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(obj_num, sub_obj_num) %>%
#         mutate(nrow = nrow(.),
#                x = nrow + 1.5,
#                xend = nrow + 1.5,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5) %>%
#         select(x, xend, y, yend) %>%
#         slice(1)
# segment_1_tbl

# # segment_2_tbl
# segment_2_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                 cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(1)
# segment_2_tbl
#                
# # segment_3_tbl
# segment_3_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(2)
# segment_3_tbl


#///////////////////////


# plot

# get order_obj_and_sub_obj_avg_first()
order_obj_and_sub_obj_avg_first <- function(var){
        
        tibble(var = var) %>% 
                left_join(., fmir %>% distinct(indicator_name, concept) %>% 
                                  add_group_index(group_vars = concept, group_name = "concept_index"), 
                          by = c("var" = "indicator_name")) %>%
                mutate(order = case_when(str_detect(string = var, pattern = "^obj_") ~ 1,
                                         str_detect(string = var, pattern = "^sub_obj_[0-9]_[0-9]$") ~ 2,
                                         TRUE ~ 3)) %>%
                arrange(order, concept_index, var) %>% pull(var)
}

# create chart
fmir_obj_1_indicator_cluster_tile_chart <- chart_data %>% 
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = factor(cluster_name), fill = values, label = round_to_digits(values, digits = 2))) + 
        geom_tile(color = "#DDDDDD") + 
        geom_text(color = "#ffffff", size = 1.75, family = "Calibri") +
        scale_y_discrete(expand = c(0, 0)) +
        # geom_segment(x = segment_1_tbl$x, xend = segment_1_tbl$xend, 
        #              y = segment_1_tbl$y, yend = segment_1_tbl$yend, color = "#ffffff", size = 3) +
        # geom_segment(x = segment_2_tbl$x, xend = segment_2_tbl$xend,
        #              y = segment_2_tbl$y, yend = segment_2_tbl$yend, color = "#CBCBCB", size = 3) +
        # geom_segment(x = segment_3_tbl$x, xend = segment_3_tbl$xend,
        #              y = segment_3_tbl$y, yend = segment_3_tbl$yend, color = "#CBCBCB", size = 3) +
        labs(x = "Indicator", y = NULL, fill = "Avg. score\n(higher is better)\n") +
        coord_fixed(ratio = 2 / 1, clip = "off") +
        theme_bw() +
        # theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD"),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "right",
                legend.key.size = unit(4, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 5, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# fmir_obj_1_indicator_cluster_tile_chart <- ggdraw(align_legend(fmir_obj_1_indicator_cluster_tile_chart))

# inspect
fmir_obj_1_indicator_cluster_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_1_indicator_cluster_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_1_indicator_cluster_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_1_indicator_tree ####

fmir_obj_1_indicator_unscaled
fmir_obj_1_indicator_unscaled %>% glimpse()

tree_data <- fmir_obj_1_indicator_unscaled %>% 
        
        # add clusters
        left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        select(-country)


#//////////////////////////


# inspect
tree_data
tree_data %>% glimpse()
tree_data %>% nrow()
tree_data %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////////


# get tree
fmir_obj_1_indicator_tree <- rpart(formula = cluster ~ .,
                                   data = tree_data, method = "class",
                                   minsplit = 2, minbucket = 1, cp = -1)


#/////////////////////


# inspect
fmir_obj_1_indicator_tree
rpart.plot(fmir_obj_1_indicator_tree, tweak = 1.3, extra = 104, nn = TRUE)

# inspect specific indicators
fmir %>% filter(year == 2020, mcp_grouping != "EU-15") %>%
        filter(indicator_name == "sub_obj_1_1_fh_local_democratic_gov") %>%
        select(country, indicator_standardized_values) %>% arrange(indicator_standardized_values) %>% 
        left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        print(n = nrow(.))


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_2_indicator_clusters ####

# set current_obj_num
current_obj_num <- "obj_2"

# get fmir_obj_2_indicator with scaled version of indicator_standardized_values
fmir_obj_2_indicator <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        distinct(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        column_to_rownames(var = "country") %>%
        mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1]))


#//////////////////////


# inspect
fmir_obj_2_indicator
fmir_obj_2_indicator %>% glimpse()
fmir_obj_2_indicator %>% nrow()
fmir_obj_2_indicator %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get fmir_obj_2_indicator_distance_matrix 
fmir_obj_2_indicator_distance_matrix <- dist(fmir_obj_2_indicator, method = "euclidean")
# fmir_obj_2_indicator_distance_matrix <- get_dist(fmir_obj_2_indicator, method = "euclidean")


#//////////////


# inspect
fmir_obj_2_indicator_distance_matrix
fmir_obj_2_indicator_distance_matrix %>% str()
fmir_obj_2_indicator_distance_matrix %>% attributes()

# inspect fmir_obj_2_indicator_distance_matrix_tbl
fmir_obj_2_indicator_distance_matrix_tbl <- as.matrix(fmir_obj_2_indicator_distance_matrix) %>% as_tibble() %>% 
        mutate(country = row.names(as.matrix(fmir_obj_2_indicator_distance_matrix))) %>%
        relocate(country, .before = everything())

fmir_obj_2_indicator_distance_matrix_tbl
fmir_obj_2_indicator_distance_matrix_tbl %>% skim() %>% as_tibble() %>% 
        select(skim_variable, numeric.mean, starts_with(match = "numeric.p")) %>%
        # arrange(numeric.mean) %>% 
        arrange(numeric.p25) %>%
        print(n = nrow(.))
fmir_obj_2_indicator_distance_matrix_tbl %>% 
        # select(country, Serbia) %>% arrange(Serbia) %>% print(n = nrow(.))
        # select(country, BiH) %>% arrange(BiH) %>% print(n = nrow(.))
        # select(country, Kosovo) %>% arrange(Kosovo) %>% print(n = nrow(.))
        # select(country, Armenia) %>% arrange(Armenia) %>% print(n = nrow(.))
        select(country, Ukraine) %>% arrange(Ukraine) %>% print(n = nrow(.))
# filter(country == "Serbia") %>% select(ends_with("stan"))

# compare indicator_standardized_values
fmir %>% filter(obj_num == "obj_3", country %in% c("Poland", "Kazakhstan"), year == 2020) %>% 
        select(country, year, indicator_name, values, values_z_std, indicator_standardized_values) %>%
        arrange(country) %>% 
        # unite(col = var, indicator_name, country) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = indicator_name, names_from = country, values_from = indicator_standardized_values)



# visualize distance
viridis_pal()(3) %>% show_col()
fviz_dist(dist.obj = fmir_obj_2_indicator_distance_matrix, gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_2_indicator_distance_matrix_tile_chart ####


# get fmir_obj_2_country_order from fmir_obj_2_fviz_dist_plot
fmir_obj_2_fviz_dist_plot <- fviz_dist(dist.obj = fmir_obj_2_indicator_distance_matrix, 
                                       gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))
fmir_obj_2_fviz_dist_plot
fmir_obj_2_fviz_dist_plot %>% attributes()
fmir_obj_2_fviz_dist_plot$data %>% as_tibble()

fmir_obj_2_country_order_tbl <- fmir_obj_2_fviz_dist_plot$data %>% as_tibble() %>% slice(1:(fmir_obj_2_indicator %>% nrow())) %>% 
        mutate(country = str_sub(string = Var1, start = 1, end = -2),
               country_order = row_number()) %>% 
        select(country, country_order)
fmir_obj_2_country_order_tbl


#/////////////////////////////


# create chart
fmir_obj_2_indicator_distance_matrix_tile_chart <- fmir_obj_2_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "var", values_to = "values") %>%
        left_join(., fmir_obj_2_country_order_tbl, by = c("var" = "country")) %>%
        rename(var_order = country_order) %>%
        left_join(., fmir_obj_2_country_order_tbl, by = "country") %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = country, .x = country_order, .fun = min), 
                                       y = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Similiarity\n(lower is more similar)", y = NULL, x = NULL) +
        coord_fixed(ratio = 1/1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 13, r = 0, b = 5, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_2_indicator_distance_matrix_tile_chart <- ggdraw(align_legend(fmir_obj_2_indicator_distance_matrix_tile_chart))


# inspect
fmir_obj_2_indicator_distance_matrix_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_2_indicator_distance_matrix_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_2_indicator_distance_matrix_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_2_indicator_similarity_table ####

fmir_obj_2_indicator_similarity_table <- fmir_obj_2_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "counterpart", values_to = "distance") %>%
        filter(country != counterpart) %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping), by = "country") %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping) %>% rename(counterpart_mcp_grouping = mcp_grouping), 
                  by = c("counterpart" = "country")) %>%
        group_by(country) %>%
        mutate(distance_in_ee = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_flag = case_when(distance == min(distance) ~ 1,
                                               TRUE ~ 0),
               max_similarity = case_when(max_similarity_flag == 1 ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_counterpart = case_when(max_similarity_flag == 1 ~ counterpart,
                                                      TRUE ~ NA_character_),
               max_similarity_in_ee_flag = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") &
                                                             distance_in_ee == min(distance_in_ee, na.rm = TRUE) ~ 1,
                                                     TRUE ~ NA_real_),
               max_similarity_in_ee = case_when(max_similarity_in_ee_flag == 1 ~ distance_in_ee,
                                                TRUE ~ NA_real_),
               max_similarity_in_ee_counterpart = case_when(max_similarity_in_ee_flag == 1 ~ counterpart,
                                                            TRUE ~ NA_character_)) %>%
        fill(max_similarity, .direction = "updown") %>%
        fill(max_similarity_counterpart, .direction = "updown") %>%
        fill(max_similarity_in_ee, .direction = "updown") %>%
        fill(max_similarity_in_ee_counterpart, .direction = "updown") %>%
        ungroup()


#//////////////////


# inspect
fmir_obj_2_indicator_similarity_table
fmir_obj_2_indicator_similarity_table %>% glimpse()
fmir_obj_2_indicator_similarity_table %>% nrow() # 756
fmir_obj_2_indicator_similarity_table %>% ncol() # 12


fmir_obj_2_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping)) %>% 
        print(n = 30)
fmir_obj_2_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping, max_similarity, max_similarity_counterpart)) %>% 
        print(n = 30)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get distinct max_similarity/_in_ee values/counterparts and text_output for table
fmir_obj_2_indicator_similarity_table <- fmir_obj_2_indicator_similarity_table %>% 
        distinct(country, mcp_grouping, max_similarity, max_similarity_counterpart,
                 max_similarity_in_ee, max_similarity_in_ee_counterpart) %>%
        mutate(max_similarity_text_output = str_c("   ", max_similarity_counterpart, " (", round_to_digits(max_similarity, digits = 2), ")"),
               max_similarity_in_ee_text_output = str_c("   ", max_similarity_in_ee_counterpart, " (", 
                                                        round_to_digits(max_similarity_in_ee, digits = 2), ")"),
               mcp_grouping_order = case_when(mcp_grouping == "E&E Balkans" ~ 1,
                                              mcp_grouping == "E&E Eurasia" ~ 2,
                                              mcp_grouping == "CARs" ~ 3,
                                              mcp_grouping == "E&E graduates" ~ 4)) %>%
        arrange(mcp_grouping_order, country) %>%
        select(country, mcp_grouping, 
               max_similarity_text_output, max_similarity, 
               max_similarity_in_ee_text_output, max_similarity_in_ee)


#///////////////////


# inspect
fmir_obj_2_indicator_similarity_table
fmir_obj_2_indicator_similarity_table %>% glimpse()
fmir_obj_2_indicator_similarity_table %>% nrow() # 28
fmir_obj_2_indicator_similarity_table %>% ncol() # 5

fmir_obj_2_indicator_similarity_table %>% print(n = nrow(.))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get table_palette


# use viridis palette
viridis_pal()(10) %>% show_col()


# historic
# note viridis is not used to avoid confusing color scales when compared to dissimilarity matrix, which has much higher
# range of values than the similarity_table
# table_palette <- colorRampPalette(c(brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(4) %>% pull(value), 
#                                      brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(9) %>% pull(value)))


#////////////////////////


# inspect
viridis_pal()(10) %>% show_col()
brewer.pal(n = 9, "Blues") %>% show_col()
table_palette(5) %>% show_col()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get values_normalized_quantile_cutoffs, creating one set of color thresholds to be applied uniformly to all obj in the table
# this just makes it easier to automatically assign color-coded quantiles for the full range of values for E&E countries
# whereas using the actual values values themselves would require checks after updates to make sure the cutoffs are still balanced
# across the full range of values values
# note that this doesn't include miri_avg, but the table does, but this is fine because miri_avg by definition must fall within
# range of max/min values, and since this only gets the color quantiles based on max/min values, it's fine to exclude miri_avg
values_normalized_quantile_cutoffs <- fmir_obj_2_indicator_similarity_table %>%
        pivot_longer(cols = c(max_similarity, max_similarity_in_ee), names_to = "var", values_to = "values") %>%
        mutate(values_min = min(values),
               values_max = max(values),
               values_normalized = (values - values_min) / (values_max - values_min),
               q0_to_q10 = .10 * (values_max - values_min) + values_min,
               q11_to_q20 = .20 * (values_max - values_min) + values_min,
               q21_to_q30 = .30 * (values_max - values_min) + values_min,
               q31_to_q40 = .40 * (values_max - values_min) + values_min,
               q41_to_q50 = .50 * (values_max - values_min) + values_min,
               q51_to_q60 = .60 * (values_max - values_min) + values_min,
               q61_to_q70 = .70 * (values_max - values_min) + values_min,
               q71_to_q80 = .80 * (values_max - values_min) + values_min,
               q81_to_q90 = .90 * (values_max - values_min) + values_min,
               q91_to_q100 = 1 * (values_max - values_min) + values_min) %>%
        select(-c(max_similarity_text_output, max_similarity_in_ee_text_output))


#/////////////////////////


# inspect
values_normalized_quantile_cutoffs %>% data.frame()
# values_normalized = (values - values_min) / (values_max - values_min)

# check q0-q10 threshold
(2.25 - 1.81) / (6.18 - 1.81)

# check q81-q90 threshold
(5.74 - 1.81) / (6.18 - 1.81)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# create table_palette_tbl
# manually apply conditional formatting in excel to highlight and change font color
# =IF(AND(B4>5, B4<10),TRUE,FALSE)
table_palette_tbl <- values_normalized_quantile_cutoffs %>% select(starts_with("q")) %>% distinct() %>%
        pivot_longer(cols = everything(), names_to = "quantile", values_to = "values_threshold") %>%
        left_join(., tibble(hex = viridis_pal()(10), 
                            quantile = c("q0_to_q10", "q11_to_q20", "q21_to_q30", "q31_to_q40", "q41_to_q50",
                                         "q51_to_q60", "q61_to_q70", "q71_to_q80", "q81_to_q90", "q91_to_q100")) %>%
                          mutate(rgb = map(.x = hex, .f = ~ col2rgb(col = .x) %>% str_c(., collapse = "."))) %>%
                          unnest(rgb),
                  by = "quantile")

# inspect
table_palette_tbl 


#//////////////////////////////////////////////////////////////////////////////////////////////////


# write to file

# note csv will be copied into xlsx shell and formatted
# old: blue header is rgb 31/73/125
# current: gray header is rgb 128/128/128
# grey grid lines are rgb 217/217/217, medium weight
# country rows on left side gets a rgb 217/217/217 heavy border
# grey top/bottom row above/under table is rgb 217/217/217, height = 15
# font is Calibri, size 15, bolded, centered horizontally/vertically in cells; 
# grey font color = rgb 51/51/51, white for column headers
# flag column is width = 15
# country name column is width = 32
# data columns are width = 28
# column header row height = 95
# data row height = 35
# data cells formatted as number with two decimal places, to avoid truncation (e.g. 0.70 converted to 0.7)

# remember to fully spell out macedonia and bih
# output is copied into powerpoint using paste special -> Picture (Enhanced Metafile), aligned center
# flags are copy/pasted in from state dept website, w locked aspect ratio, and resized to 70%

# save
# fmir_obj_2_indicator_similarity_table %>% write_csv(file = "output/charts/fmir_obj_2_indicator_similarity_table.csv")


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# get cluster tendency to see if data is likely to have clusters that are more meaningful/tight than uniform random sample
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# note that the hopkins stat uses .5 as a threshold, with higher values indicating there is a clustering tendency; 
# a hopkins stat over .75 is 90% confidence that there is clustering tendency that is real, and not random
cluster_tendency <- get_clust_tendency(fmir_obj_2_indicator, n = nrow(fmir_obj_2_indicator) - 1, graph = FALSE)


#///////////////////


# inspect
cluster_tendency
cluster_tendency$hopkins_stat


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters
fmir_obj_2_indicator_hclust <- hclust(d = fmir_obj_2_indicator_distance_matrix, method = "complete" )


#/////////////////////


# inspect
fmir_obj_2_indicator_hclust
fmir_obj_2_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_2_indicator_hclust$merge
fmir_obj_2_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_2_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters with agnes, which allows enhanced fviz visualizations
fmir_obj_2_indicator_agnes <- agnes(x = fmir_obj_2_indicator, method = "complete")


#//////////////////


# inspect
fmir_obj_2_indicator_agnes
fmir_obj_2_indicator_agnes %>% attributes()

# get agglomerative coefficient from agnes object 
# "Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 suggest a 
# more balanced clustering structure such as the complete linkage and Wards method dendrograms in Figure 21.3. 
# Values closer to 0 suggest less well-formed clusters such as the single linkage dendrogram in Figure 21.3. 
# However, the AC tends to become larger as n increases, so it should not be used to compare across data sets of very different sizes."
fmir_obj_2_indicator_agnes$ac

# plot
pltree(fmir_obj_2_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 



#//////////////////////////////////////////////////////////////////////////////////////////////////


# use agnes to look for clustering method with best agglomerative coefficient 
method_list <- c( "average", "single", "complete", "ward")
names(method_list) <- c( "average", "single", "complete", "ward")
method_list

# function to compute coefficient
cluster_methods <- function(current_method) {
        print(current_method)
        agnes(fmir_obj_2_indicator_distance_matrix, method = current_method)$ac
}

# run cluster_methods
# note that wards method has the highest/best agglomerative coefficient
map(.x = method_list, .f = ~ cluster_methods(current_method = .x))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with agnes, which had highest (best) agglomerative coefficient
fmir_obj_2_indicator_agnes <- agnes(fmir_obj_2_indicator, method = "ward")
fmir_obj_2_indicator_agnes


#///////////////////


# visualize clustering 
pltree(fmir_obj_2_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
fviz_dend(fmir_obj_2_indicator_agnes, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

# inspect optimal number of clusters
# note the gap statistic is recently developed by tibshirani et al
fviz_nbclust(fmir_obj_2_indicator, FUN = hcut, method = "wss")
fviz_nbclust(fmir_obj_2_indicator, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(fmir_obj_2_indicator, FUN = hcut, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# set optimal k
fmir_obj_2_indicator_k <- 3


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with hclust
fmir_obj_2_indicator_hclust <- hclust(d = fmir_obj_2_indicator_distance_matrix, method = "ward.D2")


#/////////////////////


# inspect
fmir_obj_2_indicator_hclust
fmir_obj_2_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_2_indicator_hclust$merge
fmir_obj_2_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_2_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut tree 

# build and cut agnes tree with hcut for enhanced fviz
fmir_obj_2_indicator_agnes_cut <- hcut(fmir_obj_2_indicator, k = fmir_obj_2_indicator_k, hc_method = "ward.D2", stand = TRUE)
fmir_obj_2_indicator_agnes_cut 

# view PCA for clusters, which requires an hcut object
fviz_cluster(object = fmir_obj_2_indicator_agnes_cut)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut hclust tree with cutree
fmir_obj_2_indicator_hclust_cut <- cutree(fmir_obj_2_indicator_hclust, k = fmir_obj_2_indicator_k)
fmir_obj_2_indicator_hclust_cut

# Number of members in each cluster
table(fmir_obj_2_indicator_hclust_cut)

# can add cutree output to df
fmir_obj_2_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value) %>%
        left_join(fmir_obj_2_indicator %>% rownames_to_column(var = "country"), ., by = "country")

# plot dendrogram with border around 4 clusters
plot(fmir_obj_2_indicator_hclust, cex = 0.6, hang = -1)
rect.hclust(fmir_obj_2_indicator_hclust, k = fmir_obj_2_indicator_k, border = 2:5)

# plot dendrogram w labels
dendrogram_obj <- as.dendrogram(fmir_obj_2_indicator_hclust)
dendrogram_obj %>% unclass() %>% str()
dendrogram_obj %>% unclass() %>% class()

subdendrogram_list <- get_subdendrograms(dendrogram_obj, k = fmir_obj_2_indicator_k)
subdendrogram_list
dendrogram_labels <- map(.x = subdendrogram_list, .f = ~ labels(.x)) %>% 
        enframe() %>% unnest(cols = value) %>% 
        rename(subdendrogram_plot_order = name, country = value) %>%
        left_join(fmir_obj_2_indicator %>% mutate(cluster = fmir_obj_2_indicator_hclust_cut) %>%
                          rownames_to_column(var = "country"),
                  ., by = "country") %>% 
        distinct(cluster, subdendrogram_plot_order) %>% 
        arrange(subdendrogram_plot_order) %>%
        mutate(cluster = str_c("\n\n\n", cluster)) %>%
        pull(cluster) %>% unname()
dendrogram_labels        

plot(fmir_obj_2_indicator_hclust, hang = -1)
rect.dendrogram(as.dendrogram(fmir_obj_2_indicator_hclust), k = fmir_obj_2_indicator_k, text = dendrogram_labels)

# visualize clusters along 2 primary PCA dimensions
fviz_cluster(list(data = fmir_obj_2_indicator, cluster = fmir_obj_2_indicator_hclust_cut))


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_2_indicator_dendrogram ####

# get fmir_obj_2_indicator_dendrogram_data
# note that dendro_data_k() takes the uncut hclust as first argument
fmir_obj_2_indicator_dendrogram_data <- dendro_data_k(hc = fmir_obj_2_indicator_hclust, 
                                                      k = fmir_obj_2_indicator_k)

# inspect
fmir_obj_2_indicator_dendrogram_data
fmir_obj_2_indicator_dendrogram_data %>% attributes()
fmir_obj_2_indicator_dendrogram_data$segments %>% head()
fmir_obj_2_indicator_dendrogram_data$labels %>% head()
fmir_obj_2_indicator_dendrogram_data$leaf_labels %>% head()
fmir_obj_2_indicator_dendrogram_data$class %>% head()


#///////////////////


# get cluster_colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
cluster_colors <- c("#333333", 
                    color_palette %>% slice(2) %>% pull(hex), 
                    color_palette %>% slice(7) %>% pull(hex),
                    color_palette %>% slice(10) %>% pull(hex))
cluster_colors


# plot
fmir_obj_2_indicator_dendrogram <- plot_ggdendro(fmir_obj_2_indicator_dendrogram_data,
                                                 direction = "tb",
                                                 scale.color = cluster_colors,
                                                 label.size = 2.5,
                                                 branch.size = 0.5,
                                                 expand.y = 0) +
        labs(y = "Similarity\n(lower groupings are more similar)", x = "") +
        coord_fixed(ratio = 1 / 1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 5, 0, 0), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(size = .5),
                axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "bold", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 0, hjust = .5),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                # axis.title.x = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 5, r = 0, b = 10, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 0, r = 15, b = 0, l = 5)),
                # axis.title.y = element_blank(),
                # axis.text.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 9, family = "Calibri", face = "plain", color = "#333333"),
                legend.text = element_text(size = 9, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# guides(color = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333", keywidth = 4),
#        linetype = guide_legend(keywidth = 4))

# inspect
fmir_obj_2_indicator_dendrogram


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_2_indicator_dendrogram)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_2_indicator_dendrogram.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_2_country_indicator_tile_chart ####

# get chart_data
# note this is 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        left_join(., fmir_obj_1_country_order_tbl, by = "country") %>%
        relocate(country_order, .after = cluster) %>% arrange(country_order) %>%
        pivot_wider(id_cols = c(country, cluster, country_order), 
                    names_from = indicator_name, values_from = indicator_standardized_values) %>%
        left_join(., fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          distinct(country, obj_num, sub_obj_num, sub_obj_avg, obj_avg) %>%
                          mutate(sub_obj_num_label = case_when(sub_obj_num == "sub_obj_2_1" ~ 
                                                                       "Sub-obj. 2.1: Trusted media and information",
                                                               sub_obj_num == "sub_obj_2_2" ~ "Sub-obj. 2.2: Media demand",
                                                               sub_obj_num == "sub_obj_2_3" ~ 
                                                                       "Sub-obj. 2.3: Media freedom")) %>%
                          select(country, obj_num, obj_avg, sub_obj_num_label, sub_obj_avg) %>%
                          pivot_wider(id_cols = c(country, obj_num, obj_avg), 
                                      names_from = sub_obj_num_label, values_from = sub_obj_avg) %>%
                          rename(!!sym("Obj. 2: Information") := obj_avg) %>% select(-obj_num), 
                  by = "country") %>%
        pivot_longer(cols = -c(country, cluster, country_order), names_to = "var", values_to = "values") %>%
        mutate(color = case_when(cluster == 1 ~ color_palette %>% slice(2) %>% pull(hex),
                                 cluster == 2 ~ color_palette %>% slice(7) %>% pull(hex),
                                 cluster == 3 ~ color_palette %>% slice(10) %>% pull(hex)),
               country_label = str_c("Group ", cluster, ": ", country),
               var_order = case_when(var == "Obj. 2: Information" ~ 1,
                                     str_detect(string = var, pattern = "^Sub-obj.") ~ 2,
                                     TRUE ~ 3))

# inspect
chart_data


#/////////////////////////////


# create chart
fmir_obj_2_country_indicator_tile_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       y = fct_reorder(.f = country_label, .x = country_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Score\n(higher is better)\n", y = NULL, x = "Indicator") +
        coord_fixed(ratio = 1 / 3, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 5, 
                                           color = chart_data %>% distinct(country, color) %>% pull(color), 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 5, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 5, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_2_country_indicator_tile_chart <- ggdraw(align_legend(fmir_obj_2_country_indicator_tile_chart))


# inspect
fmir_obj_2_country_indicator_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_2_country_indicator_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_2_country_indicator_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_2_indicator_cluster_boxplot_chart ####

# get chart_data
# note that this uses the original ISV values from the index, not the scaled ISV used in fmir_obj_2_indicator for clustering; 
# original/unscaled ISV avg is better for plot because values are then interpretable in the framework of the index (e.g. ISV btw 0-1)
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_2_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        
        
        # add cluster avg of sub_obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          left_join(., fmir_obj_2_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_2_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -c(country, cluster), names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster),
               color_bin = as.character(cluster_name),
               color = case_when(color_bin == "Group 1" ~ color_palette %>% slice(2) %>% pull(hex),
                                 color_bin == "Group 2" ~ color_palette %>% slice(7) %>% pull(hex),
                                 color_bin == "Group 3" ~ color_palette %>% slice(10) %>% pull(hex)))

# inspect
chart_data
chart_data %>% glimpse()

# create color_list for to pass to scale_color_manual
chart_data_color_list <- chart_data %>% count(color_bin, color) %>% pull(color)
names(chart_data_color_list) <- chart_data %>% count(color_bin, color) %>% pull(color_bin)
chart_data_color_list


#//////////////////////////


# create chart
fmir_obj_2_indicator_cluster_boxplot_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = values, fill = color_bin)) + 
        geom_boxplot(width = .5, lwd = .01, outlier.size = .5) + facet_wrap(facets = vars(cluster_name)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = .2), limits = c(0, 1), expand = c(0, 0)) +
        scale_fill_manual(values = chart_data_color_list) + 
        # scale_color_manual(values = chart_data_color_list) + 
        coord_fixed(ratio = 10/.6, clip = "off") +
        labs(x = "Indicator", y = "Score (higher is better)", fill = NULL, color = NULL) +
        theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                # strip.background = element_blank(),
                strip.background = element_rect(fill = "#DDDDDD", color = NA),
                strip.text = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333"),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = -15, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 5, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = 1,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 


# inspect
fmir_obj_2_indicator_cluster_boxplot_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_2_indicator_cluster_boxplot_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_2_indicator_cluster_boxplot_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_2_indicator_cluster_tile_chart ####

# get chart_data 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>% 
        
        
        # add clusters
        left_join(., fmir_obj_2_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        group_by(cluster) %>%
        mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
        ungroup() %>%
        select(-country) %>% distinct() %>%
        
        
        # add cluster avg of sub_obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          left_join(., fmir_obj_2_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_2_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -cluster, names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster))

# inspect
chart_data


#///////////////


# get obj/sub_obj segment coordinates

# note that in the entire fmir dataset, the maximum sub_obj per obj is 3, so will have 3 segment_tbls,
# and for those obj with less than 3 sub_obj, the segment_tbl will just have NA coordinates and no output in ggplot
fmir %>% distinct(obj_num, sub_obj_num) %>% 
        add_count(obj_num, name = "sub_obj_count") %>% 
        distinct(obj_num, sub_obj_count) %>%
        arrange(desc(sub_obj_count)) 

# segment_1_tbl
# segment_1_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(obj_num, sub_obj_num) %>%
#         mutate(nrow = nrow(.),
#                x = nrow + 1.5,
#                xend = nrow + 1.5,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5) %>%
#         select(x, xend, y, yend) %>%
#         slice(1)
# segment_1_tbl

# # segment_2_tbl
# segment_2_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                 cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(1)
# segment_2_tbl
#                
# # segment_3_tbl
# segment_3_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(2)
# segment_3_tbl


#///////////////////////


# plot

# get order_obj_and_sub_obj_avg_first()
order_obj_and_sub_obj_avg_first <- function(var){
        
        tibble(var = var) %>% 
                left_join(., fmir %>% distinct(indicator_name, concept) %>% 
                                  add_group_index(group_vars = concept, group_name = "concept_index"), 
                          by = c("var" = "indicator_name")) %>%
                mutate(order = case_when(str_detect(string = var, pattern = "^obj_") ~ 1,
                                         str_detect(string = var, pattern = "^sub_obj_[0-9]_[0-9]$") ~ 2,
                                         TRUE ~ 3)) %>%
                arrange(order, concept_index, var) %>% pull(var)
}

# create chart
fmir_obj_2_indicator_cluster_tile_chart <- chart_data %>% 
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = factor(cluster_name), fill = values, label = round_to_digits(values, digits = 2))) + 
        geom_tile(color = "#DDDDDD") + 
        geom_text(color = "#ffffff", size = 1.75, family = "Calibri") +
        scale_y_discrete(expand = c(0, 0)) +
        # geom_segment(x = segment_1_tbl$x, xend = segment_1_tbl$xend, 
        #              y = segment_1_tbl$y, yend = segment_1_tbl$yend, color = "#ffffff", size = 3) +
        # geom_segment(x = segment_2_tbl$x, xend = segment_2_tbl$xend,
        #              y = segment_2_tbl$y, yend = segment_2_tbl$yend, color = "#CBCBCB", size = 3) +
        # geom_segment(x = segment_3_tbl$x, xend = segment_3_tbl$xend,
        #              y = segment_3_tbl$y, yend = segment_3_tbl$yend, color = "#CBCBCB", size = 3) +
        labs(x = "Indicator", y = NULL, fill = "Avg. score\n(higher is better)\n") +
        coord_fixed(ratio = 2 / 1, clip = "off") +
        theme_bw() +
        # theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD"),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "right",
                legend.key.size = unit(4, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 5, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# fmir_obj_2_indicator_cluster_tile_chart <- ggdraw(align_legend(fmir_obj_2_indicator_cluster_tile_chart))

# inspect
fmir_obj_2_indicator_cluster_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_2_indicator_cluster_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_2_indicator_cluster_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_3_indicator_clusters ####

# set current_obj_num
current_obj_num <- "obj_3"

# get fmir_obj_3_indicator with scaled version of indicator_standardized_values
fmir_obj_3_indicator <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        distinct(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        column_to_rownames(var = "country") %>%
        mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1]))


#//////////////////////


# inspect
fmir_obj_3_indicator
fmir_obj_3_indicator %>% glimpse()
fmir_obj_3_indicator %>% nrow()
fmir_obj_3_indicator %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get fmir_obj_3_indicator_distance_matrix 
fmir_obj_3_indicator_distance_matrix <- dist(fmir_obj_3_indicator, method = "euclidean")
# distance_matrix <- get_dist(fmir_obj_3_indicator, method = "euclidean")


#//////////////


# inspect
fmir_obj_3_indicator_distance_matrix 
fmir_obj_3_indicator_distance_matrix %>% str()
fmir_obj_3_indicator_distance_matrix %>% attributes()

# inspect fmir_obj_3_distance_matrix_tbl
fmir_obj_3_indicator_distance_matrix_tbl <- as.matrix(fmir_obj_3_indicator_distance_matrix) %>% as_tibble() %>% 
        mutate(country = row.names(as.matrix(fmir_obj_3_indicator_distance_matrix))) %>%
        relocate(country, .before = everything())

fmir_obj_3_indicator_distance_matrix_tbl
fmir_obj_3_indicator_distance_matrix_tbl %>% skim() %>% as_tibble() %>% 
        select(skim_variable, numeric.mean, starts_with(match = "numeric.p")) %>%
        # arrange(numeric.mean) %>% 
        arrange(numeric.p25) %>%
        print(n = nrow(.))
fmir_obj_3_indicator_distance_matrix_tbl %>% 
        # select(country, Serbia) %>% arrange(Serbia) %>% print(n = nrow(.))
        # select(country, BiH) %>% arrange(BiH) %>% print(n = nrow(.))
        # select(country, Kosovo) %>% arrange(Kosovo) %>% print(n = nrow(.))
        # select(country, Armenia) %>% arrange(Armenia) %>% print(n = nrow(.))
        select(country, Poland) %>% arrange(Poland) %>% print(n = nrow(.))
# filter(country == "Serbia") %>% select(ends_with("stan"))

# compare indicator_standardized_values
fmir %>% filter(obj_num == "obj_3", country %in% c("Poland", "Kazakhstan"), year == 2020) %>% 
        select(country, year, indicator_name, values, values_z_std, indicator_standardized_values) %>%
        arrange(country) %>% 
        # unite(col = var, indicator_name, country) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = indicator_name, names_from = country, values_from = indicator_standardized_values)

# visualize distance
viridis_pal()(3) %>% show_col()
fviz_dist(dist.obj = fmir_obj_3_indicator_distance_matrix, gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_3_indicator_distance_matrix_tile_chart ####


# get country_order from fviz_dist_plot
fmir_obj_3_fviz_dist_plot <- fviz_dist(dist.obj = fmir_obj_3_indicator_distance_matrix, 
                                       gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))
fmir_obj_3_fviz_dist_plot
fmir_obj_3_fviz_dist_plot %>% attributes()
fmir_obj_3_fviz_dist_plot$data %>% as_tibble()

fmir_obj_3_country_order_tbl <- fmir_obj_3_fviz_dist_plot$data %>% as_tibble() %>% slice(1:(fmir_obj_3_indicator %>% nrow())) %>% 
        mutate(country = str_sub(string = Var1, start = 1, end = -2),
               country_order = row_number()) %>% 
        select(country, country_order)
fmir_obj_3_country_order_tbl


#/////////////////////////////


# create chart
fmir_obj_3_indicator_distance_matrix_tile_chart <- fmir_obj_3_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "var", values_to = "values") %>%
        left_join(., fmir_obj_3_country_order_tbl, by = c("var" = "country")) %>%
        rename(var_order = country_order) %>%
        left_join(., fmir_obj_3_country_order_tbl, by = "country") %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = country, .x = country_order, .fun = min), 
                                       y = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Similiarity\n(lower is more similar)", y = NULL, x = NULL) +
        coord_fixed(ratio = 1/1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 13, r = 0, b = 5, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_3_indicator_distance_matrix_tile_chart <- ggdraw(align_legend(fmir_obj_3_indicator_distance_matrix_tile_chart))


# inspect
fmir_obj_3_indicator_distance_matrix_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_3_indicator_distance_matrix_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_3_indicator_distance_matrix_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_3_indicator_similarity_table ####

fmir_obj_3_indicator_similarity_table <- fmir_obj_3_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "counterpart", values_to = "distance") %>%
        filter(country != counterpart) %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping), by = "country") %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping) %>% rename(counterpart_mcp_grouping = mcp_grouping), 
                  by = c("counterpart" = "country")) %>%
        group_by(country) %>%
        mutate(distance_in_ee = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_flag = case_when(distance == min(distance) ~ 1,
                                               TRUE ~ 0),
               max_similarity = case_when(max_similarity_flag == 1 ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_counterpart = case_when(max_similarity_flag == 1 ~ counterpart,
                                                      TRUE ~ NA_character_),
               max_similarity_in_ee_flag = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") &
                                                             distance_in_ee == min(distance_in_ee, na.rm = TRUE) ~ 1,
                                                     TRUE ~ NA_real_),
               max_similarity_in_ee = case_when(max_similarity_in_ee_flag == 1 ~ distance_in_ee,
                                                TRUE ~ NA_real_),
               max_similarity_in_ee_counterpart = case_when(max_similarity_in_ee_flag == 1 ~ counterpart,
                                                            TRUE ~ NA_character_)) %>%
        fill(max_similarity, .direction = "updown") %>%
        fill(max_similarity_counterpart, .direction = "updown") %>%
        fill(max_similarity_in_ee, .direction = "updown") %>%
        fill(max_similarity_in_ee_counterpart, .direction = "updown") %>%
        ungroup()


#//////////////////


# inspect
fmir_obj_3_indicator_similarity_table
fmir_obj_3_indicator_similarity_table %>% glimpse()
fmir_obj_3_indicator_similarity_table %>% nrow() # 756
fmir_obj_3_indicator_similarity_table %>% ncol() # 12


fmir_obj_3_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping)) %>% 
        print(n = 30)
fmir_obj_3_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping, max_similarity, max_similarity_counterpart)) %>% 
        print(n = 30)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get distinct max_similarity/_in_ee values/counterparts and text_output for table
fmir_obj_3_indicator_similarity_table <- fmir_obj_3_indicator_similarity_table %>% 
        distinct(country, mcp_grouping, max_similarity, max_similarity_counterpart,
                 max_similarity_in_ee, max_similarity_in_ee_counterpart) %>%
        mutate(max_similarity_text_output = str_c("   ", max_similarity_counterpart, " (", round_to_digits(max_similarity, digits = 2), ")"),
               max_similarity_in_ee_text_output = str_c("   ", max_similarity_in_ee_counterpart, " (", 
                                                        round_to_digits(max_similarity_in_ee, digits = 2), ")"),
               mcp_grouping_order = case_when(mcp_grouping == "E&E Balkans" ~ 1,
                                              mcp_grouping == "E&E Eurasia" ~ 2,
                                              mcp_grouping == "CARs" ~ 3,
                                              mcp_grouping == "E&E graduates" ~ 4)) %>%
        arrange(mcp_grouping_order, country) %>%
        select(country, mcp_grouping, 
               max_similarity_text_output, max_similarity, 
               max_similarity_in_ee_text_output, max_similarity_in_ee)


#///////////////////


# inspect
fmir_obj_3_indicator_similarity_table
fmir_obj_3_indicator_similarity_table %>% glimpse()
fmir_obj_3_indicator_similarity_table %>% nrow() # 28
fmir_obj_3_indicator_similarity_table %>% ncol() # 5

fmir_obj_3_indicator_similarity_table %>% print(n = nrow(.))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get table_palette


# use viridis palette
viridis_pal()(10) %>% show_col()


# historic
# note viridis is not used to avoid confusing color scales when compared to dissimilarity matrix, which has much higher
# range of values than the similarity_table
# table_palette <- colorRampPalette(c(brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(4) %>% pull(value), 
#                                      brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(9) %>% pull(value)))


#////////////////////////


# inspect
viridis_pal()(10) %>% show_col()
brewer.pal(n = 9, "Blues") %>% show_col()
table_palette(5) %>% show_col()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get values_normalized_quantile_cutoffs, creating one set of color thresholds to be applied uniformly to all obj in the table
# this just makes it easier to automatically assign color-coded quantiles for the full range of values for E&E countries
# whereas using the actual values values themselves would require checks after updates to make sure the cutoffs are still balanced
# across the full range of values values
# note that this doesn't include miri_avg, but the table does, but this is fine because miri_avg by definition must fall within
# range of max/min values, and since this only gets the color quantiles based on max/min values, it's fine to exclude miri_avg
values_normalized_quantile_cutoffs <- fmir_obj_3_indicator_similarity_table %>%
        pivot_longer(cols = c(max_similarity, max_similarity_in_ee), names_to = "var", values_to = "values") %>%
        mutate(values_min = min(values),
               values_max = max(values),
               values_normalized = (values - values_min) / (values_max - values_min),
               q0_to_q10 = .10 * (values_max - values_min) + values_min,
               q11_to_q20 = .20 * (values_max - values_min) + values_min,
               q21_to_q30 = .30 * (values_max - values_min) + values_min,
               q31_to_q40 = .40 * (values_max - values_min) + values_min,
               q41_to_q50 = .50 * (values_max - values_min) + values_min,
               q51_to_q60 = .60 * (values_max - values_min) + values_min,
               q61_to_q70 = .70 * (values_max - values_min) + values_min,
               q71_to_q80 = .80 * (values_max - values_min) + values_min,
               q81_to_q90 = .90 * (values_max - values_min) + values_min,
               q91_to_q100 = 1 * (values_max - values_min) + values_min) %>%
        select(-c(max_similarity_text_output, max_similarity_in_ee_text_output))


#/////////////////////////


# inspect
values_normalized_quantile_cutoffs %>% data.frame()
# values_normalized = (values - values_min) / (values_max - values_min)

# check q0-q10 threshold
(2.25 - 1.81) / (6.18 - 1.81)

# check q81-q90 threshold
(5.74 - 1.81) / (6.18 - 1.81)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# create table_palette_tbl
# manually apply conditional formatting in excel to highlight and change font color
# =IF(AND(B4>5, B4<10),TRUE,FALSE)
table_palette_tbl <- values_normalized_quantile_cutoffs %>% select(starts_with("q")) %>% distinct() %>%
        pivot_longer(cols = everything(), names_to = "quantile", values_to = "values_threshold") %>%
        left_join(., tibble(hex = viridis_pal()(10), 
                            quantile = c("q0_to_q10", "q11_to_q20", "q21_to_q30", "q31_to_q40", "q41_to_q50",
                                         "q51_to_q60", "q61_to_q70", "q71_to_q80", "q81_to_q90", "q91_to_q100")) %>%
                          mutate(rgb = map(.x = hex, .f = ~ col2rgb(col = .x) %>% str_c(., collapse = "."))) %>%
                          unnest(rgb),
                  by = "quantile")

# inspect
table_palette_tbl 


#//////////////////////////////////////////////////////////////////////////////////////////////////


# write to file

# note csv will be copied into xlsx shell and formatted
# old: blue header is rgb 31/73/125
# current: gray header is rgb 128/128/128
# grey grid lines are rgb 217/217/217, medium weight
# country rows on left side gets a rgb 217/217/217 heavy border
# grey top/bottom row above/under table is rgb 217/217/217, height = 15
# font is Calibri, size 15, bolded, centered horizontally/vertically in cells; 
# grey font color = rgb 51/51/51, white for column headers
# flag column is width = 15
# country name column is width = 32
# data columns are width = 28
# column header row height = 95
# data row height = 35
# data cells formatted as number with two decimal places, to avoid truncation (e.g. 0.70 converted to 0.7)

# remember to fully spell out macedonia and bih
# output is copied into powerpoint using paste special -> Picture (Enhanced Metafile), aligned center
# flags are copy/pasted in from state dept website, w locked aspect ratio, and resized to 70%

# save
# fmir_obj_3_indicator_similarity_table %>% write_csv(file = "output/charts/fmir_obj_3_indicator_similarity_table.csv")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# get cluster tendency to see if data is likely to have clusters that are more meaningful/tight than uniform random sample
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# note that the hopkins stat uses .5 as a threshold, with higher values indicating there is a clustering tendency; 
# a hopkins stat over .75 is 90% confidence that there is clustering tendency that is real, and not random
cluster_tendency <- get_clust_tendency(fmir_obj_3_indicator, n = nrow(fmir_obj_3_indicator) - 1, graph = FALSE)


#///////////////////


# inspect
cluster_tendency
cluster_tendency$hopkins_stat


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters
fmir_obj_3_indicator_hclust <- hclust(d = fmir_obj_3_indicator_distance_matrix, method = "complete" )


#/////////////////////


# inspect
fmir_obj_3_indicator_hclust
fmir_obj_3_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_3_indicator_hclust$merge
fmir_obj_3_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_3_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters with agnes, which allows enhanced fviz visualizations
fmir_obj_3_indicator_agnes <- agnes(x = fmir_obj_3_indicator, method = "complete")


#//////////////////


# inspect
fmir_obj_3_indicator_agnes
fmir_obj_3_indicator_agnes %>% attributes()

# get agglomerative coefficient from agnes object 
# "Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 suggest a 
# more balanced clustering structure such as the complete linkage and Wards method dendrograms in Figure 21.3. 
# Values closer to 0 suggest less well-formed clusters such as the single linkage dendrogram in Figure 21.3. 
# However, the AC tends to become larger as n increases, so it should not be used to compare across data sets of very different sizes."
fmir_obj_3_indicator_agnes$ac

# plot
pltree(fmir_obj_3_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 



#//////////////////////////////////////////////////////////////////////////////////////////////////


# use agnes to look for clustering method with best agglomerative coefficient 
method_list <- c( "average", "single", "complete", "ward")
names(method_list) <- c( "average", "single", "complete", "ward")
method_list

# function to compute coefficient
cluster_methods <- function(current_method) {
        print(current_method)
        agnes(fmir_obj_3_indicator_distance_matrix, method = current_method)$ac
}

# run cluster_methods
# note that wards method has the highest/best agglomerative coefficient
map(.x = method_list, .f = ~ cluster_methods(current_method = .x))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with agnes, which had highest (best) agglomerative coefficient
fmir_obj_3_indicator_agnes <- agnes(fmir_obj_3_indicator, method = "ward")
fmir_obj_3_indicator_agnes


#///////////////////


# visualize clustering 
pltree(fmir_obj_3_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
fviz_dend(fmir_obj_3_indicator_agnes, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

# inspect optimal number of clusters
# note the gap statistic is recently developed by tibshirani et al
fviz_nbclust(fmir_obj_3_indicator, FUN = hcut, method = "wss")
fviz_nbclust(fmir_obj_3_indicator, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(fmir_obj_3_indicator, FUN = hcut, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# set optimal k
fmir_obj_3_indicator_k <- 3


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with hclust
fmir_obj_3_indicator_hclust <- hclust(d = fmir_obj_3_indicator_distance_matrix, method = "ward.D2")


#/////////////////////


# inspect
fmir_obj_3_indicator_hclust
fmir_obj_3_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_3_indicator_hclust$merge
fmir_obj_3_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_3_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut tree 

# build and cut agnes tree with hcut for enhanced fviz
fmir_obj_3_indicator_agnes_cut <- hcut(fmir_obj_3_indicator, k = fmir_obj_3_indicator_k, hc_method = "ward.D2", stand = TRUE)
fmir_obj_3_indicator_agnes_cut 

# view PCA for clusters, which requires an hcut object
fviz_cluster(object = fmir_obj_3_indicator_agnes_cut)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut hclust tree with cutree
fmir_obj_3_indicator_hclust_cut <- cutree(fmir_obj_3_indicator_hclust, k = fmir_obj_3_indicator_k)
fmir_obj_3_indicator_hclust_cut

# Number of members in each cluster
table(fmir_obj_3_indicator_hclust_cut)

# can add cutree output to df
fmir_obj_3_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value) %>%
        left_join(fmir_obj_3_indicator %>% rownames_to_column(var = "country"), ., by = "country")

# plot dendrogram with border around 4 clusters
plot(fmir_obj_3_indicator_hclust, cex = 0.6, hang = -1)
rect.hclust(fmir_obj_3_indicator_hclust, k = fmir_obj_3_indicator_k, border = 2:5)

# plot dendrogram w labels
dendrogram_obj <- as.dendrogram(fmir_obj_3_indicator_hclust)
dendrogram_obj %>% unclass() %>% str()
dendrogram_obj %>% unclass() %>% class()

subdendrogram_list <- get_subdendrograms(dendrogram_obj, k = fmir_obj_3_indicator_k)
subdendrogram_list
dendrogram_labels <- map(.x = subdendrogram_list, .f = ~ labels(.x)) %>% 
        enframe() %>% unnest(cols = value) %>% 
        rename(subdendrogram_plot_order = name, country = value) %>%
        left_join(fmir_obj_3_indicator %>% mutate(cluster = fmir_obj_3_indicator_hclust_cut) %>%
                          rownames_to_column(var = "country"),
                  ., by = "country") %>% 
        distinct(cluster, subdendrogram_plot_order) %>% 
        arrange(subdendrogram_plot_order) %>%
        mutate(cluster = str_c("\n\n\n", cluster)) %>%
        pull(cluster) %>% unname()
dendrogram_labels        

plot(fmir_obj_3_indicator_hclust, hang = -1)
rect.dendrogram(as.dendrogram(fmir_obj_3_indicator_hclust), k = fmir_obj_3_indicator_k, text = dendrogram_labels)

# visualize clusters along 2 primary PCA dimensions
fviz_cluster(list(data = fmir_obj_3_indicator, cluster = fmir_obj_3_indicator_hclust_cut))


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_3_indicator_dendrogram ####

# get fmir_obj_3_indicator_dendrogram_data
# note that dendro_data_k() takes the uncut hclust as first argument
fmir_obj_3_indicator_dendrogram_data <- dendro_data_k(hc = fmir_obj_3_indicator_hclust, 
                                                      k = fmir_obj_3_indicator_k)

# inspect
fmir_obj_3_indicator_dendrogram_data
fmir_obj_3_indicator_dendrogram_data %>% attributes()
fmir_obj_3_indicator_dendrogram_data$segments %>% head()
fmir_obj_3_indicator_dendrogram_data$labels %>% head()
fmir_obj_3_indicator_dendrogram_data$leaf_labels %>% head()
fmir_obj_3_indicator_dendrogram_data$class %>% head()


#///////////////////


# get cluster_colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
cluster_colors <- c("#333333", 
                    color_palette %>% slice(2) %>% pull(hex), 
                    color_palette %>% slice(7) %>% pull(hex),
                    color_palette %>% slice(10) %>% pull(hex))
cluster_colors


# plot
fmir_obj_3_indicator_dendrogram <- plot_ggdendro(fmir_obj_3_indicator_dendrogram_data,
                                                 direction = "tb",
                                                 scale.color = cluster_colors,
                                                 label.size = 2.5,
                                                 branch.size = 0.5,
                                                 expand.y = 0) +
        labs(y = "Similarity\n(lower groupings are more similar)", x = "") +
        coord_fixed(ratio = 2 / 1.75, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 5, 0, 0), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(size = .5),
                axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "bold", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 0, hjust = .5),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                # axis.title.x = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 5, r = 0, b = 10, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 0, r = 15, b = 0, l = 5)),
                # axis.title.y = element_blank(),
                # axis.text.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 9, family = "Calibri", face = "plain", color = "#333333"),
                legend.text = element_text(size = 9, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# guides(color = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333", keywidth = 4),
#        linetype = guide_legend(keywidth = 4))

# inspect
fmir_obj_3_indicator_dendrogram


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_3_indicator_dendrogram)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_3_indicator_dendrogram.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_3_country_indicator_tile_chart ####

# get chart_data
# note this is 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        left_join(., fmir_obj_1_country_order_tbl, by = "country") %>%
        relocate(country_order, .after = cluster) %>% arrange(country_order) %>%
        pivot_wider(id_cols = c(country, cluster, country_order), 
                    names_from = indicator_name, values_from = indicator_standardized_values) %>%
        left_join(., fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          distinct(country, obj_num, sub_obj_num, sub_obj_avg, obj_avg) %>%
                          mutate(sub_obj_num_label = case_when(sub_obj_num == "sub_obj_3_1" ~ 
                                                                       "Sub-obj. 3.1: Energy security",
                                                               sub_obj_num == "sub_obj_3_2" ~ "Sub-obj. 3.2: Independence from Russian energy",
                                                               sub_obj_num == "sub_obj_3_3" ~ 
                                                                       "Sub-obj. 3.3: Energy regulation")) %>%
                          select(country, obj_num, obj_avg, sub_obj_num_label, sub_obj_avg) %>%
                          pivot_wider(id_cols = c(country, obj_num, obj_avg), 
                                      names_from = sub_obj_num_label, values_from = sub_obj_avg) %>%
                          rename(!!sym("Obj. 3: Energy") := obj_avg) %>% select(-obj_num), 
                  by = "country") %>%
        pivot_longer(cols = -c(country, cluster, country_order), names_to = "var", values_to = "values") %>%
        mutate(color = case_when(cluster == 1 ~ color_palette %>% slice(2) %>% pull(hex),
                                 cluster == 2 ~ color_palette %>% slice(7) %>% pull(hex),
                                 cluster == 3 ~ color_palette %>% slice(10) %>% pull(hex)),
               country_label = str_c("Group ", cluster, ": ", country),
               var_order = case_when(var == "Obj. 3: Energy" ~ 1,
                                     str_detect(string = var, pattern = "^Sub-obj.") ~ 2,
                                     TRUE ~ 3))

# inspect
chart_data


#/////////////////////////////


# create chart
fmir_obj_3_country_indicator_tile_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       y = fct_reorder(.f = country_label, .x = country_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Score\n(higher is better)\n", y = NULL, x = "Indicator") +
        coord_fixed(ratio = 1 / 7, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 5, 
                                           color = chart_data %>% distinct(country, color) %>% pull(color), 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 5, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 5, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_3_country_indicator_tile_chart <- ggdraw(align_legend(fmir_obj_3_country_indicator_tile_chart))


# inspect
fmir_obj_3_country_indicator_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_3_country_indicator_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_3_country_indicator_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_3_indicator_cluster_boxplot_chart ####

# get chart_data
# note that this uses the original ISV values from the index, not the scaled ISV used in fmir_obj_3_indicator for clustering; 
# original/unscaled ISV avg is better for plot because values are then interpretable in the framework of the index (e.g. ISV btw 0-1)
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_3_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        
        
        # add cluster avg of sub_obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          left_join(., fmir_obj_3_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_3_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -c(country, cluster), names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster),
               color_bin = as.character(cluster_name),
               color = case_when(color_bin == "Group 1" ~ color_palette %>% slice(2) %>% pull(hex),
                                 color_bin == "Group 2" ~ color_palette %>% slice(7) %>% pull(hex),
                                 color_bin == "Group 3" ~ color_palette %>% slice(10) %>% pull(hex)))

# inspect
chart_data
chart_data %>% glimpse()

# create color_list for to pass to scale_color_manual
chart_data_color_list <- chart_data %>% count(color_bin, color) %>% pull(color)
names(chart_data_color_list) <- chart_data %>% count(color_bin, color) %>% pull(color_bin)
chart_data_color_list


#//////////////////////////


# create chart
fmir_obj_3_indicator_cluster_boxplot_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = values, fill = color_bin)) + 
        geom_boxplot(width = .5, lwd = .01, outlier.size = .5) + facet_wrap(facets = vars(cluster_name)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = .2), limits = c(0, 1), expand = c(0, 0)) +
        scale_fill_manual(values = chart_data_color_list) + 
        # scale_color_manual(values = chart_data_color_list) + 
        coord_fixed(ratio = 10/1.5, clip = "off") +
        labs(x = "Indicator", y = "Score (higher is better)", fill = NULL, color = NULL) +
        theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                # strip.background = element_blank(),
                strip.background = element_rect(fill = "#DDDDDD", color = NA),
                strip.text = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333"),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = -15, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 5, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = 1,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 


# inspect
fmir_obj_3_indicator_cluster_boxplot_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_3_indicator_cluster_boxplot_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_3_indicator_cluster_boxplot_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_3_indicator_cluster_tile_chart ####

# get chart_data 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_3_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        group_by(cluster) %>%
        mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
        ungroup() %>%
        select(-country) %>% distinct() %>%
        
        
        # add cluster avg of sub_obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          left_join(., fmir_obj_3_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_3_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -cluster, names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster))

# inspect
chart_data


#///////////////


# get obj/sub_obj segment coordinates

# note that in the entire fmir dataset, the maximum sub_obj per obj is 3, so will have 3 segment_tbls,
# and for those obj with less than 3 sub_obj, the segment_tbl will just have NA coordinates and no output in ggplot
fmir %>% distinct(obj_num, sub_obj_num) %>% 
        add_count(obj_num, name = "sub_obj_count") %>% 
        distinct(obj_num, sub_obj_count) %>%
        arrange(desc(sub_obj_count)) 

# segment_1_tbl
# segment_1_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(obj_num, sub_obj_num) %>%
#         mutate(nrow = nrow(.),
#                x = nrow + 1.5,
#                xend = nrow + 1.5,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5) %>%
#         select(x, xend, y, yend) %>%
#         slice(1)
# segment_1_tbl

# # segment_2_tbl
# segment_2_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                 cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(1)
# segment_2_tbl
#                
# # segment_3_tbl
# segment_3_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(2)
# segment_3_tbl


#///////////////////////


# plot

# get order_obj_and_sub_obj_avg_first()
order_obj_and_sub_obj_avg_first <- function(var){
        
        tibble(var = var) %>% 
                left_join(., fmir %>% distinct(indicator_name, concept) %>% 
                                  add_group_index(group_vars = concept, group_name = "concept_index"), 
                          by = c("var" = "indicator_name")) %>%
                mutate(order = case_when(str_detect(string = var, pattern = "^obj_") ~ 1,
                                         str_detect(string = var, pattern = "^sub_obj_[0-9]_[0-9]$") ~ 2,
                                         TRUE ~ 3)) %>%
                arrange(order, concept_index, var) %>% pull(var)
}

# create chart
fmir_obj_3_indicator_cluster_tile_chart <- chart_data %>% 
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = factor(cluster_name), fill = values, label = round_to_digits(values, digits = 2))) + 
        geom_tile(color = "#DDDDDD") + 
        geom_text(color = "#ffffff", size = 1.75, family = "Calibri") +
        scale_y_discrete(expand = c(0, 0)) +
        # geom_segment(x = segment_1_tbl$x, xend = segment_1_tbl$xend, 
        #              y = segment_1_tbl$y, yend = segment_1_tbl$yend, color = "#ffffff", size = 3) +
        # geom_segment(x = segment_2_tbl$x, xend = segment_2_tbl$xend,
        #              y = segment_2_tbl$y, yend = segment_2_tbl$yend, color = "#CBCBCB", size = 3) +
        # geom_segment(x = segment_3_tbl$x, xend = segment_3_tbl$xend,
        #              y = segment_3_tbl$y, yend = segment_3_tbl$yend, color = "#CBCBCB", size = 3) +
        labs(x = "Indicator", y = NULL, fill = "Avg. score\n(higher is better)\n") +
        coord_fixed(ratio = 2 / 3, clip = "off") +
        theme_bw() +
        # theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD"),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "right",
                legend.key.size = unit(4, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 5, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# fmir_obj_3_indicator_cluster_tile_chart <- ggdraw(align_legend(fmir_obj_3_indicator_cluster_tile_chart))

# inspect
fmir_obj_3_indicator_cluster_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_3_indicator_cluster_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_3_indicator_cluster_tile_chart.docx")



#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_4_indicator_clusters ####

# set current_obj_num
current_obj_num <- "obj_4"

# get fmir_obj_4_indicator with scaled version of indicator_standardized_values
fmir_obj_4_indicator <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        distinct(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        column_to_rownames(var = "country") %>%
        mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1]))


#//////////////////////


# inspect
fmir_obj_4_indicator
fmir_obj_4_indicator %>% glimpse()
fmir_obj_4_indicator %>% nrow()
fmir_obj_4_indicator %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get fmir_obj_4_indicator_distance_matrix 
fmir_obj_4_indicator_distance_matrix <- dist(fmir_obj_4_indicator, method = "euclidean")
# fmir_obj_4_indicator_distance_matrix <- get_dist(fmir_obj_4_indicator, method = "euclidean")


#//////////////


# inspect
fmir_obj_4_indicator_distance_matrix
fmir_obj_4_indicator_distance_matrix %>% str()
fmir_obj_4_indicator_distance_matrix %>% attributes()

# inspect fmir_obj_4_indicator_distance_matrix_tbl
fmir_obj_4_indicator_distance_matrix_tbl <- as.matrix(fmir_obj_4_indicator_distance_matrix) %>% as_tibble() %>% 
        mutate(country = row.names(as.matrix(fmir_obj_4_indicator_distance_matrix))) %>%
        relocate(country, .before = everything())

fmir_obj_4_indicator_distance_matrix_tbl
fmir_obj_4_indicator_distance_matrix_tbl %>% skim() %>% as_tibble() %>% 
        select(skim_variable, numeric.mean, starts_with(match = "numeric.p")) %>%
        # arrange(numeric.mean) %>% 
        arrange(numeric.p25) %>%
        print(n = nrow(.))
fmir_obj_4_indicator_distance_matrix_tbl %>% 
        # select(country, Serbia) %>% arrange(Serbia) %>% print(n = nrow(.))
        # select(country, BiH) %>% arrange(BiH) %>% print(n = nrow(.))
        # select(country, Kosovo) %>% arrange(Kosovo) %>% print(n = nrow(.))
        # select(country, Armenia) %>% arrange(Armenia) %>% print(n = nrow(.))
        select(country, Poland) %>% arrange(Poland) %>% print(n = nrow(.))
# filter(country == "Serbia") %>% select(ends_with("stan"))

# compare indicator_standardized_values
fmir %>% filter(obj_num == "obj_3", country %in% c("Poland", "Kazakhstan"), year == 2020) %>% 
        select(country, year, indicator_name, values, values_z_std, indicator_standardized_values) %>%
        arrange(country) %>% 
        # unite(col = var, indicator_name, country) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = indicator_name, names_from = country, values_from = indicator_standardized_values)

# visualize distance
viridis_pal()(3) %>% show_col()
fviz_dist(dist.obj = fmir_obj_4_indicator_distance_matrix, gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_4_indicator_distance_matrix_tile_chart ####


# get country_order from fviz_dist_plot
fmir_obj_4_fviz_dist_plot <- fviz_dist(dist.obj = fmir_obj_4_indicator_distance_matrix, 
                                       gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))
fmir_obj_4_fviz_dist_plot
fmir_obj_4_fviz_dist_plot %>% attributes()
fmir_obj_4_fviz_dist_plot$data %>% as_tibble()

fmir_obj_4_country_order_tbl <- fmir_obj_4_fviz_dist_plot$data %>% as_tibble() %>% slice(1:(fmir_obj_4_indicator %>% nrow())) %>% 
        mutate(country = str_sub(string = Var1, start = 1, end = -2),
               country_order = row_number()) %>% 
        select(country, country_order)
fmir_obj_4_country_order_tbl


#/////////////////////////////


# create chart
fmir_obj_4_indicator_distance_matrix_tile_chart <- fmir_obj_4_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "var", values_to = "values") %>%
        left_join(., fmir_obj_4_country_order_tbl, by = c("var" = "country")) %>%
        rename(var_order = country_order) %>%
        left_join(., fmir_obj_4_country_order_tbl, by = "country") %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = country, .x = country_order, .fun = min), 
                                       y = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Similiarity\n(lower is more similar)", y = NULL, x = NULL) +
        coord_fixed(ratio = 1/1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 13, r = 0, b = 5, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_4_indicator_distance_matrix_tile_chart <- ggdraw(align_legend(fmir_obj_4_indicator_distance_matrix_tile_chart))


# inspect
fmir_obj_4_indicator_distance_matrix_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_4_indicator_distance_matrix_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_4_indicator_distance_matrix_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_4_indicator_similarity_table ####

fmir_obj_4_indicator_similarity_table <- fmir_obj_4_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "counterpart", values_to = "distance") %>%
        filter(country != counterpart) %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping), by = "country") %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping) %>% rename(counterpart_mcp_grouping = mcp_grouping), 
                  by = c("counterpart" = "country")) %>%
        group_by(country) %>%
        mutate(distance_in_ee = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_flag = case_when(distance == min(distance) ~ 1,
                                               TRUE ~ 0),
               max_similarity = case_when(max_similarity_flag == 1 ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_counterpart = case_when(max_similarity_flag == 1 ~ counterpart,
                                                      TRUE ~ NA_character_),
               max_similarity_in_ee_flag = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") &
                                                             distance_in_ee == min(distance_in_ee, na.rm = TRUE) ~ 1,
                                                     TRUE ~ NA_real_),
               max_similarity_in_ee = case_when(max_similarity_in_ee_flag == 1 ~ distance_in_ee,
                                                TRUE ~ NA_real_),
               max_similarity_in_ee_counterpart = case_when(max_similarity_in_ee_flag == 1 ~ counterpart,
                                                            TRUE ~ NA_character_)) %>%
        fill(max_similarity, .direction = "updown") %>%
        fill(max_similarity_counterpart, .direction = "updown") %>%
        fill(max_similarity_in_ee, .direction = "updown") %>%
        fill(max_similarity_in_ee_counterpart, .direction = "updown") %>%
        ungroup()


#//////////////////


# inspect
fmir_obj_4_indicator_similarity_table
fmir_obj_4_indicator_similarity_table %>% glimpse()
fmir_obj_4_indicator_similarity_table %>% nrow() # 756
fmir_obj_4_indicator_similarity_table %>% ncol() # 12


fmir_obj_4_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping)) %>% 
        print(n = 30)
fmir_obj_4_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping, max_similarity, max_similarity_counterpart)) %>% 
        print(n = 30)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get distinct max_similarity/_in_ee values/counterparts and text_output for table
fmir_obj_4_indicator_similarity_table <- fmir_obj_4_indicator_similarity_table %>% 
        distinct(country, mcp_grouping, max_similarity, max_similarity_counterpart,
                 max_similarity_in_ee, max_similarity_in_ee_counterpart) %>%
        mutate(max_similarity_text_output = str_c("   ", max_similarity_counterpart, " (", round_to_digits(max_similarity, digits = 2), ")"),
               max_similarity_in_ee_text_output = str_c("   ", max_similarity_in_ee_counterpart, " (", 
                                                        round_to_digits(max_similarity_in_ee, digits = 2), ")"),
               mcp_grouping_order = case_when(mcp_grouping == "E&E Balkans" ~ 1,
                                              mcp_grouping == "E&E Eurasia" ~ 2,
                                              mcp_grouping == "CARs" ~ 3,
                                              mcp_grouping == "E&E graduates" ~ 4)) %>%
        arrange(mcp_grouping_order, country) %>%
        select(country, mcp_grouping, 
               max_similarity_text_output, max_similarity, 
               max_similarity_in_ee_text_output, max_similarity_in_ee)


#///////////////////


# inspect
fmir_obj_4_indicator_similarity_table
fmir_obj_4_indicator_similarity_table %>% glimpse()
fmir_obj_4_indicator_similarity_table %>% nrow() # 28
fmir_obj_4_indicator_similarity_table %>% ncol() # 5

fmir_obj_4_indicator_similarity_table %>% print(n = nrow(.))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get table_palette


# use viridis palette
viridis_pal()(10) %>% show_col()


# historic
# note viridis is not used to avoid confusing color scales when compared to dissimilarity matrix, which has much higher
# range of values than the similarity_table
# table_palette <- colorRampPalette(c(brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(4) %>% pull(value), 
#                                      brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(9) %>% pull(value)))


#////////////////////////


# inspect
viridis_pal()(10) %>% show_col()
brewer.pal(n = 9, "Blues") %>% show_col()
table_palette(5) %>% show_col()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get values_normalized_quantile_cutoffs, creating one set of color thresholds to be applied uniformly to all obj in the table
# this just makes it easier to automatically assign color-coded quantiles for the full range of values for E&E countries
# whereas using the actual values values themselves would require checks after updates to make sure the cutoffs are still balanced
# across the full range of values values
# note that this doesn't include miri_avg, but the table does, but this is fine because miri_avg by definition must fall within
# range of max/min values, and since this only gets the color quantiles based on max/min values, it's fine to exclude miri_avg
values_normalized_quantile_cutoffs <- fmir_obj_4_indicator_similarity_table %>%
        pivot_longer(cols = c(max_similarity, max_similarity_in_ee), names_to = "var", values_to = "values") %>%
        mutate(values_min = min(values),
               values_max = max(values),
               values_normalized = (values - values_min) / (values_max - values_min),
               q0_to_q10 = .10 * (values_max - values_min) + values_min,
               q11_to_q20 = .20 * (values_max - values_min) + values_min,
               q21_to_q30 = .30 * (values_max - values_min) + values_min,
               q31_to_q40 = .40 * (values_max - values_min) + values_min,
               q41_to_q50 = .50 * (values_max - values_min) + values_min,
               q51_to_q60 = .60 * (values_max - values_min) + values_min,
               q61_to_q70 = .70 * (values_max - values_min) + values_min,
               q71_to_q80 = .80 * (values_max - values_min) + values_min,
               q81_to_q90 = .90 * (values_max - values_min) + values_min,
               q91_to_q100 = 1 * (values_max - values_min) + values_min) %>%
        select(-c(max_similarity_text_output, max_similarity_in_ee_text_output))


#/////////////////////////


# inspect
values_normalized_quantile_cutoffs %>% data.frame()
# values_normalized = (values - values_min) / (values_max - values_min)

# check q0-q10 threshold
(2.25 - 1.81) / (6.18 - 1.81)

# check q81-q90 threshold
(5.74 - 1.81) / (6.18 - 1.81)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# create table_palette_tbl
# manually apply conditional formatting in excel to highlight and change font color
# =IF(AND(B4>5, B4<10),TRUE,FALSE)
table_palette_tbl <- values_normalized_quantile_cutoffs %>% select(starts_with("q")) %>% distinct() %>%
        pivot_longer(cols = everything(), names_to = "quantile", values_to = "values_threshold") %>%
        left_join(., tibble(hex = viridis_pal()(10), 
                            quantile = c("q0_to_q10", "q11_to_q20", "q21_to_q30", "q31_to_q40", "q41_to_q50",
                                         "q51_to_q60", "q61_to_q70", "q71_to_q80", "q81_to_q90", "q91_to_q100")) %>%
                          mutate(rgb = map(.x = hex, .f = ~ col2rgb(col = .x) %>% str_c(., collapse = "."))) %>%
                          unnest(rgb),
                  by = "quantile")

# inspect
table_palette_tbl 


#//////////////////////////////////////////////////////////////////////////////////////////////////


# write to file

# note csv will be copied into xlsx shell and formatted
# old: blue header is rgb 31/73/125
# current: gray header is rgb 128/128/128
# grey grid lines are rgb 217/217/217, medium weight
# country rows on left side gets a rgb 217/217/217 heavy border
# grey top/bottom row above/under table is rgb 217/217/217, height = 15
# font is Calibri, size 15, bolded, centered horizontally/vertically in cells; 
# grey font color = rgb 51/51/51, white for column headers
# flag column is width = 15
# country name column is width = 32
# data columns are width = 28
# column header row height = 95
# data row height = 35
# data cells formatted as number with two decimal places, to avoid truncation (e.g. 0.70 converted to 0.7)

# remember to fully spell out macedonia and bih
# output is copied into powerpoint using paste special -> Picture (Enhanced Metafile), aligned center
# flags are copy/pasted in from state dept website, w locked aspect ratio, and resized to 70%

# save
# fmir_obj_4_indicator_similarity_table %>% write_csv(file = "output/charts/fmir_obj_4_indicator_similarity_table.csv")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# get cluster tendency to see if data is likely to have clusters that are more meaningful/tight than uniform random sample
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# note that the hopkins stat uses .5 as a threshold, with higher values indicating there is a clustering tendency; 
# a hopkins stat over .75 is 90% confidence that there is clustering tendency that is real, and not random
cluster_tendency <- get_clust_tendency(fmir_obj_4_indicator, n = nrow(fmir_obj_4_indicator) - 1, graph = FALSE)


#///////////////////


# inspect
cluster_tendency
cluster_tendency$hopkins_stat


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters
fmir_obj_4_indicator_hclust <- hclust(d = fmir_obj_4_indicator_distance_matrix, method = "complete" )


#/////////////////////


# inspect
fmir_obj_4_indicator_hclust
fmir_obj_4_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_4_indicator_hclust$merge
fmir_obj_4_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_4_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters with agnes, which allows enhanced fviz visualizations
fmir_obj_4_indicator_agnes <- agnes(x = fmir_obj_4_indicator, method = "complete")


#//////////////////


# inspect
fmir_obj_4_indicator_agnes
fmir_obj_4_indicator_agnes %>% attributes()

# get agglomerative coefficient from agnes object 
# "Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 suggest a 
# more balanced clustering structure such as the complete linkage and Wards method dendrograms in Figure 21.3. 
# Values closer to 0 suggest less well-formed clusters such as the single linkage dendrogram in Figure 21.3. 
# However, the AC tends to become larger as n increases, so it should not be used to compare across data sets of very different sizes."
fmir_obj_4_indicator_agnes$ac

# plot
pltree(fmir_obj_4_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 



#//////////////////////////////////////////////////////////////////////////////////////////////////


# use agnes to look for clustering method with best agglomerative coefficient 
method_list <- c( "average", "single", "complete", "ward")
names(method_list) <- c( "average", "single", "complete", "ward")
method_list

# function to compute coefficient
cluster_methods <- function(current_method) {
        print(current_method)
        agnes(fmir_obj_4_indicator_distance_matrix, method = current_method)$ac
}

# run cluster_methods
# note that wards method has the highest/best agglomerative coefficient
map(.x = method_list, .f = ~ cluster_methods(current_method = .x))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with agnes, which had highest (best) agglomerative coefficient
fmir_obj_4_indicator_agnes <- agnes(fmir_obj_4_indicator, method = "ward")
fmir_obj_4_indicator_agnes


#///////////////////


# visualize clustering 
pltree(fmir_obj_4_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
fviz_dend(fmir_obj_4_indicator_agnes, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

# inspect optimal number of clusters
# note the gap statistic is recently developed by tibshirani et al
fviz_nbclust(fmir_obj_4_indicator, FUN = hcut, method = "wss")
fviz_nbclust(fmir_obj_4_indicator, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(fmir_obj_4_indicator, FUN = hcut, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# set optimal k
fmir_obj_4_indicator_k <- 2


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with hclust
fmir_obj_4_indicator_hclust <- hclust(d = fmir_obj_4_indicator_distance_matrix, method = "ward.D2")


#/////////////////////


# inspect
fmir_obj_4_indicator_hclust
fmir_obj_4_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_4_indicator_hclust$merge
fmir_obj_4_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_4_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut tree 

# build and cut agnes tree with hcut for enhanced fviz
fmir_obj_4_indicator_agnes_cut <- hcut(fmir_obj_4_indicator, k = fmir_obj_4_indicator_k, hc_method = "ward.D2", stand = TRUE)
fmir_obj_4_indicator_agnes_cut 

# view PCA for clusters, which requires an hcut object
fviz_cluster(object = fmir_obj_4_indicator_agnes_cut)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut hclust tree with cutree
fmir_obj_4_indicator_hclust_cut <- cutree(fmir_obj_4_indicator_hclust, k = fmir_obj_4_indicator_k)
fmir_obj_4_indicator_hclust_cut

# Number of members in each cluster
table(fmir_obj_4_indicator_hclust_cut)

# can add cutree output to df
fmir_obj_4_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value) %>%
        left_join(fmir_obj_4_indicator %>% rownames_to_column(var = "country"), ., by = "country")

# plot dendrogram with border around 4 clusters
plot(fmir_obj_4_indicator_hclust, cex = 0.6, hang = -1)
rect.hclust(fmir_obj_4_indicator_hclust, k = fmir_obj_4_indicator_k, border = 2:5)

# plot dendrogram w labels
dendrogram_obj <- as.dendrogram(fmir_obj_4_indicator_hclust)
dendrogram_obj %>% unclass() %>% str()
dendrogram_obj %>% unclass() %>% class()

subdendrogram_list <- get_subdendrograms(dendrogram_obj, k = fmir_obj_4_indicator_k)
subdendrogram_list
dendrogram_labels <- map(.x = subdendrogram_list, .f = ~ labels(.x)) %>% 
        enframe() %>% unnest(cols = value) %>% 
        rename(subdendrogram_plot_order = name, country = value) %>%
        left_join(fmir_obj_4_indicator %>% mutate(cluster = fmir_obj_4_indicator_hclust_cut) %>%
                          rownames_to_column(var = "country"),
                  ., by = "country") %>% 
        distinct(cluster, subdendrogram_plot_order) %>% 
        arrange(subdendrogram_plot_order) %>%
        mutate(cluster = str_c("\n\n\n", cluster)) %>%
        pull(cluster) %>% unname()
dendrogram_labels        

plot(fmir_obj_4_indicator_hclust, hang = -1)
rect.dendrogram(as.dendrogram(fmir_obj_4_indicator_hclust), k = fmir_obj_4_indicator_k, text = dendrogram_labels)

# visualize clusters along 2 primary PCA dimensions
fviz_cluster(list(data = fmir_obj_4_indicator, cluster = fmir_obj_4_indicator_hclust_cut))


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_4_indicator_dendrogram ####

# get fmir_obj_4_indicator_dendrogram_data
# note that dendro_data_k() takes the uncut hclust as first argument
fmir_obj_4_indicator_dendrogram_data <- dendro_data_k(hc = fmir_obj_4_indicator_hclust, 
                                                      k = fmir_obj_4_indicator_k)

# inspect
fmir_obj_4_indicator_dendrogram_data
fmir_obj_4_indicator_dendrogram_data %>% attributes()
fmir_obj_4_indicator_dendrogram_data$segments %>% head()
fmir_obj_4_indicator_dendrogram_data$labels %>% head()
fmir_obj_4_indicator_dendrogram_data$leaf_labels %>% head()
fmir_obj_4_indicator_dendrogram_data$class %>% head()


#///////////////////


# get cluster_colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
cluster_colors <- c(color_palette %>% slice(2) %>% pull(hex), 
                    color_palette %>% slice(7) %>% pull(hex))
cluster_colors


# plot
fmir_obj_4_indicator_dendrogram <- plot_ggdendro(fmir_obj_4_indicator_dendrogram_data,
                                                 direction = "tb",
                                                 scale.color = cluster_colors,
                                                 label.size = 2.5,
                                                 branch.size = 0.5,
                                                 expand.y = 0) +
        labs(y = "Similarity\n(lower groupings are more similar)", x = "") +
        coord_fixed(ratio = 1 / 1.25, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 5, 0, 0), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(size = .5),
                axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "bold", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 0, hjust = .5),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                # axis.title.x = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 5, r = 0, b = 10, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 0, r = 15, b = 0, l = 5)),
                # axis.title.y = element_blank(),
                # axis.text.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 9, family = "Calibri", face = "plain", color = "#333333"),
                legend.text = element_text(size = 9, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# guides(color = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333", keywidth = 4),
#        linetype = guide_legend(keywidth = 4))

# inspect
fmir_obj_4_indicator_dendrogram


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_4_indicator_dendrogram)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_4_indicator_dendrogram.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_4_country_indicator_tile_chart ####

# get chart_data
# note this is 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        left_join(., fmir_obj_1_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        left_join(., fmir_obj_1_country_order_tbl, by = "country") %>%
        relocate(country_order, .after = cluster) %>% arrange(country_order) %>%
        pivot_wider(id_cols = c(country, cluster, country_order), 
                    names_from = indicator_name, values_from = indicator_standardized_values) %>%
        left_join(., fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          distinct(country, obj_num, sub_obj_num, sub_obj_avg, obj_avg) %>%
                          mutate(sub_obj_num_label = case_when(sub_obj_num == "sub_obj_4_1" ~ 
                                                                       "Sub-obj. 4.1: Economic competitiveness",
                                                               sub_obj_num == "sub_obj_4_2" ~ "Sub-obj. 4.2: Financial health")) %>%
                          select(country, obj_num, obj_avg, sub_obj_num_label, sub_obj_avg) %>%
                          pivot_wider(id_cols = c(country, obj_num, obj_avg), 
                                      names_from = sub_obj_num_label, values_from = sub_obj_avg) %>%
                          rename(!!sym("Obj. 4: Economic") := obj_avg) %>% select(-obj_num), 
                  by = "country") %>%
        pivot_longer(cols = -c(country, cluster, country_order), names_to = "var", values_to = "values") %>%
        mutate(color = case_when(cluster == 1 ~ color_palette %>% slice(2) %>% pull(hex),
                                 cluster == 2 ~ color_palette %>% slice(7) %>% pull(hex),
                                 cluster == 3 ~ color_palette %>% slice(10) %>% pull(hex)),
               country_label = str_c("Group ", cluster, ": ", country),
               var_order = case_when(var == "Obj. 4: Economic" ~ 1,
                                     str_detect(string = var, pattern = "^Sub-obj.") ~ 2,
                                     TRUE ~ 3))

# inspect
chart_data


#/////////////////////////////


# create chart
fmir_obj_4_country_indicator_tile_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       y = fct_reorder(.f = country_label, .x = country_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Score\n(higher is better)\n", y = NULL, x = "Indicator") +
        coord_fixed(ratio = 1 / 5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 5, 
                                           color = chart_data %>% distinct(country, color) %>% pull(color), 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 5, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 5, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_4_country_indicator_tile_chart <- ggdraw(align_legend(fmir_obj_4_country_indicator_tile_chart))


# inspect
fmir_obj_4_country_indicator_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_4_country_indicator_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_4_country_indicator_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_4_indicator_cluster_boxplot_chart ####

# get chart_data
# note that this uses the original ISV values from the index, not the scaled ISV used in fmir_obj_4_indicator for clustering; 
# original/unscaled ISV avg is better for plot because values are then interpretable in the framework of the index (e.g. ISV btw 0-1)
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_4_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        
        
        # add cluster avg of sub_obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          left_join(., fmir_obj_4_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_4_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -c(country, cluster), names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster),
               color_bin = as.character(cluster_name),
               color = case_when(color_bin == "Group 1" ~ color_palette %>% slice(2) %>% pull(hex),
                                 color_bin == "Group 2" ~ color_palette %>% slice(7) %>% pull(hex),
                                 color_bin == "Group 3" ~ color_palette %>% slice(10) %>% pull(hex)))

# inspect
chart_data
chart_data %>% glimpse()

# create color_list for to pass to scale_color_manual
chart_data_color_list <- chart_data %>% count(color_bin, color) %>% pull(color)
names(chart_data_color_list) <- chart_data %>% count(color_bin, color) %>% pull(color_bin)
chart_data_color_list


#//////////////////////////


# create chart
fmir_obj_4_indicator_cluster_boxplot_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = values, fill = color_bin)) + 
        geom_boxplot(width = .5, lwd = .01, outlier.size = .5) + facet_wrap(facets = vars(cluster_name)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = .2), limits = c(0, 1), expand = c(0, 0)) +
        scale_fill_manual(values = chart_data_color_list) + 
        # scale_color_manual(values = chart_data_color_list) + 
        coord_fixed(ratio = 10/1.25, clip = "off") +
        labs(x = "Indicator", y = "Score (higher is better)", fill = NULL, color = NULL) +
        theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                # strip.background = element_blank(),
                strip.background = element_rect(fill = "#DDDDDD", color = NA),
                strip.text = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333"),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 5, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = -15, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 5, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = 1,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 


# inspect
fmir_obj_4_indicator_cluster_boxplot_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_4_indicator_cluster_boxplot_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_4_indicator_cluster_boxplot_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_4_indicator_cluster_tile_chart ####

# get chart_data 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_4_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        group_by(cluster) %>%
        mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
        ungroup() %>%
        select(-country) %>% distinct() %>%
        
        
        # add cluster avg of sub_obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          left_join(., fmir_obj_4_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_4_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -cluster, names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster))

# inspect
chart_data


#///////////////


# get obj/sub_obj segment coordinates

# note that in the entire fmir dataset, the maximum sub_obj per obj is 3, so will have 3 segment_tbls,
# and for those obj with less than 3 sub_obj, the segment_tbl will just have NA coordinates and no output in ggplot
fmir %>% distinct(obj_num, sub_obj_num) %>% 
        add_count(obj_num, name = "sub_obj_count") %>% 
        distinct(obj_num, sub_obj_count) %>%
        arrange(desc(sub_obj_count)) 

# segment_1_tbl
# segment_1_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(obj_num, sub_obj_num) %>%
#         mutate(nrow = nrow(.),
#                x = nrow + 1.5,
#                xend = nrow + 1.5,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5) %>%
#         select(x, xend, y, yend) %>%
#         slice(1)
# segment_1_tbl

# # segment_2_tbl
# segment_2_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                 cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(1)
# segment_2_tbl
#                
# # segment_3_tbl
# segment_3_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(2)
# segment_3_tbl


#///////////////////////


# plot

# get order_obj_and_sub_obj_avg_first()
order_obj_and_sub_obj_avg_first <- function(var){
        
        tibble(var = var) %>% 
                left_join(., fmir %>% distinct(indicator_name, concept) %>% 
                                  add_group_index(group_vars = concept, group_name = "concept_index"), 
                          by = c("var" = "indicator_name")) %>%
                mutate(order = case_when(str_detect(string = var, pattern = "^obj_") ~ 1,
                                         str_detect(string = var, pattern = "^sub_obj_[0-9]_[0-9]$") ~ 2,
                                         TRUE ~ 3)) %>%
                arrange(order, concept_index, var) %>% pull(var)
}

# create chart
fmir_obj_4_indicator_cluster_tile_chart <- chart_data %>% 
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = factor(cluster_name), fill = values, label = round_to_digits(values, digits = 2))) + 
        geom_tile(color = "#DDDDDD") + 
        geom_text(color = "#ffffff", size = 1.75, family = "Calibri") +
        scale_y_discrete(expand = c(0, 0)) +
        # geom_segment(x = segment_1_tbl$x, xend = segment_1_tbl$xend, 
        #              y = segment_1_tbl$y, yend = segment_1_tbl$yend, color = "#ffffff", size = 3) +
        # geom_segment(x = segment_2_tbl$x, xend = segment_2_tbl$xend,
        #              y = segment_2_tbl$y, yend = segment_2_tbl$yend, color = "#CBCBCB", size = 3) +
        # geom_segment(x = segment_3_tbl$x, xend = segment_3_tbl$xend,
        #              y = segment_3_tbl$y, yend = segment_3_tbl$yend, color = "#CBCBCB", size = 3) +
        labs(x = "Indicator", y = NULL, fill = "Avg. score\n(higher is better)\n") +
        coord_fixed(ratio = 2 / 1, clip = "off") +
        theme_bw() +
        # theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD"),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "right",
                legend.key.size = unit(4, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 5, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# fmir_obj_4_indicator_cluster_tile_chart <- ggdraw(align_legend(fmir_obj_4_indicator_cluster_tile_chart))

# inspect
fmir_obj_4_indicator_cluster_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_4_indicator_cluster_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_4_indicator_cluster_tile_chart.docx")



#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_c_indicator_clusters ####

# set current_obj_num
current_obj_num <- "obj_c"

# get fmir_obj_c_indicator with scaled version of indicator_standardized_values
fmir_obj_c_indicator <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        distinct(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        column_to_rownames(var = "country") %>%
        mutate(across(.cols = everything(), .fns = ~ scale(.x, center = TRUE, scale = TRUE)[ , 1]))


#//////////////////////


# inspect
fmir_obj_c_indicator
fmir_obj_c_indicator %>% glimpse()
fmir_obj_c_indicator %>% nrow()
fmir_obj_c_indicator %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get fmir_obj_c_indicator_distance_matrix 
fmir_obj_c_indicator_distance_matrix <- dist(fmir_obj_c_indicator, method = "euclidean")
# fmir_obj_c_indicator_distance_matrix <- get_dist(fmir_obj_c_indicator, method = "euclidean")


#//////////////


# inspect
fmir_obj_c_indicator_distance_matrix
fmir_obj_c_indicator_distance_matrix %>% str()
fmir_obj_c_indicator_distance_matrix %>% attributes()

# inspect fmir_obj_c_indicator_distance_matrix_tbl
fmir_obj_c_indicator_distance_matrix_tbl <- as.matrix(fmir_obj_c_indicator_distance_matrix) %>% as_tibble() %>% 
        mutate(country = row.names(as.matrix(fmir_obj_c_indicator_distance_matrix))) %>%
        relocate(country, .before = everything())

fmir_obj_c_indicator_distance_matrix_tbl
fmir_obj_c_indicator_distance_matrix_tbl %>% skim() %>% as_tibble() %>% 
        select(skim_variable, numeric.mean, starts_with(match = "numeric.p")) %>%
        # arrange(numeric.mean) %>% 
        arrange(numeric.p25) %>%
        print(n = nrow(.))
fmir_obj_c_indicator_distance_matrix_tbl %>% 
        # select(country, Serbia) %>% arrange(Serbia) %>% print(n = nrow(.))
        # select(country, BiH) %>% arrange(BiH) %>% print(n = nrow(.))
        # select(country, Kosovo) %>% arrange(Kosovo) %>% print(n = nrow(.))
        # select(country, Armenia) %>% arrange(Armenia) %>% print(n = nrow(.))
        select(country, Poland) %>% arrange(Poland) %>% print(n = nrow(.))
# filter(country == "Serbia") %>% select(ends_with("stan"))

# compare indicator_standardized_values
fmir %>% filter(obj_num == "obj_3", country %in% c("Poland", "Kazakhstan"), year == 2020) %>% 
        select(country, year, indicator_name, values, values_z_std, indicator_standardized_values) %>%
        arrange(country) %>% 
        # unite(col = var, indicator_name, country) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = indicator_name, names_from = country, values_from = indicator_standardized_values)

# visualize distance
viridis_pal()(3) %>% show_col()
fviz_dist(dist.obj = fmir_obj_c_indicator_distance_matrix, gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_c_indicator_distance_matrix_tile_chart ####


# get country_order from fviz_dist_plot
fmir_obj_c_fviz_dist_plot <- fviz_dist(dist.obj = fmir_obj_c_indicator_distance_matrix, 
                                       gradient = list(low = "#440154FF", mid = "#21908CFF", high = "#FDE725FF"))
fmir_obj_c_fviz_dist_plot
fmir_obj_c_fviz_dist_plot %>% attributes()
fmir_obj_c_fviz_dist_plot$data %>% as_tibble()

fmir_obj_c_country_order_tbl <- fmir_obj_c_fviz_dist_plot$data %>% as_tibble() %>% slice(1:(fmir_obj_c_indicator %>% nrow())) %>% 
        mutate(country = str_sub(string = Var1, start = 1, end = -2),
               country_order = row_number()) %>% 
        select(country, country_order)
fmir_obj_c_country_order_tbl


#/////////////////////////////


# create chart
fmir_obj_c_indicator_distance_matrix_tile_chart <- fmir_obj_c_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "var", values_to = "values") %>%
        left_join(., fmir_obj_c_country_order_tbl, by = c("var" = "country")) %>%
        rename(var_order = country_order) %>%
        left_join(., fmir_obj_c_country_order_tbl, by = "country") %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = country, .x = country_order, .fun = min), 
                                       y = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Similiarity\n(lower is more similar)", y = NULL, x = NULL) +
        coord_fixed(ratio = 1/1.5, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 13, r = 0, b = 5, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_c_indicator_distance_matrix_tile_chart <- ggdraw(align_legend(fmir_obj_c_indicator_distance_matrix_tile_chart))


# inspect
fmir_obj_c_indicator_distance_matrix_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_c_indicator_distance_matrix_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_c_indicator_distance_matrix_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_c_indicator_similarity_table ####

fmir_obj_c_indicator_similarity_table <- fmir_obj_c_indicator_distance_matrix_tbl %>% 
        pivot_longer(cols = -country, names_to = "counterpart", values_to = "distance") %>%
        filter(country != counterpart) %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping), by = "country") %>%
        left_join(., country_crosswalk %>% select(country, mcp_grouping) %>% rename(counterpart_mcp_grouping = mcp_grouping), 
                  by = c("counterpart" = "country")) %>%
        group_by(country) %>%
        mutate(distance_in_ee = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_flag = case_when(distance == min(distance) ~ 1,
                                               TRUE ~ 0),
               max_similarity = case_when(max_similarity_flag == 1 ~ distance,
                                          TRUE ~ NA_real_),
               max_similarity_counterpart = case_when(max_similarity_flag == 1 ~ counterpart,
                                                      TRUE ~ NA_character_),
               max_similarity_in_ee_flag = case_when(counterpart_mcp_grouping %in% c("E&E Balkans", "E&E Eurasia") &
                                                             distance_in_ee == min(distance_in_ee, na.rm = TRUE) ~ 1,
                                                     TRUE ~ NA_real_),
               max_similarity_in_ee = case_when(max_similarity_in_ee_flag == 1 ~ distance_in_ee,
                                                TRUE ~ NA_real_),
               max_similarity_in_ee_counterpart = case_when(max_similarity_in_ee_flag == 1 ~ counterpart,
                                                            TRUE ~ NA_character_)) %>%
        fill(max_similarity, .direction = "updown") %>%
        fill(max_similarity_counterpart, .direction = "updown") %>%
        fill(max_similarity_in_ee, .direction = "updown") %>%
        fill(max_similarity_in_ee_counterpart, .direction = "updown") %>%
        ungroup()


#//////////////////


# inspect
fmir_obj_c_indicator_similarity_table
fmir_obj_c_indicator_similarity_table %>% glimpse()
fmir_obj_c_indicator_similarity_table %>% nrow() # 756
fmir_obj_c_indicator_similarity_table %>% ncol() # 12


fmir_obj_c_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping)) %>% 
        print(n = 30)
fmir_obj_c_indicator_similarity_table %>% 
        select(-c(mcp_grouping, counterpart_mcp_grouping, max_similarity, max_similarity_counterpart)) %>% 
        print(n = 30)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get distinct max_similarity/_in_ee values/counterparts and text_output for table
fmir_obj_c_indicator_similarity_table <- fmir_obj_c_indicator_similarity_table %>% 
        distinct(country, mcp_grouping, max_similarity, max_similarity_counterpart,
                 max_similarity_in_ee, max_similarity_in_ee_counterpart) %>%
        mutate(max_similarity_text_output = str_c("   ", max_similarity_counterpart, " (", round_to_digits(max_similarity, digits = 2), ")"),
               max_similarity_in_ee_text_output = str_c("   ", max_similarity_in_ee_counterpart, " (", 
                                                        round_to_digits(max_similarity_in_ee, digits = 2), ")"),
               mcp_grouping_order = case_when(mcp_grouping == "E&E Balkans" ~ 1,
                                              mcp_grouping == "E&E Eurasia" ~ 2,
                                              mcp_grouping == "CARs" ~ 3,
                                              mcp_grouping == "E&E graduates" ~ 4)) %>%
        arrange(mcp_grouping_order, country) %>%
        select(country, mcp_grouping, 
               max_similarity_text_output, max_similarity, 
               max_similarity_in_ee_text_output, max_similarity_in_ee)


#///////////////////


# inspect
fmir_obj_c_indicator_similarity_table
fmir_obj_c_indicator_similarity_table %>% glimpse()
fmir_obj_c_indicator_similarity_table %>% nrow() # 28
fmir_obj_c_indicator_similarity_table %>% ncol() # 5

fmir_obj_c_indicator_similarity_table %>% print(n = nrow(.))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get table_palette


# use viridis palette
viridis_pal()(10) %>% show_col()


# historic
# note viridis is not used to avoid confusing color scales when compared to dissimilarity matrix, which has much higher
# range of values than the similarity_table
# table_palette <- colorRampPalette(c(brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(4) %>% pull(value), 
#                                      brewer.pal(n = 9, "Blues") %>% as_tibble() %>% slice(9) %>% pull(value)))


#////////////////////////


# inspect
viridis_pal()(10) %>% show_col()
brewer.pal(n = 9, "Blues") %>% show_col()
table_palette(5) %>% show_col()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get values_normalized_quantile_cutoffs, creating one set of color thresholds to be applied uniformly to all obj in the table
# this just makes it easier to automatically assign color-coded quantiles for the full range of values for E&E countries
# whereas using the actual values values themselves would require checks after updates to make sure the cutoffs are still balanced
# across the full range of values values
# note that this doesn't include miri_avg, but the table does, but this is fine because miri_avg by definition must fall within
# range of max/min values, and since this only gets the color quantiles based on max/min values, it's fine to exclude miri_avg
values_normalized_quantile_cutoffs <- fmir_obj_c_indicator_similarity_table %>%
        pivot_longer(cols = c(max_similarity, max_similarity_in_ee), names_to = "var", values_to = "values") %>%
        mutate(values_min = min(values),
               values_max = max(values),
               values_normalized = (values - values_min) / (values_max - values_min),
               q0_to_q10 = .10 * (values_max - values_min) + values_min,
               q11_to_q20 = .20 * (values_max - values_min) + values_min,
               q21_to_q30 = .30 * (values_max - values_min) + values_min,
               q31_to_q40 = .40 * (values_max - values_min) + values_min,
               q41_to_q50 = .50 * (values_max - values_min) + values_min,
               q51_to_q60 = .60 * (values_max - values_min) + values_min,
               q61_to_q70 = .70 * (values_max - values_min) + values_min,
               q71_to_q80 = .80 * (values_max - values_min) + values_min,
               q81_to_q90 = .90 * (values_max - values_min) + values_min,
               q91_to_q100 = 1 * (values_max - values_min) + values_min) %>%
        select(-c(max_similarity_text_output, max_similarity_in_ee_text_output))


#/////////////////////////


# inspect
values_normalized_quantile_cutoffs %>% data.frame()
# values_normalized = (values - values_min) / (values_max - values_min)

# check q0-q10 threshold
(2.25 - 1.81) / (6.18 - 1.81)

# check q81-q90 threshold
(5.74 - 1.81) / (6.18 - 1.81)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# create table_palette_tbl
# manually apply conditional formatting in excel to highlight and change font color
# =IF(AND(B4>5, B4<10),TRUE,FALSE)
table_palette_tbl <- values_normalized_quantile_cutoffs %>% select(starts_with("q")) %>% distinct() %>%
        pivot_longer(cols = everything(), names_to = "quantile", values_to = "values_threshold") %>%
        left_join(., tibble(hex = viridis_pal()(10), 
                            quantile = c("q0_to_q10", "q11_to_q20", "q21_to_q30", "q31_to_q40", "q41_to_q50",
                                         "q51_to_q60", "q61_to_q70", "q71_to_q80", "q81_to_q90", "q91_to_q100")) %>%
                          mutate(rgb = map(.x = hex, .f = ~ col2rgb(col = .x) %>% str_c(., collapse = "."))) %>%
                          unnest(rgb),
                  by = "quantile")

# inspect
table_palette_tbl 


#//////////////////////////////////////////////////////////////////////////////////////////////////


# write to file

# note csv will be copied into xlsx shell and formatted
# old: blue header is rgb 31/73/125
# current: gray header is rgb 128/128/128
# grey grid lines are rgb 217/217/217, medium weight
# country rows on left side gets a rgb 217/217/217 heavy border
# grey top/bottom row above/under table is rgb 217/217/217, height = 15
# font is Calibri, size 15, bolded, centered horizontally/vertically in cells; 
# grey font color = rgb 51/51/51, white for column headers
# flag column is width = 15
# country name column is width = 32
# data columns are width = 28
# column header row height = 95
# data row height = 35
# data cells formatted as number with two decimal places, to avoid truncation (e.g. 0.70 converted to 0.7)

# remember to fully spell out macedonia and bih
# output is copied into powerpoint using paste special -> Picture (Enhanced Metafile), aligned center
# flags are copy/pasted in from state dept website, w locked aspect ratio, and resized to 70%

# save
# fmir_obj_c_indicator_similarity_table %>% write_csv(file = "output/charts/fmir_obj_c_indicator_similarity_table.csv")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# get cluster tendency to see if data is likely to have clusters that are more meaningful/tight than uniform random sample
# https://www.datanovia.com/en/lessons/assessing-clustering-tendency/
# note that the hopkins stat uses .5 as a threshold, with higher values indicating there is a clustering tendency; 
# a hopkins stat over .75 is 90% confidence that there is clustering tendency that is real, and not random
cluster_tendency <- get_clust_tendency(fmir_obj_c_indicator, n = nrow(fmir_obj_c_indicator) - 1, graph = FALSE)


#///////////////////


# inspect
cluster_tendency
cluster_tendency$hopkins_stat


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters
fmir_obj_c_indicator_hclust <- hclust(d = fmir_obj_c_indicator_distance_matrix, method = "complete" )


#/////////////////////


# inspect
fmir_obj_c_indicator_hclust
fmir_obj_c_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_c_indicator_hclust$merge
fmir_obj_c_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_c_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get hierarchical clusters with agnes, which allows enhanced fviz visualizations
fmir_obj_c_indicator_agnes <- agnes(x = fmir_obj_c_indicator, method = "complete")


#//////////////////


# inspect
fmir_obj_c_indicator_agnes
fmir_obj_c_indicator_agnes %>% attributes()

# get agglomerative coefficient from agnes object 
# "Generally speaking, the AC describes the strength of the clustering structure. Values closer to 1 suggest a 
# more balanced clustering structure such as the complete linkage and Wards method dendrograms in Figure 21.3. 
# Values closer to 0 suggest less well-formed clusters such as the single linkage dendrogram in Figure 21.3. 
# However, the AC tends to become larger as n increases, so it should not be used to compare across data sets of very different sizes."
fmir_obj_c_indicator_agnes$ac

# plot
pltree(fmir_obj_c_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 



#//////////////////////////////////////////////////////////////////////////////////////////////////


# use agnes to look for clustering method with best agglomerative coefficient 
method_list <- c( "average", "single", "complete", "ward")
names(method_list) <- c( "average", "single", "complete", "ward")
method_list

# function to compute coefficient
cluster_methods <- function(current_method) {
        print(current_method)
        agnes(fmir_obj_c_indicator_distance_matrix, method = current_method)$ac
}

# run cluster_methods
# note that wards method has the highest/best agglomerative coefficient
map(.x = method_list, .f = ~ cluster_methods(current_method = .x))


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with agnes, which had highest (best) agglomerative coefficient
fmir_obj_c_indicator_agnes <- agnes(fmir_obj_c_indicator, method = "ward")
fmir_obj_c_indicator_agnes


#///////////////////


# visualize clustering 
pltree(fmir_obj_c_indicator_agnes, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
fviz_dend(fmir_obj_c_indicator_agnes, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07"))

# inspect optimal number of clusters
# note the gap statistic is recently developed by tibshirani et al
fviz_nbclust(fmir_obj_c_indicator, FUN = hcut, method = "wss")
fviz_nbclust(fmir_obj_c_indicator, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(fmir_obj_c_indicator, FUN = hcut, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# set optimal k
fmir_obj_c_indicator_k <- 3


#//////////////////////////////////////////////////////////////////////////////////////////////////


# use ward's method with hclust
fmir_obj_c_indicator_hclust <- hclust(d = fmir_obj_c_indicator_distance_matrix, method = "ward.D2")


#/////////////////////


# inspect
fmir_obj_c_indicator_hclust
fmir_obj_c_indicator_hclust %>% attributes()

# inspect hclust object
?hclust
fmir_obj_c_indicator_hclust$merge
fmir_obj_c_indicator_hclust$height

# plot dendrogram
plot(fmir_obj_c_indicator_hclust, cex = 0.6, hang = -1)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut tree 

# build and cut agnes tree with hcut for enhanced fviz
fmir_obj_c_indicator_agnes_cut <- hcut(fmir_obj_c_indicator, k = fmir_obj_c_indicator_k, hc_method = "ward.D2", stand = TRUE)
fmir_obj_c_indicator_agnes_cut 

# view PCA for clusters, which requires an hcut object
fviz_cluster(object = fmir_obj_c_indicator_agnes_cut)


#//////////////////////////////////////////////////////////////////////////////////////////////////


# cut hclust tree with cutree
fmir_obj_c_indicator_hclust_cut <- cutree(fmir_obj_c_indicator_hclust, k = fmir_obj_c_indicator_k)
fmir_obj_c_indicator_hclust_cut

# Number of members in each cluster
table(fmir_obj_c_indicator_hclust_cut)

# can add cutree output to df
fmir_obj_c_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value) %>%
        left_join(fmir_obj_c_indicator %>% rownames_to_column(var = "country"), ., by = "country")

# plot dendrogram with border around 4 clusters
plot(fmir_obj_c_indicator_hclust, cex = 0.6, hang = -1)
rect.hclust(fmir_obj_c_indicator_hclust, k = fmir_obj_c_indicator_k, border = 2:5)

# plot dendrogram w labels
dendrogram_obj <- as.dendrogram(fmir_obj_c_indicator_hclust)
dendrogram_obj %>% unclass() %>% str()
dendrogram_obj %>% unclass() %>% class()

subdendrogram_list <- get_subdendrograms(dendrogram_obj, k = fmir_obj_c_indicator_k)
subdendrogram_list
dendrogram_labels <- map(.x = subdendrogram_list, .f = ~ labels(.x)) %>% 
        enframe() %>% unnest(cols = value) %>% 
        rename(subdendrogram_plot_order = name, country = value) %>%
        left_join(fmir_obj_c_indicator %>% mutate(cluster = fmir_obj_c_indicator_hclust_cut) %>%
                          rownames_to_column(var = "country"),
                  ., by = "country") %>% 
        distinct(cluster, subdendrogram_plot_order) %>% 
        arrange(subdendrogram_plot_order) %>%
        mutate(cluster = str_c("\n\n\n", cluster)) %>%
        pull(cluster) %>% unname()
dendrogram_labels        

plot(fmir_obj_c_indicator_hclust, hang = -1)
rect.dendrogram(as.dendrogram(fmir_obj_c_indicator_hclust), k = fmir_obj_c_indicator_k, text = dendrogram_labels)

# visualize clusters along 2 primary PCA dimensions
fviz_cluster(list(data = fmir_obj_c_indicator, cluster = fmir_obj_c_indicator_hclust_cut))


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_c_indicator_dendrogram ####

# get fmir_obj_c_indicator_dendrogram_data
# note that dendro_data_k() takes the uncut hclust as first argument
fmir_obj_c_indicator_dendrogram_data <- dendro_data_k(hc = fmir_obj_c_indicator_hclust, 
                                                      k = fmir_obj_c_indicator_k)

# inspect
fmir_obj_c_indicator_dendrogram_data
fmir_obj_c_indicator_dendrogram_data %>% attributes()
fmir_obj_c_indicator_dendrogram_data$segments %>% head()
fmir_obj_c_indicator_dendrogram_data$labels %>% head()
fmir_obj_c_indicator_dendrogram_data$leaf_labels %>% head()
fmir_obj_c_indicator_dendrogram_data$class %>% head()


#///////////////////


# get cluster_colors
# note that colors are passed in a vector that is k + 1 in length
# the first value is the neutral color for branches above the point at which the tree has been cut up to the root
# the second/third/fourth etc values are indexed to the cluster number; 
# so cluster 1 is second value in colors vector, cluster 2 is third value in colors vector, cluster 3 is fourth value in colors vector
cluster_colors <- c("#333333", 
                    color_palette %>% slice(2) %>% pull(hex), 
                    color_palette %>% slice(7) %>% pull(hex),
                    color_palette %>% slice(10) %>% pull(hex))
cluster_colors


# plot
fmir_obj_c_indicator_dendrogram <- plot_ggdendro(fmir_obj_c_indicator_dendrogram_data,
                                                 direction = "tb",
                                                 scale.color = cluster_colors,
                                                 label.size = 2.5,
                                                 branch.size = 0.5,
                                                 expand.y = 0) +
        labs(y = "Similarity\n(lower groupings are more similar)", x = "") +
        coord_fixed(ratio = 1 / .75, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 5, 0, 0), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(size = .5),
                axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "bold", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 0, hjust = .5),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                # axis.title.x = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 5, r = 0, b = 10, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333",
                                            margin = margin(t = 0, r = 15, b = 0, l = 5)),
                # axis.title.y = element_blank(),
                # axis.text.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 9, family = "Calibri", face = "plain", color = "#333333"),
                legend.text = element_text(size = 9, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# guides(color = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333", keywidth = 4),
#        linetype = guide_legend(keywidth = 4))

# inspect
fmir_obj_c_indicator_dendrogram


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_c_indicator_dendrogram)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_c_indicator_dendrogram.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_c_country_indicator_tile_chart ####

# get chart_data
# note this is 
chart_data <- fmir %>% filter(year == 2020) %>%
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        left_join(., fmir_obj_c_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>%
        relocate(cluster, .after = country) %>%
        left_join(., fmir_obj_c_country_order_tbl, by = "country") %>%
        relocate(country_order, .after = cluster) %>% arrange(country_order) %>%
        pivot_wider(id_cols = c(country, cluster, country_order),
                    names_from = indicator_name, values_from = indicator_standardized_values) %>%
        left_join(., fmir %>%
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          distinct(country, obj_num, sub_obj_num, sub_obj_avg, obj_avg) %>%
                          pivot_wider(id_cols = c(country, obj_num, obj_avg), names_from = sub_obj_num, values_from = sub_obj_avg) %>%
                          rename(!!sym("Cross-cutting obj: Corruption") := obj_avg) %>% select(-obj_num) %>% select(-sub_obj_c),
                  by = "country") %>%
        pivot_longer(cols = -c(country, cluster, country_order), names_to = "var", values_to = "values") %>%
        mutate(color = case_when(cluster == 1 ~ color_palette %>% slice(2) %>% pull(hex),
                                 cluster == 2 ~ color_palette %>% slice(7) %>% pull(hex),
                                 cluster == 3 ~ color_palette %>% slice(10) %>% pull(hex)),
               country_label = str_c("Group ", cluster, ": ", country),
               var_order = case_when(var == "Cross-cutting obj: Corruption" ~ 1,
                                     str_detect(string = var, pattern = "^sub_obj_[0-9]_[0-9]$") ~ 2,
                                     TRUE ~ 3))

# inspect
chart_data


#/////////////////////////////


# create chart
fmir_obj_c_country_indicator_tile_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_reorder(.f = var, .x = var_order, .fun = min), 
                                       y = fct_reorder(.f = country_label, .x = country_order, .fun = min), 
                                       fill = values)) +
        geom_tile(color = "#DDDDDD") + scale_fill_viridis() +
        labs(fill = "Score\n(higher is better)\n", y = NULL, x = "Indicator") +
        coord_fixed(ratio = 1 / 12, clip = "off") +
        theme_bw() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD", size = .25),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                panel.spacing.x = unit(1, "lines"),
                panel.spacing.y = unit(1, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333"),
                # axis.ticks.y = element_blank(),
                axis.ticks.y = element_line(color = "#333333", size = .25),
                # axis.ticks.x = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                axis.ticks.length.y.left = unit(.1, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 6, color = "#333333",
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                # axis.text.x = element_blank(),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 6, 
                                           color = chart_data %>% distinct(country, color) %>% pull(color), 
                                           margin = margin(t = 0, r = 5, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 6, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 6, color = "#333333", 
                                            margin = margin(t = 0, r = 5, b = 0, l = 2)),
                plot.title = element_text(size = 7, face = "bold", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")),
                legend.position = "right",
                # legend.key.size = unit(2, "mm"), 
                legend.title = element_text(size = 6, family = "Calibri", face = "plain", color = "#333333", hjust = .5),
                legend.text = element_text(size = 6, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333"),
                # legend.spacing.x = unit(1.0, 'cm')
                # legend.spacing.y = unit(5.5, "cm")
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
                legend.key.height = unit(1, "line")
        ) 
# guides(fill = guide_legend(nrow = 2, byrow = TRUE, label.hjust = 0, color = "#333333"))
# guides(color = guide_legend(override.aes = list(size = 6)))
fmir_obj_c_country_indicator_tile_chart <- ggdraw(align_legend(fmir_obj_c_country_indicator_tile_chart))


# inspect
fmir_obj_c_country_indicator_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_c_country_indicator_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_c_country_indicator_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_c_indicator_cluster_boxplot_chart ####

# get chart_data
# note that this uses the original ISV values from the index, not the scaled ISV used in fmir_obj_c_indicator for clustering; 
# original/unscaled ISV avg is better for plot because values are then interpretable in the framework of the index (e.g. ISV btw 0-1)
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_c_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        
        
        # # add cluster avg of sub_obj_avg 
        # left_join(fmir %>% 
        #                   filter(year == 2020,
        #                          mcp_grouping != "EU-15",
        #                          obj_num == current_obj_num) %>%
        #                   select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
        #                   arrange(country) %>%
        #                   select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
        #                   distinct() %>%
        #                   pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
        #                   left_join(., fmir_obj_c_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
        #                             by = "country") %>% 
        #                   relocate(cluster, .after = country),
        #           ., by = c("country", "cluster")) %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_c_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country),
                  ., by = c("country", "cluster")) %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -c(country, cluster), names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster),
               color_bin = as.character(cluster_name),
               color = case_when(color_bin == "Group 1" ~ color_palette %>% slice(2) %>% pull(hex),
                                 color_bin == "Group 2" ~ color_palette %>% slice(7) %>% pull(hex),
                                 color_bin == "Group 3" ~ color_palette %>% slice(10) %>% pull(hex)))

# inspect
chart_data
chart_data %>% glimpse()

# create color_list for to pass to scale_color_manual
chart_data_color_list <- chart_data %>% count(color_bin, color) %>% pull(color)
names(chart_data_color_list) <- chart_data %>% count(color_bin, color) %>% pull(color_bin)
chart_data_color_list


#//////////////////////////


# create chart
fmir_obj_c_indicator_cluster_boxplot_chart <- chart_data %>%
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = values, fill = color_bin)) + 
        geom_boxplot(width = .5, lwd = .01, outlier.size = .5) + facet_wrap(facets = vars(cluster_name)) +
        scale_y_continuous(breaks = seq(from = 0, to = 1, by = .2), limits = c(0, 1), expand = c(0, 0)) +
        scale_fill_manual(values = chart_data_color_list) + 
        # scale_color_manual(values = chart_data_color_list) + 
        coord_fixed(ratio = 10/2.25, clip = "off") +
        labs(x = "Indicator", y = "Score (higher is better)", fill = NULL, color = NULL) +
        theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "#DDDDDD"),
                # panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                # strip.background = element_blank(),
                strip.background = element_rect(fill = "#DDDDDD", color = NA),
                strip.text = element_text(family = "Calibri", face = "plain", size = 8, color = "#333333"),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 5, r = 0, b = -15, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 5, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "bottom",
                legend.key.size = unit(5, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = 1,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 


# inspect
fmir_obj_c_indicator_cluster_boxplot_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_c_indicator_cluster_boxplot_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_c_indicator_cluster_boxplot_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_obj_c_indicator_cluster_tile_chart ####

# get chart_data 
chart_data <- fmir %>% filter(year == 2020) %>% 
        filter(mcp_grouping != "EU-15") %>%
        filter(obj_num == current_obj_num) %>%
        select(country, indicator_name, indicator_standardized_values) %>%
        pivot_wider(id_cols = country, names_from = indicator_name, values_from = indicator_standardized_values) %>%
        
        
        # add clusters
        left_join(., fmir_obj_c_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                  by = "country") %>% 
        relocate(cluster, .after = country) %>%
        group_by(cluster) %>%
        mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
        ungroup() %>%
        select(-country) %>% distinct() %>%
        
        
        # # add cluster avg of sub_obj_avg 
        # left_join(fmir %>% 
        #                   filter(year == 2020,
        #                          mcp_grouping != "EU-15",
        #                          obj_num == current_obj_num) %>%
        #                   select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
        #                   arrange(country) %>%
        #                   select(-c(year, obj_num, indicator_name, obj_avg)) %>% 
        #                   distinct() %>%
        #                   pivot_wider(id_cols = country, names_from = sub_obj_num, values_from = sub_obj_avg) %>%
        #                   left_join(., fmir_obj_c_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
        #                             by = "country") %>% 
        #                   relocate(cluster, .after = country) %>%
        #                   group_by(cluster) %>%
        #                   mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
        #                   ungroup() %>%
        #                   select(-country) %>% distinct(),
        #           ., by = "cluster") %>%
        
        
        # add cluster avg of obj_avg 
        left_join(fmir %>% 
                          filter(year == 2020,
                                 mcp_grouping != "EU-15",
                                 obj_num == current_obj_num) %>%
                          select(country, year, obj_num, sub_obj_num, indicator_name, sub_obj_avg, obj_avg) %>%
                          arrange(country) %>%
                          select(-c(year, sub_obj_num, indicator_name, sub_obj_avg)) %>% 
                          distinct() %>%
                          pivot_wider(id_cols = country, names_from = obj_num, values_from = obj_avg) %>%
                          left_join(., fmir_obj_c_indicator_hclust_cut %>% enframe() %>% rename(country = name, cluster = value),
                                    by = "country") %>% 
                          relocate(cluster, .after = country) %>%
                          group_by(cluster) %>%
                          mutate(across(.cols = -country, .fns = ~ mean(.x), .names = "{.col}")) %>%
                          ungroup() %>%
                          select(-country) %>% distinct(),
                  ., by = "cluster") %>%
        
        
        # get distinct values for each cluster
        pivot_longer(cols = -cluster, names_to = "var", values_to = "values") %>%
        mutate(cluster_name = str_c("Group ", cluster))

# inspect
chart_data


#///////////////


# get obj/sub_obj segment coordinates

# note that in the entire fmir dataset, the maximum sub_obj per obj is 3, so will have 3 segment_tbls,
# and for those obj with less than 3 sub_obj, the segment_tbl will just have NA coordinates and no output in ggplot
fmir %>% distinct(obj_num, sub_obj_num) %>% 
        add_count(obj_num, name = "sub_obj_count") %>% 
        distinct(obj_num, sub_obj_count) %>%
        arrange(desc(sub_obj_count)) 

# segment_1_tbl
# segment_1_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(obj_num, sub_obj_num) %>%
#         mutate(nrow = nrow(.),
#                x = nrow + 1.5,
#                xend = nrow + 1.5,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5) %>%
#         select(x, xend, y, yend) %>%
#         slice(1)
# segment_1_tbl

# # segment_2_tbl
# segment_2_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                 cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(1)
# segment_2_tbl
#                
# # segment_3_tbl
# segment_3_tbl <- fmir %>% filter(obj_num == current_obj_num) %>%
#         distinct(sub_obj_num, indicator_name) %>% count(sub_obj_num, name = "indicator_count") %>%
#         mutate(sub_obj_count = nrow(.),
#                cum_indicator_count = cumsum(indicator_count),
#                x = 1.5 + sub_obj_count + cum_indicator_count,
#                xend = 1.5 + sub_obj_count + cum_indicator_count,
#                y = .5,
#                yend = (chart_data %>% count(cluster) %>% nrow()) + .5,
#                sub_obj_index = row_number()) %>%
#         left_join(tibble(segment = c(1, 2, 3)),
#                   ., by = c("segment" = "sub_obj_index")) %>%
#         slice(2)
# segment_3_tbl


#///////////////////////


# plot

# get order_obj_and_sub_obj_avg_first()
order_obj_and_sub_obj_avg_first <- function(var){
        
        tibble(var = var) %>% 
                left_join(., fmir %>% distinct(indicator_name, concept) %>% 
                                  add_group_index(group_vars = concept, group_name = "concept_index"), 
                          by = c("var" = "indicator_name")) %>%
                mutate(order = case_when(str_detect(string = var, pattern = "^obj_") ~ 1,
                                         str_detect(string = var, pattern = "^sub_obj_[0-9]_[0-9]$") ~ 2,
                                         TRUE ~ 3)) %>%
                arrange(order, concept_index, var) %>% pull(var)
}

# create chart
fmir_obj_c_indicator_cluster_tile_chart <- chart_data %>% 
        ggplot(data = ., mapping = aes(x = fct_relevel(.f = var, order_obj_and_sub_obj_avg_first), 
                                       y = factor(cluster_name), fill = values, label = round_to_digits(values, digits = 2))) + 
        geom_tile(color = "#DDDDDD") + 
        geom_text(color = "#ffffff", size = 1.75, family = "Calibri") +
        scale_y_discrete(expand = c(0, 0)) +
        # geom_segment(x = segment_1_tbl$x, xend = segment_1_tbl$xend, 
        #              y = segment_1_tbl$y, yend = segment_1_tbl$yend, color = "#ffffff", size = 3) +
        # geom_segment(x = segment_2_tbl$x, xend = segment_2_tbl$xend,
        #              y = segment_2_tbl$y, yend = segment_2_tbl$yend, color = "#CBCBCB", size = 3) +
        # geom_segment(x = segment_3_tbl$x, xend = segment_3_tbl$xend,
        #              y = segment_3_tbl$y, yend = segment_3_tbl$yend, color = "#CBCBCB", size = 3) +
        labs(x = "Indicator", y = NULL, fill = "Avg. score\n(higher is better)\n") +
        coord_fixed(ratio = 2 / 5, clip = "off") +
        theme_bw() +
        # theme_minimal() +
        theme(
                # plot.background = element_rect(fill = "blue"),
                # plot.background = element_rect(colour = "#DDDDDD", size = .5),
                plot.margin = unit(c(0, 0, 0, 5), "mm"),
                plot.caption = element_text(hjust = 0, size = 11, face = "plain", family = "Calibri", 
                                            color = "#595959", margin = margin(t = 4, r = 0, b = 0, l = 0)),
                # text = element_text(family = "Calibri", size = 46, face = "plain", color = "#000000"),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                # panel.grid.major.y = element_line(color = "#DDDDDD"),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                # panel.border = element_rect(color = "#DDDDDD", size = .5),
                # panel.grid = element_blank(),
                # line = element_blank(),
                # rect = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_line(color = "#333333", size = .25),
                # axis.ticks.length.y.left = unit(.2, "cm"),
                axis.ticks.length.x.bottom = unit(.1, "cm"),
                axis.text.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                           margin = margin(t = 3, r = 0, b = 0, l = 0), angle = 45, hjust = 1),
                axis.text.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
                # axis.line.x.bottom = element_line(color = "#333333"),
                axis.line.x.bottom = element_blank(),
                axis.line.y.left = element_blank(),
                axis.title.x = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333", 
                                            margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(family = "Calibri", face = "plain", size = 7, color = "#333333",
                                            margin = margin(t = 0, r = 10, b = 0, l = 0)),
                # axis.title.y = element_blank(),
                plot.title = element_text(size = 7, face = "plain", hjust = .5, family = "Calibri", color = "#333333", 
                                          margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "pt")),
                legend.position = "right",
                legend.key.size = unit(4, "mm"),
                legend.title = element_text(size = 7, family = "Calibri", face = "plain", color = "#333333", hjust = .5,
                                            margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
                legend.text = element_text(size = 7, family = "Calibri", margin(t = 0, r = 0, b = 5, l = 0, unit = "pt"), 
                                           hjust = .5, color = "#333333")
                # legend.spacing.y = unit(5.5, "cm"),
                # legend.key = element_rect(size = 5),
                # legend.key.size = unit(2, 'lines')
        ) 
# fmir_obj_c_indicator_cluster_tile_chart <- ggdraw(align_legend(fmir_obj_c_indicator_cluster_tile_chart))

# inspect
fmir_obj_c_indicator_cluster_tile_chart


#//////////////////////////////////


# save chart as emf
filename <- tempfile(fileext = ".emf")
emf(file = filename)
print(fmir_obj_c_indicator_cluster_tile_chart)
dev.off()

# add emf to word doc
read_docx() %>%
        body_add_img(src = filename, width = 6, height = 6) %>%
        print(target = "output/charts/fmir_obj_c_indicator_cluster_tile_chart.docx")


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////


# fmir_sub_obj_1_1_trajectory_cluster ####

# https://cran.r-project.org/web/packages/traj/vignettes/trajVignette.pdf

# set current_sub_obj_num
current_sub_obj_num <- "sub_obj_1_1"
current_sub_obj_num

# get fmir_sub_obj_1_1_trajectory
fmir_sub_obj_1_1_trajectory <- fmir %>% filter(sub_obj_num == current_sub_obj_num,
                                               mcp_grouping != "EU-15") %>% 
        distinct(country, year, sub_obj_num, sub_obj_avg) %>%
        arrange(country, year) %>%
        unite(sub_obj_num_x_year, sub_obj_num, year, sep = "_x_") %>%
        pivot_wider(id_cols = country, names_from = sub_obj_num_x_year, values_from = sub_obj_avg) %>%
        arrange(country) %>%
        mutate(ID = row_number()) %>% relocate(ID, .before = everything()) %>% select(-country)


#////////////////////////


# inspect
fmir_sub_obj_1_1_trajectory
fmir_sub_obj_1_1_trajectory %>% glimpse()
fmir_sub_obj_1_1_trajectory %>% nrow() # 28
fmir_sub_obj_1_1_trajectory %>% ncol() # 12


#//////////////////////////////////////////////////////////////////////////////////////////////////


# get fmir_sub_obj_1_1_time
fmir_sub_obj_1_1_time <- map(.x = 1:((fmir_sub_obj_1_1_trajectory %>% ncol) - 1), 
                             .f = ~ tibble(!!sym(str_c("time_", .x)) := rep(.x, times = (fmir_sub_obj_1_1_trajectory %>% nrow())))) %>%
        bind_cols() %>%
        mutate(ID = row_number()) %>% relocate(ID, .before = everything())


#//////////////////


# inspect
fmir_sub_obj_1_1_time
fmir_sub_obj_1_1_time %>% glimpse()
fmir_sub_obj_1_1_time %>% nrow()
fmir_sub_obj_1_1_time %>% ncol()


#//////////////////////////////////////////////////////////////////////////////////////////////////


# step 1: get 24 summary measures for the trajectories
step_1_measures <- step1measures(Data = fmir_sub_obj_1_1_trajectory, Time = fmir_sub_obj_1_1_time, ID = TRUE)
step_1_measures
step_1_measures %>% attributes()

# step 2: use factor analysis to get subset of 24 measures that best describes the main features of the trajectories
step_2_factors <- step2factors(trajMeasures = step_1_measures)
step_2_factors
step_2_factors %>% attributes()
step_2_factors$factors %>% head()

# step 3: cluster based on the selected subset of trajectory measures
# step_3_clusters <- step3clusters(step_2_factors)
step_3_clusters <- step3clusters(step_2_factors, nclusters = 3)
step_3_clusters
step_3_clusters %>% attributes()


#///////////////////////////


# plot trajectory clusters
?plot.traj
plot(step_3_clusters, num.samples = 10)
plot(step_3_clusters, num.samples = 2, color.vect = c("#ff0000", "#00ff00"))
plot(step_3_clusters, num.samples = 10, clust.num = 2)

# plot trajectory mean
?plotMeanTraj
?plotCom
plotMeanTraj(step_3_clusters, clust.num = NULL)
plotCombTraj(step_3_clusters)

# plot trajectory distribution
?plotMedTraj
plotMedTraj(step_3_clusters)
plotBoxplotTraj(step_3_clusters)


#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////




