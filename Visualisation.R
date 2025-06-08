# Visualisation of Real-World Data using MDS, FAMD, and T-SNE

# Libraries Required
library(cluster)
library(MASS)
library(kmed)
library(lsa)
library(nomclust)
library(UBL)
library(stats)
library(datasets)
library(manydist)
library(ggplot2)
library(Rtsne)
library(FactoMineR)
library(tidyverse)
set.seed(2025)

# functions and classes
source('classes.R')
source('sim_study_funcs.R')

# Import data and clean
cust_data <- read.csv('../data/customer_segmentation.csv', sep=",")
cust_data <- na.omit(cust_data)
cust_data <- cust_data[-c(which(cust_data$Ever_Married == "")),]
cust_data <- cust_data[-c(which(cust_data$Graduated == "")),]
cust_data <- cust_data[-c(which(cust_data$Profession == "")),]
cust_data <- cust_data[-c(which(cust_data$Var_1 == "")),]

# Take random subset of data (computationally quicker)
rand_cust = sample(nrow(cust_data), nrow(cust_data)%/%10)
cust_df <- cust_data[,-c(1,11)][rand_cust,]
cust_clus <- cust_data[,11][rand_cust]

# take only distinct values
cust_clus <- cust_clus[!duplicated(cust_df)]
cust_df <- cust_df[!duplicated(cust_df),]

for (x in c(1:2,4:5,7,9)) {
  cust_df[,x] <- factor(cust_df[,x])
}

# Gower + MDS
gower_dist_cust <- daisy(cust_df, metric="gower", 
                         type=list(asymm=c(1:2,4), factor=c(5,7,9)))
gower_mds_cust <- isoMDS(gower_dist_cust)$points

# Gap Statistic to find number of clusters
gap_cust = clusGap(cust_num_df, pam, K.max=12, B=50)
clusters_globalSEmax <- maxSE(gap_cust$Tab[, "gap"], gap_cust$Tab[, "SE.sim"], method = "globalSEmax")
clusters_globalmax <- maxSE(gap_cust$Tab[, "gap"], gap_cust$Tab[, "SE.sim"], method = "globalmax")

labels_globalSEmax = pam(gower_dist_cust, clusters_globalSEmax, diss = TRUE)$clustering
labels_globalmax = pam(gower_dist_cust, clusters_globalmax, diss = TRUE)$clustering


# Convert to data frame suitable for ggplot2
mds_SEmax <- data.frame(X1 = gower_mds_cust[,1], 
                     X2 = gower_mds_cust[,2], 
                     Label = as.factor(labels_globalSEmax))
names(mds_SEmax) <- c("Dim1", "Dim2", "Label")

# Plot and save
plot_SEmax <- ggplot(mds_SEmax, aes(x = Dim1, y = Dim2, color = Label)) + 
          geom_point() +
          labs(title='Gower MDS: Customer Segmentation Data', x = "Dim 1", y = "Dim 2") +
          scale_color_brewer(palette = "Paired")
ggsave('cust_mds_SEmax.png', plot_SEmax)

# Repeat with other labels
mds_max <- data.frame(X1 = gower_mds_cust[,1], 
                     X2 = gower_mds_cust[,2], 
                     Label = as.factor(labels_globalmax))
names(mds_max) <- c("Dim1", "Dim2", "Label")

# Plot and save
plot_max <- ggplot(mds_max, aes(x = Dim1, y = Dim2, color = Label)) + 
  geom_point() +
  labs(title='Gower MDS: Customer Segmentation Data', x = "Dim 1", y = "Dim 2") +
  scale_color_brewer(palette = "Paired")
ggsave('cust_mds_max.png', plot_max)



# FAMD
famd_res <- FAMD(cust_df, graph = FALSE)

# First set of labels
famd_SEmax <- data.frame(X1 = famd_res$ind$coord[,1], 
                      X2 = famd_res$ind$coord[,2], 
                      Label = as.factor(labels_globalSEmax))
names(famd_SEmax) <- c("Dim1", "Dim2", "Label")

plot_FAMD_SEmax <- ggplot(famd_SEmax, aes(x = Dim1, y = Dim2, color = Label)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "FAMD Projection: Customer Segmentation Data", x = "Dim 1", y = "Dim 2")
ggsave('cust_FAMD_SEmax.png', plot_FAMD_SEmax)

# Second set of labels
famd_max <- data.frame(X1 = famd_res$ind$coord[,1], 
                      X2 = famd_res$ind$coord[,2], 
                      Label = as.factor(labels_globalmax))
names(famd_max) <- c("Dim1", "Dim2", "Label")

plot_FAMD_max <- ggplot(famd_max, aes(x = Dim1, y = Dim2, color = Label)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "FAMD Projection: Customer Segmentation Data", x = "Dim 1", y = "Dim 2")
ggsave('cust_FAMD_max.png', plot_FAMD_max)



#T-SNE 
tsne_res <- Rtsne(gower_dist_cust, is_distance = TRUE)

# First set of labels
tsne_SEmax <- data.frame(X1 = tsne_res$Y[,1], 
                         X2 = tsne_res$Y[,2], 
                         Label = as.factor(labels_globalSEmax))
names(tsne_SEmax) <- c("Dim1", "Dim2", "Label")

plot_TSNE_SEmax <- ggplot(aes(x = Dim1, y = Dim2), data = tsne_SEmax) +
  geom_point(aes(color = Label)) +
  labs(title = "T-SNE Projection: Customer Segmentation Data", x = "Dim 1", y = "Dim 2")
ggsave('cust_TSNE_SEmax.png', plot_TSNE_SEmax)


tsne_max <- data.frame(X1 = tsne_res$Y[,1], 
                       X2 = tsne_res$Y[,2], 
                       Label = as.factor(labels_globalmax))
names(tsne_max) <- c("Dim1", "Dim2", "Label")

plot_TSNE_max <- ggplot(aes(x = Dim1, y = Dim2), data = tsne_max) +
  geom_point(aes(color = Label)) +
  labs(title = "T-SNE Projection: Customer Segmentation Data", x = "Dim 1", y = "Dim 2")
ggsave('cust_TSNE_max.png', plot_TSNE_max)


