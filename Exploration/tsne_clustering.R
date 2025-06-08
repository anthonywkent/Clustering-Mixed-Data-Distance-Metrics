#import and preprocess data
cust_data <- read.csv('../data/customer_segmentation.csv', sep=",")
cust_data <- na.omit(cust_data)

cust_data <- cust_data[-c(which(cust_data$Ever_Married == "")),]
cust_data <- cust_data[-c(which(cust_data$Graduated == "")),]
cust_data <- cust_data[-c(which(cust_data$Profession == "")),]
cust_data <- cust_data[-c(which(cust_data$Var_1 == "")),]

rand_cust = sample(nrow(cust_data), nrow(cust_data))
cust_df <- cust_data[,-c(1,11)][rand_cust,]
cust_clus <- cust_data[,11][rand_cust]

cust_clus <- cust_clus[!duplicated(cust_df)]
cust_df <- cust_df[!duplicated(cust_df),]

for (x in c(1:2,4:5,7,9)) {
  cust_df[,x] <- factor(cust_df[,x])
}


#compute tsne projections
tsne_res1 <- Rtsne(cust_df, dims=1)
tsne_res2 <- Rtsne(cust_df, 2)
tsne_res3 <- Rtsne(cust_df, dims=3)

#extract embeddings
data_tsne1 = tsne_res1$Y
data_tsne2 = tsne_res2$Y
data_tsne3 = tsne_res3$Y


# cluster data using kmeans
km_tsne1_labels_globalSEmax = kmeans(data_tsne1, clusters_globalSEmax)$cluster
km_tsne1_labels_globalmax = kmeans(data_tsne1, clusters_globalmax)$cluster

km_tsne2_labels_globalSEmax = kmeans(data_tsne2, clusters_globalSEmax)$cluster
km_tsne2_labels_globalmax = kmeans(data_tsne2, clusters_globalmax)$cluster

km_tsne3_labels_globalSEmax = kmeans(data_tsne3, clusters_globalSEmax)$cluster
km_tsne3_labels_globalmax = kmeans(data_tsne3, clusters_globalmax)$cluster

#evaluate stats
dim2_names = c('Dimensions', 'ARI', 'AMI', 'Sil')
tsne_evals1 = array(dim=c(3,4), dimnames = list(c(1:3), dim2_names))
tsne_evals2 = array(dim=c(3,4), dimnames = list(c(1:3), dim2_names))

tsne_evals1[1,1] = 1
tsne_evals1[1,2] = ARI(km_tsne1_labels_globalSEmax, cust_clus)
tsne_evals1[1,3] = AMI(km_tsne1_labels_globalSEmax, cust_clus)
tsne_evals1[1,4] = mean(silhouette(km_tsne1_labels_globalSEmax, dist(data_tsne1, method="euclidean"))[,3])

tsne_evals1[2,1] = 2
tsne_evals1[2,2] = ARI(km_tsne2_labels_globalSEmax, cust_clus)
tsne_evals1[2,3] = AMI(km_tsne2_labels_globalSEmax, cust_clus)
tsne_evals1[2,4] = mean(silhouette(km_tsne2_labels_globalSEmax, dist(data_tsne2, method="euclidean"))[,3])

tsne_evals1[3,1] = 3
tsne_evals1[3,2] = ARI(km_tsne3_labels_globalSEmax, cust_clus)
tsne_evals1[3,3] = AMI(km_tsne3_labels_globalSEmax, cust_clus)
tsne_evals1[3,4] = mean(silhouette(km_tsne3_labels_globalSEmax, dist(data_tsne3, method="euclidean"))[,3])

tsne_evals2[1,1] = 1
tsne_evals2[1,2] = ARI(km_tsne1_labels_globalmax, cust_clus)
tsne_evals2[1,3] = AMI(km_tsne1_labels_globalmax, cust_clus)
tsne_evals2[1,4] = mean(silhouette(km_tsne1_labels_globalmax, dist(data_tsne1, method="euclidean"))[,3])

tsne_evals2[2,1] = 2
tsne_evals2[2,2] = ARI(km_tsne2_labels_globalmax, cust_clus)
tsne_evals2[2,3] = AMI(km_tsne2_labels_globalmax, cust_clus)
tsne_evals2[2,4] = mean(silhouette(km_tsne2_labels_globalmax, dist(data_tsne2, method="euclidean"))[,3])

tsne_evals2[3,1] = 3
tsne_evals2[3,2] = ARI(km_tsne3_labels_globalmax, cust_clus)
tsne_evals2[3,3] = AMI(km_tsne3_labels_globalmax, cust_clus)
tsne_evals2[3,4] = mean(silhouette(km_tsne3_labels_globalSEmax, dist(data_tsne3, method="euclidean"))[,3])



saveRDS(tsne_evals1, "tsne_evals_SEmax")
saveRDS(tsne_evals2, "tsne_evals_max")


