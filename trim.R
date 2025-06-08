source("classes.R")
source("sim_study_funcs.R")

bank_data <- read.csv("data/banking.csv", sep = ";")

percentiles <- c(0.01, 0.05, 0.10, 0.20)
dist_names = c('Unbiased Dependant', 'Euclidean', 'Gower', 'Cosine', 'HVDM', 'Ahmad & Dey')
meth_names = c('pam', 'kmedoids')
score_names = c('ARI', 'AMI', 'Sil', 'Sil SD')
nr = length(dist_names)
nc = length(meth_names)
ns = length(score_names)

total_evals = array(c(1:144)*0,dim=c(nr,nc,ns, 4), dimnames = list(dist_names, meth_names, score_names, percentiles))

for (n in 1:200){
  print(n)
  rand_bank <- sample(nrow(bank_data), nrow(bank_data)%/%100)
  bank_samp <- bank_data[rand_bank, ]
  bank_df <- bank_samp[, -17]
  bank_clus <- bank_samp[, 17]

  bank_evals <- array(dim=c(nr,nc,ns, 4), dimnames = list(dist_names, meth_names, score_names, percentiles))
  bank_Dist <- Distance$new(bank_samp, c(1,6,12:15), c(5,7:8), c(2:4,9:11,16), target_col=17)
  num_clust <- 2

  for (i in 1:nr){
    dist = dist_names[i]
    dist_mat = as.matrix(bank_Dist$dist_matrix(dist))
    avg_dist <- rowMeans(dist_mat)
    
    for (k in 1:4){
      x <- percentiles[k]
      # find bottom and top 10% values
      lower <- quantile(avg_dist, x)
      upper <- quantile(avg_dist, 1-x)
      
      # trim bank_df based on average distance
      trimmed_dist_mat <- dist_mat[avg_dist >= lower & avg_dist <= upper, avg_dist >= lower & avg_dist <= upper]
      trimmed_bank_clus <- bank_clus[avg_dist >= lower & avg_dist <= upper]
      
      for (j in 1:nc){
        meth <- meth_names[j]
        
        bank_Clustering <- Clustering$new(trimmed_dist_mat, num_clust)
        bank_labels <- bank_Clustering$cluster(meth)
        
        bank_Scores <- Scores$new(true_clus=trimmed_bank_clus, labels=bank_labels, dist_mat=trimmed_dist_mat)
        
        bank_evals[i, j, 1, k] <- bank_Scores$ari()
        bank_evals[i, j, 2, k] <- bank_Scores$ami()
        bank_evals[i, j, 3, k] <- mean(bank_Scores$sil())
        bank_evals[i, j, 4, k] <- sd(bank_Scores$sil())

      }
    }
  }
  total_evals <- total_evals + bank_evals
}

total_evals <- total_evals / 200

saveRDS(total_evals, "test_evals_trim2.rds")