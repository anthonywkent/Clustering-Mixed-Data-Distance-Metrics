library(cluster)
library(MASS)
library(arules)
library(aricode)
library(knitr)
library(kmed)
library(lsa)
library(nomclust)
library(UBL)
library(stats)
library(datasets)
library(manydist)
library(R6)
library(FactoMineR)
library(MixSim)

# Categorized numerical variables: function
intv <- function(vec, class) {
  nbase <- (1:(class-1))/class
  nq <- numeric(length(nbase))
  for (i in 1:length(nq)) {
    nq[i] <- quantile(vec, nbase[i])
  }
  res <- c(min(vec), nq, max(vec)) 
  res[1] <- res[1]-1
  for (i in 2:length(res)){
    if (res[i-1]==res[i]){
      res[i] <- res[i]+2e-15
    }
  }
  return(res)
}

# Simulate dataset with 'true' clustering: function
sim_data <- function(clusters=2, overlap=0.02, rows=400, columns=14, pi_val=0.3, cat_ratio=0.4){
  # Construct artificial data set - sphericity set to FALSE by default
  mixsimaux <- MixSim(BarOmega = overlap, PiLow=pi_val,
                      K = clusters, p = columns, resN = 1000000)
  mixdtaux <- simdataset(n = rows, Pi = mixsimaux$Pi, Mu = mixsimaux$Mu, S = mixsimaux$S)
  # Discretise first half attributes using 4 levels
  for (k in 1:(round(columns*cat_ratio))){
    mixdtaux$X[,k] <- 
      as.factor(cut(mixdtaux$X[,k], intv(mixdtaux$X[,k], 4), labels = (1:4)))
  }
  mixdt1df <- as.data.frame(mixdtaux$X)
  id <- mixdtaux$id
  for (k in 1:(round(columns*cat_ratio))){
    mixdt1df[,k] <- as.factor(mixdt1df[,k])
  }
  # this temporarily fixes a nasty bug of cluspcamix
  colnames(mixdt1df) <- sprintf("a%d", 1:columns)
  
  # Return data and ids
  mixdtaux$X <- mixdt1df
  mixdtaux}


# Calculate evaluation metrics: function
run_evals <- function(args, nreps=1){
  clusters = as.integer(args[1])
  overlap = as.double(args[2])
  rows = as.integer(args[3])
  columns = as.integer(args[4])
  pi_val = as.double(args[5])
  cat_ratio = as.double(args[6])
  
  
  dist_names = c('Unbiased Dependant', 'Euclidean', 'Gower', 'Cosine', 'HVDM', 'Ahmad & Dey')
  meth_names = c('pam', 'kmedoids')
  score_names = c('ARI', 'AMI', 'Silhouette')
  
  nr = length(dist_names)
  nc = length(meth_names)
  ns = length(score_names)
  
  evals = array(dim=c(nr,nc,ns, nreps), dimnames = list(dist_names, meth_names, score_names, c(1:nreps)))
  for (n in 1:nreps){
    # set seed
    set.seed(1234+n)
    # generate data
    mixdtaux <- sim_data(clusters, overlap, rows, columns, pi_val, cat_ratio)
    # Distance metric
    Dist = Distance$new(mixdtaux$X, num_cols=c(as.integer(round(cat_ratio*columns)+1):columns), cat_cols=c(1:round(cat_ratio*columns)))
    
      for (i in 1:nr){
        try(
        for (dummy in 1:1){
          dist = dist_names[i]
          dist_mat = Dist$dist_matrix(dist)
          
          for (j in 1:nc){
            meth = meth_names[j]
            Clust = Clustering$new(dist_mat, clusters)
            labels = Clust$cluster(meth)
            
            Score = Scores$new(true_clus=mixdtaux$id, labels=labels, dist_mat=dist_mat)
            
            evals[i,j,1,n] = Score$ari()
            evals[i,j,2,n] = Score$ami()
            evals[i,j,3,n] = mean(Score$sil())
          }
        }      
      )
    }
  }
  evals
}

# Run the study: function
study <- function(char, values, nreps){
  dist_names = c('Unbiased Dependant', 'Euclidean', 'Gower', 'Cosine', 'HVDM', 'Ahmad & Dey')
  meth_names = c('pam', 'kmedoids')
  score_names = c('ARI', 'AMI', 'Silhouette')
  
  
  nv = length(values)
  store = array(dim=c(6,2,3,2,nv), dimnames = list(dist_names, meth_names, score_names, c('mean', 'sd'), c(1:nv)))
  
  args = list(clusters=2, overlap=0.02, rows=400, columns=14, pi_val=0.3, cat_ratio=0.4)
  
  for (i in 1:nv){
    args[[char]] = values[i]
    evals = run_evals(args, nreps=nreps)
    evals_mean = apply(evals, c(1,2,3), FUN=mean)
    evals_sd = apply(evals, c(1,2,3), FUN=sd)
    
    store[,,,1,i] = evals_mean
    store[,,,2,i] = evals_sd
  }
  store
}
