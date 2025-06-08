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

Distance <- R6Class(
  "Distance",
  public = list(
    data = NULL,
    num_cols = NULL,
    bin_cols = NULL,
    cat_cols = NULL,
    target_col = NULL,
    clus = NULL,
    
    # Constructor
    initialize = function(data, num_cols=NULL, bin_cols=NULL, cat_cols=NULL, target_col=NULL, clus=NULL) {
      self$data <- data
      self$num_cols <- num_cols
      self$bin_cols <- bin_cols
      self$cat_cols <- cat_cols
      self$target_col <- target_col
      self$clus <- clus
      self
    },
    
    sep_target = function() {
      if (is.null(self$target_col)){}
      else {
        data <- self$data
        self$clus <- data[,self$target_col]
        self$data <- data[,-self$target_col]
        self$target_col <- NULL
      }
    },
    
    
    # calculate distance
    dist_matrix = function(dist) {
      self$sep_target()
      data <- self$data
      if (dist=='Euclidean'){
        for (x in c(1:ncol(data))){data[,x] <- as.numeric(data[,x])}
        
        dist(data, method="euclidean")
      }
      else if (dist=='Gower'){
        for (x in c(self$bin_cols,self$cat_cols)) {data[,x] <- factor(data[,x])}
        
        if (is.null(self$bin_cols)){daisy(data, metric='gower')}
        else{daisy(data, metric='gower', type=list(asymm=self$bin_cols))}
      }
      else if (dist=='Cosine'){
        for (x in c(self$bin_cols,self$cat_cols)) {data[,x] <- factor(data[,x])}
        for (x in c(1:ncol(data))){data[,x] <- as.numeric(data[,x])}
        
        1 - cosine(t(data))
      }
      else if (dist=='HVDM') {
        if (!(is.null(self$bin_cols))){for (x in c(self$bin_cols)) {data[,x] <- factor(data[,x])}}
        if (!(is.null(self$cat_cols))){for (x in c(self$cat_cols)) {data[,x] <- factor(data[,x])}}
        distances(1, data, 'HVDM')
      }
      else if (dist=='Ahmad & Dey'){
        distmix(data, method='ahmad', idnum=self$num_cols, idbin=self$bin_cols, idcat=self$cat_cols)
      }
      else if (dist=='Unbiased Dependant'){
        mdist(data, preset = 'unbiased_dependent')
      }
      else {
        print('Method not supported! Current acceptable methods are: "Euclidean", "Gower","Cosine", "HVDM", "Ahmad & Dey", "Unbiased Dependant".')
      }
    }
  )
)


Clustering <- R6Class(
  "Clustering",
  public = list(
    dist_mat = NULL,
    clusters = NULL,
    func_res = NULL,
    
    
    # Constructor
    initialize = function(dist_mat, clusters=NULL) {
      self$dist_mat = dist_mat
      self$clusters = clusters
    },
    
    cluster = function(meth='pam') {
      if (meth=='pam'){
        self$func_res <- pam(self$dist_mat, self$clusters, diss=TRUE)
        (self$func_res)$clustering
        
      }
      else if (meth=='kmedoids'){
        self$func_res <- skm(self$dist_mat, self$clusters)
        (self$func_res)$cluster
      }
      else{
        print('Method not supported! Current acceptable methods are: "pam", "kmedoids".')
      }
    }
  )
)


Scores <- R6Class(
  "Scores",
  public = list(
    dist_mat = NULL,
    labels = NULL,
    true_clus = NULL,
    
    
    # Constructor
    initialize = function(true_clus=NULL, labels=NULL, dist_mat=NULL) {
      self$true_clus = true_clus
      self$labels = labels
      self$dist_mat = dist_mat
      
    },
    
    ari = function() {
      if (is.null(self$true_clus) | is.null(self$labels)){
        print('Missing information! Need true clustering and numerical clustering labels!')
      }
      else{ARI(self$true_clus, self$labels)}
    },
    ami = function() {
      if (is.null(self$true_clus) | is.null(self$labels)){
        print('Missing information! Need true clustering and numerical clustering labels!')
      }
      else{AMI(self$true_clus, self$labels)}
    },
    sil = function(avg=FALSE){
      if (is.null(self$true_clus) | is.null(self$labels)){
        print('Missing information! Need true clustering and numerical clustering labels!')
      }
      else{
        sil_score <- silhouette(self$labels, self$dist_mat)
        if (avg==TRUE){
          mean(sil_score[,3])
        }
        else{
          sil_score[,3]
        }
      }
    }
  )
)