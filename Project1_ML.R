# loading required libraries
library(dplyr)
library(mvtnorm)

# Getting the dataset into R
proj_data <-  read.csv("/Users/rohan/Desktop/MLProject1/cluster.csv",sep = ",", header = F)

# Exploring the dataset
head(proj_data)
dim(proj_data)
summary(proj_data)


# Removing the first column from the dataset
proj_data <- proj_data[-1]


# Checking for missing/NA values
sum(is.na(proj_data))

# Removing na values from the dataset
proj_data <- na.omit(proj_data)
dim(proj_data)


# Setting the number of clusters for K-means
n_clusters <- 2

# Getting the dimensions of the dataset
n_dimensions <- ncol(proj_data)


# K means clustering to get the initial parameters (mu_initial, sigma_initial)
cluster_pred <- kmeans(x = proj_data, centers = n_clusters)


# Getting initial mu values
mu_initials <- cluster_pred[2]
mu_initials <- mu_initials$centers

# Adding the cluster prediction to the dataset
proj_data_df <- data.frame(cbind(proj_data, cluster = cluster_pred$cluster))


# Getting the initial covariance array
sigma_initials <- array(rep(0),dim = c(n_dimensions, n_dimensions, n_clusters))
for (i in 1:n_clusters){
  sigma_initials[,,i] <- cov(proj_data_df[proj_data_df$cluster == i, -71])
}


# Getting the initial Pi_k values
proj_data_df_by_cluster <- proj_data_df %>% group_by(cluster)
Pi_k_initials <- proj_data_df_by_cluster %>% summarise(size = n()) %>% mutate(Pi_k = size / sum(size))
Pi_k_initials <- Pi_k_initials$Pi_k


# empty data structures for comp & responsibility
comp <- matrix(rep(0), nrow = nrow(proj_data), ncol = n_clusters)
responsibility <- matrix(rep(0), nrow = nrow(proj_data), ncol = n_clusters)
#dim(comp)




# _______________ Expectation Step _______________ #
# Defining expectation function
expectation <- function(x, mu_initials, sigma_initials, Pi_k_initials){
  
  for (j in 1:n_clusters){
    comp[,j] <- dmvnorm(x, as.matrix(mu_initials[j,]), as.matrix(sigma_initials[,,j])) * Pi_k_initials[j]
  }
  
  sum_of_comps <- rowSums(comp)
  responsibility <- comp / sum_of_comps
  
  ln_sum_of_comps <- log(sum_of_comps, base = exp(1))
  sum_of_ln_sum_of_comps <- sum(ln_sum_of_comps)
  
  list("log_likelihood" = sum_of_ln_sum_of_comps,
       "responsibility_df" = responsibility)
}




# _______________ Maximisation Step  _______________ #
# Creating empty data structures to store the values
Mu_new <- matrix(rep(0), nrow = n_clusters, ncol = n_dimensions)
N_k_new <- rep(0, n_clusters)
Pi_k_new <- rep(0, n_clusters)
sigma_t <- matrix(0, nrow = n_dimensions, ncol = n_dimensions)
sigma_t1 <- matrix(0, nrow = n_dimensions, ncol = n_dimensions)
Sigma_new <- array(rep(0),dim = c(n_dimensions, n_dimensions, n_clusters))


# Defining maximisation function
maximisation <- function(x, responsibility_df){
  
  for (k in 1:n_clusters){
    N_k_new[k] <- sum(responsibility_df[,k])
    Pi_k_new[k] <- N_k_new[k] / nrow(x)
  }
  
  for (k in 1:n_clusters){
    Mu_t <- ((1/N_k_new[k]) * (responsibility_df[,k] * x))
    for (i in 1:n_dimensions){
      Mu_new [k, i] <- sum(Mu_t[,i])
    }
  }  
  
  for (k in 1:n_clusters){
    for (j in 1:nrow(x)){
      difference_transposed <- t(x[j,] - Mu_new[k,])
      difference <- x[j,] - Mu_new[k,]
      co_efficients <- responsibility_df[j, k] / N_k_new[k]
      sigma_t <- co_efficients * (difference_transposed %*% as.matrix(difference))
      sigma_t1 <- sigma_t1 + sigma_t
    }
    Sigma_new[,,k] <- sigma_t1
  }
  
  list("Mu_new" = Mu_new,
       "Sig_new" = Sigma_new,
       "Pi_k_New" = Pi_k_new)
}




# _______________ Convergence Step _______________ #
  for (i in 1:20){
    if (i == 1){
      # Initialization step
      e_step <- expectation(proj_data, mu_initials, sigma_initials, Pi_k_initials)
      m_step <- maximisation(proj_data, e_step[["responsibility_df"]])
      current_log_likelihood <- e_step[["log_likelihood"]]
      loglikelihood_vector <- e_step[["log_likelihood"]]
    } 
    else{
      # Repeatation for E and M steps until convergence
      e_step <- expectation(proj_data, m_step[["Mu_new"]], m_step[["Sig_new"]], m_step[["Pi_k_New"]])
      m_step <- maximisation(proj_data, e_step[["responsibility_df"]])
      loglikelihood_vector <- c(loglikelihood_vector, e_step[["log_likelihood"]])
      loglikelihood_diff <- abs((current_log_likelihood - e_step[["log_likelihood"]]))
      
      # Check for Convergence
      if(loglikelihood_diff < 0.0001){
        break
      } 
      else{
        current_log_likelihood <- e_step[["log_likelihood"]]
      }
    }
  }