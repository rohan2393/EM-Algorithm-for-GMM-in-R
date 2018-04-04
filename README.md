# EM-Algorithm-for-GMM-in-R
Manually coded in R using initial values from K-means algorithm

## Output chart for convergence:
We ran the EM algorithm for convergence threshold equal to 0.000001, 0.0001 and 0.01 for k = 2 to 10. And we came to a conclusion that with a 4 digit convergence threshold i.e. (0.0001) there are 9 clusters(k=9) with maximum iterations 20.
The converged log-likelihood value for 9 clusters is 140623.3978

## Log-likelihood for each iteration are as follows

![Alt text](https://github.com/rohan2393/EM-Algorithm-for-GMM-in-R/blob/master/Log-likelihood_Chart.png)
 
## Summary of the Project
The main objective of this project was to build an Expectation Maximization Algorithm manually using R.  
The process to build the Expectation Maximization algorithm involved the following steps:
•	Running K-means Algorithm
•	Obtaining Initial Parameter values from the K-means Algorithm for the Initialization step in EM algorithm
•	Executing the EM Algorithm

## About the Dataset and Cleaning the Dataset
We are using cluster.csv dataset for our project. The size of the dataset is 1077 observations and 71 dimensions.  Before executing the K-means algorithms we first cleaned the dataset. The dataset contained some NA values which were removed from the dataset by completely eliminating the row. The resulting clean dataset’s size is 1047 observations and 71 dimensions. Later on dimension ID was removed which resulted in 70 dimensions.

## Running K-means Algorithm
After cleaning the dataset, we first initialized the number of clusters we want to start with in the K-means Algorithm. To begin with we started with k=2(in our case we defined it as n_clusters = 2). Then we ran the K-means with n_clusters = 2. We obtained the Cluster mean values for each cluster in each dimension. This is the basis for the EM algorithm.

## Obtaining Initial Parameter values from the K-means Algorithm for the Initialization step in EM algorithm
To begin with the EM algorithm we need some initial values for the parameters which were obtained from the K-means Algorithm. 
•	 μk  (in our case we defined it as mu_initials)
•	 Covariance Σk  (in our case we defined it as sigma_initials)
•	  πk  (in our case we defined it as Pi_k_initials)
•	 Nk  (which is nothing but the total data points for k cluster)

Note: k stands for n_clusters

## Executing the EM Algorithm
The Expectation-Maximization algorithm for Gaussian mixture is an iterative algorithm that starts from some initial parameter values and then proceeds to iteratively update either until convergence is detected or out of max iterations. Each iteration consists of an E-step and an M-step.
In our case convergence would be reached if log_likelihood is lower than the threshold value that is, in our case, 0.0001. Keeping the time constraint in mind we are iterating the model 20 times. 

## EM Algorithm has 4 steps:
1.	Determining the Initial GMM Parameters (μk, Covariance Σk  , πk , Nk )
2.	Expectation Step
3.	Maximization Step
4.	Check for Convergence
