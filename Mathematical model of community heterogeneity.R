# model29_competition&immigration&heterogenerate&dispersion circulation

rm(list=ls(all=TRUE))
require(tidyverse)
library(gridExtra)
library(e1071)
library(beepr)

library(ggplot2)
library(dplyr)
library(tidyr)
library(moments)
library(ggridges)
library(vegan)
library(MASS)


#################################################################
#### Take into account random immigration rate, growth rate, microbial interactions
#################################################################


num_iterations <- 100  # Number of iterations for the loop
otu_table_list <- list()  # Create an empty list to store randomly generated OTU tables
res_ra_list <- list()  # Create an empty list to store simulated communities from the LV model
mean_dist_vector <- numeric(num_iterations) 
cv_mean_vector <- numeric(num_iterations)
Heterogeneity_vector <- numeric(num_iterations)

# Define the dispersion parameter values to test
dispersion_values <- c(0.08, 0.18, 0.34, 0.60, 1.00, 1.70, 3.50, 8.00, 18.00)
# A similar approach can be used to test the effects of sigma_mig and sigma_grate on heterogeneity

# Loop over each dispersion value
for (dispersion in dispersion_values)
{

 for (a in 1:num_iterations) 
 {
 ###Generate a random community matrix with a given heterogeneity level  
 num_samples <- 100  # Number of samples
 num_otus <- 100  # Number of OTUs
 
 mean_size <- 100  # Mean size parameter for the negative binomial distribution
 target_mean_min <- 0.1  # Minimum target mean Bray-Curtis distance
 target_mean_max <- 0.9  # Maximum target mean Bray-Curtis distance
 
 otu_table <- matrix(rnegbin(n = num_samples * num_otus, mu = mean_size, theta = dispersion), nrow = num_samples, ncol = num_otus)
 otu_table_ra <- otu_table / rowSums(otu_table)
  
  diss_matrix <- vegdist(otu_table_ra, method = "bray")
  
  # Calculate the mean of the distance matrix and check if it falls within the target range
  mean_dist <- mean(diss_matrix)
  if (mean_dist >= target_mean_min && mean_dist <= target_mean_max)
    {
    # The generated matrix meets the target mean range; save the relative abundance matrix
    otu_table_list[[a]] <- otu_table_ra
    
    std <-apply(otu_table_ra, 2, sd)
    mean_vals <- colMeans(otu_table_ra)
    cv <- std / mean_vals #Coefficient of variation for each OTU across particles; cv is a vector of length num_otus
    cv_mean <- mean(cv)
    cv_mean_vector[a] <- cv_mean
    mean_dist_vector[a]<- mean_dist
    }
  
 
 ### Simulate microbial community with the generalized Lotka-Volterra model include dispersal from a species pool
 np <- 100  # Number of particles
 nm <- 100  # Number of species
 
 mean_mig = 1  # Mean immigration coefficient for species i
 mean_grate = 1 # Mean growth rate for species i
 
 # After confirming that sigma_mig, sigma_grate, and sigma_K have little effect，we set them as constants for plotting
 sigma_mig = 0.1   # Standard deviation of immigration rate λ among species
 sigma_grate = 0.1   # Standard deviation of growth rate among species
 
 # For the standard deviation of immigration rate across particles, we use the cv vector from the generated random OTU matrix
 cv_partmrate = cv  # cv of immigration rate for each species across particles
 # Differences in growth rate across particles for the same species are negligible, for simplicity, set to 0.1
 cv_partgrate = 0.1  # Standard deviation of growth rate for the same species across particles, When particles are degradable, the sequential arrival of decomposers causes the coefficient of variation of this parameter to be very large, making it an important factor
 
 # produce variation here
 # The columns lambda and com are generated from normal distributions
 # The resulting parameter data frame 'pars' contains columns mag, lambda, and rm
 
 pars <- list( mag = 1:nm, lambda = mean_mig*rlnorm( nm,0, sigma_mig )*exp(-sigma_mig^2 / 2),
              rm =  mean_grate*rlnorm( nm,0, sigma_grate )*exp(-sigma_grate^2 / 2)) %>% 
   as.data.frame()
 
 tf <- 100  #  the duration of the experiment t

 # Generate an np × nm immigration rate matrix; rows are particles, columns are species
 lambda <- matrix(0, nrow = np, ncol = nm)
 for (i in 1:nm)
     {
  # Different particle environments correspond to different immigration rates

    lambda[ ,i] <- as.numeric(pars[i, 2]) * rlnorm(np, 0, cv_partmrate) / tf  
     }
 
 
 # Generate an np × nm growth rate matrix; rows are particles, columns are species
 r <- matrix(0, nrow = np, ncol = nm)
 for (i in 1:nm)
      {
   # Different particles have different growth rates
   
     r[ ,i] <- as.numeric(pars[i, 3]) * rlnorm(np, 0, cv_partgrate) / tf 
      }

 res <- matrix(0, nrow = np, ncol = nm) 

 # The outermost loop iterates over particles, np times
 for (p in 1:np ) 
       {

   ## This block of code implements the analytical succession of the generalized Lotka-Volterra model over time T for a single particle
    
     T <- 1000 # Number of iterations
  
     step <- tf / T # Iteration step size
    
     A_mean_range <- c(0.1, 0.1)
     A_mean <- seq(A_mean_range[1], A_mean_range[2], length.out = nm)
     AA <- matrix(runif(nm * nm) * A_mean * 2, nrow = nm, ncol = nm)
     AA <- AA - diag(diag(AA)) + diag(nm) 
     
     N_mean <- 1 / (nm * A_mean * sqrt(pi/2))
     sigma3 <- N_mean / 12
     NN <- matrix(0, nrow = T, ncol = nm)  
    
     # Initial species abundances
     N0 <- abs(rnorm(nm, N_mean, sigma3))
     NN[1, ] <- N0
    
     result <- matrix(0, nrow = T, ncol = nm)
     for (t in 2:T) {
       for (j in 1:nm) {
         k1 <- r[p,j] * NN[(t - 1), j] * (1 - sum(AA[j, ] * NN[(t - 1), ])) * step + lambda[p,j]
         k2 <- r[p,j] * (NN[(t - 1), j] + k1 / 2 ) * (1 - sum(AA[j, ] * (NN[(t - 1), ]))) * step + lambda[p,j]
         k3 <- r[p,j] * (NN[(t - 1), j] + k2 / 2 ) * (1 - sum(AA[j, ] * (NN[(t - 1), ]))) * step + lambda[p,j]
         k4 <- r[p,j] * (NN[(t - 1), j] + k3) * (1 - sum(AA[j, ] * (NN[(t - 1), ]))) * step + lambda[p,j]
         NN[t, j] <- NN[(t - 1), j] + 1/6 * (k1 + 2*k2 + 2*k3 + k4)
                       }
        result <- NN
        res[p, ] <- result[T,]
                    }
        }
 
 res_fil <- res[complete.cases(res), ]
 #write.csv(res_fil, file = "res_fil.csv", row.names = TRUE)

 # Convert to relative abundance
 res_ra <- res_fil / rowSums(res_fil)
 # Set a threshold of 1e-5
 res_ra[res_ra < 0.00001] <- 0
 #write.csv(res_ra, file = "res_ra.csv", row.names = TRUE)
 # Because a row in res_fil may become all zeros, computing relative abundance could yield NA; filter those out
 res_ra <- res_ra[complete.cases(res_ra), ]
 # Output all simulated communities
 res_ra_list[[a]] <- res_ra 
 
 ##Function to export a distance matrix to a data frame
 dist.to.df <- function(d)
         {
   size <- attr(d, "Size")
   return(
     data.frame(
       subset(expand.grid(row=2:size, col=1:(size-1)), row > col),
       distance=as.numeric(d),
       row.names = NULL  )
      )  }

  d <- vegdist(res_ra, "bray")
  d <- dist.to.df(d)
  Heterogeneity_mean <- mean(d[,3])
  Heterogeneity_vector[a] <- Heterogeneity_mean
  }

# Create data frames for each vector
mean_dist_df <- data.frame(mean_dist = mean_dist_vector)
cv_mean_df <- data.frame(cv_mean = cv_mean_vector)
heterogeneity_df <- data.frame(Heterogeneity = Heterogeneity_vector)

# Combine the data frames into a single data frame
result_df <- cbind(mean_dist_df, cv_mean_df, heterogeneity_df)

# file based on the value of dispersion
output_file <- paste0("output_", dispersion, ".csv")
write.csv(result_df, output_file, row.names = FALSE)

}


