# Generate the weekly base forecasts parameters
# CHANGE THE WORKING DIRECTORY BEFORE RUNNING
rm(list=ls())
library(bayesRecon)

# Generate A matrix for weekly hierarchy
A <- .gen_weekly()
# Save matrix to file
write.table(A,file="./Weekly-Gaussian_A.csv",row.names = FALSE,sep=',',col.names = FALSE,quote = FALSE)

# Generate randomly bottom means
set.seed(1)
bottom_means <- runif(52,min=5,max=10)

# The upper means are generated by 
# aggregating the bottom means (according to A) with an epsilon of incoherence
eps = 0.1
upper_means <- (1+eps)*A%*%bottom_means
upper_means <- as.vector(upper_means)

# Create the matrix to save to csv
base_forecasts_out <- cbind(c(upper_means,bottom_means),c(rep(3,length(upper_means)),rep(2,length(bottom_means))))
write.table(base_forecasts_out,file="./Weekly-Gaussian_basef.csv",row.names = FALSE,sep=',',col.names = FALSE,quote = FALSE)


