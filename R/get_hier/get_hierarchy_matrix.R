library(lpSolve)


#-------------------------------------------------------------------------------
# FUNCTION TO FIND THE ROWS OF A
# Returns a vector with length equal to the number of rows of A
# Each entry is 1 if the corresponding row has to be picked, and 0 otherwise

get_hier_rows <- function(A) {
  
  k <- nrow(A)
  m <- ncol(A)
  
  # matrix C of the coefficients of the non-linear problem 
  # (k variables, 1 constraint)
  
  C <- A %*% t(A) 
  
  for (i in 1:k) {
    for (j in 1:k) {
      C[i,j] <- C[i,j] * sum((A[i,] - A[j,]) * (A[i,] - A[j,] - 1)) *
        sum((A[j,] - A[i,]) * (A[j,] - A[i,] - 1))
    }
  }
  
  #-------------------
  # LINEARIZED PROBLEM
  # Number of variables: k + k^2
  # Number of constraints: 1 + N^2 + N^2 + N^2
  
  
  # Set coefficients of the objective function
  f.obj <- c(rep(1, k), rep(0, k^2))
  
  
  # Set matrix corresponding to coefficients of constraints by rows
  
  coeff <- c(rep(0, k), as.vector(C)) #first constraint
  
  M1 <- matrix(0, k^2, k + k^2)
  for (i in 1:k) {
    temp <- matrix(0, k, k)
    temp[i,] <- rep(-1, k)
    M1[,i] <- as.vector(temp)
  }
  M1[, (k+1):(k+k^2)] <- diag(1, k^2)
  
  M2 <- matrix(0, k^2, k + k^2)
  for (i in 1:k) {
    temp <- matrix(0, k, k)
    temp[,i] <- rep(-1, k)
    M2[,i] <- as.vector(temp)
  }
  M2[, (k+1):(k+k^2)] <- diag(1, k^2)
  
  M3 <- matrix(0, k^2, k + k^2)
  for (i in 1:k) {
    temp <- matrix(0, k, k)
    temp[,i] <- rep(-1, k)
    M3[((i-1)*k + 1) : (i*k),] <- temp - diag(1, k)
  }
  M3[, (k+1):(k+k^2)] <- diag(1, k^2)
  
  f.con <- rbind(coeff, M1, M2, M3)
  
  
  # Set unequality/equality signs
  f.dir <- c("=", rep("<=", k^2), rep("<=", k^2), rep(">=", k^2))
  
  
  # Set right hand side coefficients
  f.rhs <- c(rep(0,1 + 2*k^2), rep(-1, k^2))
  
  
  
  #---------------------
  # Solve the LP problem
  
  # Variables final values
  indices_sol <- lp("max", f.obj, f.con, f.dir, f.rhs, all.bin = TRUE)$solution[1:k]
  
  return(indices_sol)
}


#-------------------------------------------------------------------------------

# Function that extract the "hierarchy rows" from A, and sorts them in the 
# correct order (i.e. bottom-up)
# Also sorts accordingly the vector v (e.g. of parameters) 

get_H <- function(A, v) {
  
  #get the indices of the "hierarchy rows" of A
  indices_sol <- get_hier_rows(A)
  
  #extract rows from A
  ind_h <- as.logical(indices_sol)
  H <- A[ind_h,]
  v_h <- v[ind_h]
  
  #sort bottom-up
  ord <- order(rowSums(H))
  H <- H[ord,]
  v_h <- v[ord]
  
  #collect remaining rows in matrix G
  ind_g <- as.logical(1 - indices_sol)
  G <- A[ind_g,]
  v_g <- v[ind_g]
  
  return(list(H, v_h, G, v_g))
  
}


#-------------------------------------------------------------------------------

# Functions to generate the monthly and weekly A matrices 


gen_monthly <- function() {
  
  H <- matrix(0, nrow=10, ncol=12)
  for (j in 1:6) {
    H[j, (2*(j-1)+1):(2*j)] <- 1
  }
  for (j in 1:3) {
    H[6+j, (4*(j-1)+1):(4*j)] <- 1
  }
  H[6+3+1,] <- 1
  
  G <- matrix(0, nrow=6, ncol=12)
  for (j in 1:4) {
    G[j, (3*(j-1)+1):(3*j)] <- 1
  }
  for (j in 1:2) {
    G[4+j, (6*(j-1)+1):(6*j)] <- 1
  }
  
  return(rbind(H,G))
  
}


gen_weekly <- function() {
  
  H <- matrix(0, nrow=40, ncol=52)
  for (j in 1:26) {
    H[j, (2*(j-1)+1):(2*j)] <- 1
  }
  for (j in 1:13) {
    H[26+j, (4*(j-1)+1):(4*j)] <- 1
  }
  H[26+13+1,] <- 1
  
  G <- matrix(0, nrow=6, ncol=52)
  for (j in 1:4) {
    G[j, (13*(j-1)+1):(13*j)] <- 1
  }
  for (j in 1:2) {
    G[4+j, (26*(j-1)+1):(26*j)] <- 1
  }
  
  return(rbind(H,G))
  
}




