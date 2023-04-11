source("hierarchy.R")


#-------------------------------------------------------------------------------

# Example 1

a1 <- c(1,1,0,0)
a2 <- c(0,0,1,1)
a3 <- c(1,1,1,0)
a4 <- c(1,1,1,1)

A <- rbind(a1, a2, a3, a4)

ind <- get_hier_rows(A)  #indices of the "hierarchy rows" of A
print(ind)


#-------------------------------------------------------------------------------
# Example 2 (monthly)

Am <- gen_monthly()

ind_monthly <- get_hier_rows(Am)

View(Am[as.logical(ind_monthly),])  #corresponds to H, but unsorted!


#-------------------------------------------------------------------------------
# Example 3 (weekly)


Aw <- gen_weekly()

k <- nrow(Aw)
set.seed(100)
v <- rnorm(k) # e.g. vector of parameters of the upper variables

l <- get_H(Aw, v)

H   <- l[[1]]
v_h <- l[[2]]
G   <- l[[3]]
v_g <- l[[4]]



#Try again, but first shuffle rows of Aw

rand <- sample(k)
Aw2 <- Aw[rand,]
v2 <- v[rand]

l2 <- get_H(Aw, v)

H2   <- l2[[1]]
v_h2 <- l2[[2]]
G2   <- l2[[3]]
v_g2 <- l2[[4]]

# check that they are the same
H == H2
G == G2
v_h == v_h2
v_g == v_g2

