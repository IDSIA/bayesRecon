test_that("get_reconc_matrices produces expected aggregation and summing matrices", {
  aggrs  <- c(1,3,6,12)
  h <- 12

  sort_aggrs <- sort(aggrs,decreasing=TRUE)
  expected_rowSumsS <- rep(sort_aggrs,h/sort_aggrs)
  expectedLenRowSumsA <- sum(h/sort_aggrs[-length(sort_aggrs)])
  expected_rowSumsA <- expected_rowSumsS[1:expectedLenRowSumsA]

  out <- get_reconc_matrices(aggrs,h)

  diff <- max(abs(expected_rowSumsA-rowSums(out$A)))+max(abs(expected_rowSumsS-rowSums(out$S)))
  expect_equal(diff, 0)
})

test_that(".get_Au and .lowest_lev produce expected outcomes", {

  A <- matrix(
    data=c(1, 1, 1, 1, 1, 1,
           1, 1, 0, 0, 0, 0,
           0, 0, 1, 1, 0, 0,
           0, 0, 0, 0, 1, 1),nrow=4,byrow = TRUE)
  
  
  expect_equal(.get_Au(A),matrix(c(1,1,1),ncol=3))
  expect_equal(.lowest_lev(A),c(2,3,4))
  
})

test_that(".get_Au behaves identically on A with different row order.",{
  A <- matrix(data=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                     1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
              nrow=10,byrow = TRUE)
  
  A1 <- matrix(data=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                      1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1),
               nrow=10,byrow = TRUE)
  
  expect_equal(.get_Au(A),.get_Au(A1))
})
