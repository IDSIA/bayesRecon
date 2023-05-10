test_that("get_reconc_matrices produces expected aggregation and summing matrices", {
  aggrs  <- c(1,3,6,12)
  bott.f <- 12
  bott.h <- 12

  sort_aggrs <- sort(aggrs,decreasing=TRUE)
  expected_rowSumsS <- rep(sort_aggrs,bott.h/sort_aggrs)
  expectedLenRowSumsA <- sum(bott.h/sort_aggrs[-length(sort_aggrs)])
  expected_rowSumsA <- expected_rowSumsS[1:expectedLenRowSumsA]

  out <- get_reconc_matrices(aggrs,bott.f,bott.h)

  diff <- max(abs(expected_rowSumsA-rowSums(out$A)))+max(abs(expected_rowSumsS-rowSums(out$S)))
  expect_equal(diff, 0)
})
