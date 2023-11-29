test_that("Test hierfamily", {
  S = matrix(
    data = c(1,1,1,1,1,1,
             1,1,0,0,0,0,
             0,0,1,1,0,0,
             0,0,0,0,1,1,
             1,0,0,0,0,0,
             0,1,0,0,0,0,
             0,0,1,0,0,0,
             0,0,0,1,0,0,
             0,0,0,0,1,0,
             0,0,0,0,0,1), ncol = 6, byrow = TRUE
  )
  Y = c(1,2,3,4,5,6,7,8,9,10)
  split_hierarchy.res = .split_hierarchy(S, Y)
  
  distr = list(
    "continuous",
    "continuous", "continuous","continuous",
    "continuous","continuous","continuous","continuous","continuous","continuous"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), 0)
  
  distr = list(
    "continuous",
    "discrete", "continuous","continuous",
    "continuous","continuous","continuous","continuous","continuous","continuous"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), -1)
  
  distr = list(
    "discrete",
    "continuous", "continuous","continuous",
    "continuous","continuous","continuous","discrete","continuous","continuous"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), -1)
  
  distr = list(
    "discrete",
    "continuous", "continuous","continuous",
    "discrete","discrete","discrete","discrete","discrete","discrete"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), 0)
  
  distr = list(
    "continuous",
    "continuous", "continuous","continuous",
    "discrete","discrete","discrete","discrete","discrete","discrete"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), 0)
  
  distr = list(
    "gaussian",
    "continuous", "continuous","continuous",
    "discrete","discrete","discrete","discrete","discrete","discrete"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), 0)
  
  distr = list(
    "gaussian",
    "gaussian", "gaussian","gaussian",
    "nbinom","poisson","poisson","nbinom","nbinom","nbinom"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), 0)
  
  distr = list(
    "gaussian",
    "gaussian", "gaussian","nbinom",
    "nbinom","poisson","poisson","nbinom","nbinom","nbinom"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), 0)
  
  distr = list(
    "gaussian",
    "gaussian", "gaussian","nbinom",
    "nbinom","poisson","poisson","nbinom","gaussian","nbinom"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), -1)
  
  distr = list(
    "gaussian",
    "gaussian", "gaussian","nbinom",
    "nbinom","poisson","poisson","nbinom","continuous","nbinom"
  )
  expect_equal(.check_hierfamily_rel(split_hierarchy.res, distr, debug = TRUE), -1)
})