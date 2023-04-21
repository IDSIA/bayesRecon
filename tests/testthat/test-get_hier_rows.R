test_that(".get_hier_rows works on a monthly hierarchy", {
  month_hier <- .gen_monthly()
  ind_best <- .get_hier_rows(month_hier)

  best_obj_fun <- max(colSums(month_hier[as.logical(ind_best),]))-sum(ind_best)

  expect_equal(best_obj_fun, -7)
})


test_that(".get_hier_rows works on a weekly hierarchy", {
  week_hier <- .gen_weekly()
  ind_best <- .get_hier_rows(week_hier)

  best_obj_fun <- max(colSums(week_hier[as.logical(ind_best),]))-sum(ind_best)

  expect_equal(best_obj_fun, -37)
})
