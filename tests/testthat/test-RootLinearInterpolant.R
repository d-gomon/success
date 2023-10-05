test_that("Check whether interpolant does what it's supposed to", {
  #We only work with increasing data in y, so no need to check when y is not increasing
  x <- c(1, 2, 3, 4)
  y <- c(3, 5, 9, 20)
  expect_equal(RootLinearInterpolant(x, y, y0 = 4), 1.5)
  expect_equal(RootLinearInterpolant(x, y, y0 = 7), 2.5)
  #check for unsorted
  rord <- sample(1:4)
  x <- x[rord]
  y <- y[rord]
  expect_equal(RootLinearInterpolant(x, y, y0 = 4), 1.5)
  expect_equal(RootLinearInterpolant(x, y, y0 = 7), 2.5)
  #Return NA when no inverse exists
  expect_equal(RootLinearInterpolant(x, y, y0 = 30), NA)
})
