test_that("Entrytime vs Entrytime + unit & inputs",{
  expect_length(arrival_rate(surgerydat), length(unique(surgerydat$unit)))
  expect_length(arrival_rate(surgerydat[, "entrytime", drop = FALSE]), 1)
  expect_warning(arrival_rate(surgerydat[, "age"]), "Provided vector is not named, assuming that specified data represents entrytimes.")
  expect_error(arrival_rate(surgerydat[, "age", drop = FALSE]), "Argument data should contain at least a numeric column named 'entrytime'")
})






test_that("Arrival rate larger when more people", {
  arr_surg <- surgerydat[, c("entrytime", "age")]
  expect_lt(arrival_rate(subset(arr_surg, age < 60)), arrival_rate(arr_surg))
})
