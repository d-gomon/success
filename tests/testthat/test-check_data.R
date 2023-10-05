test_that("Input checks", {
  matcheck <- matrix(rnorm(12), ncol = 3)
  colnames(matcheck) <- c("entrytime", "survtime", "censorid")
  expect_warning(check_data(matcheck), "Provided data is not a data frame, attempting to convert.")
  dfcheck <- data.frame(survtime = rep(3, 5), censorid = rep(1, 5))
  expect_error(check_data(dfcheck), "Entry time missing for subjects. Please specify them as named column
        'entrytime' in your data frame.")
  colnames(dfcheck)[1] <- "entrytime"
  expect_error(check_data(dfcheck), "Survival time missing for subjects. Please specify them as named
          column 'survtime' in your data frame.")
  colnames(dfcheck)[2] <- "survtime"
  expect_warning(check_data(dfcheck), "No censoring mechanism specified. Assuming data is uncensored.
            May lead to an increased amount of signals!")
  suppressWarnings(checked_data <- check_data(dfcheck))
  #Expect a column of 1's to be created.
  expect_true(all(checked_data$censorid == 1))
  checked_data[1, 1] <- NA
  expect_error(check_data(checked_data), "Please make sure 'data' has no missing values.")
})
