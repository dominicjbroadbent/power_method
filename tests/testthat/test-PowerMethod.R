test_that("power.method works", {
  expect_equal(power.method(rbind(c(1,1),c(1,1)))$eigenvalues[2],2)
  expect_true(power.method(rbind(c(1,2),c(1,1)),dominant =FALSE)$eigenvalues[10] - -0.4142137 < 0.001)
})
