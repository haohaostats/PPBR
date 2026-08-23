test_that("NeuSTART analysis matches locked reference values", {
  design <- ppbr_design(c(1, 3, 6, 8, 10), 0.10, 3, 33, 0.20)
  fit <- ppbr(design, c(3, 10, 12, 8, 0), c(0, 0, 2, 0, 0), 8)
  expect_equal(fit$final_dose, 8)
  expect_equal(fit$targets$probability,
               c(0.00007885037464, 0.0009807559411, 0.01703417899,
                 0.05891449436, 0.3976018854, 0.5253898349),
               tolerance = 1e-8)
  expect_equal(fit$brackets$lower_favored[4], 0.9368, tolerance = 1e-4)
  expect_false(ppbr_next(fit)$continue)
})

test_that("an empty trial starts at the configured dose", {
  design <- ppbr_design(1:4, 0.25, 3, 24, start_dose = 2)
  fit <- ppbr(design, integer(4), integer(4))
  decision <- ppbr_next(fit)
  expect_true(decision$continue)
  expect_equal(decision$next_dose, 2)
})

test_that("invalid dose-level data are rejected", {
  design <- ppbr_design(1:3, 0.25, 3, 18)
  expect_error(ppbr(design, c(3, 3, 0), c(0, 4, 0)), "no greater")
  expect_error(ppbr(design, c(3, 0), c(0, 0)), "3 non-negative")
})

test_that("simulation accepts user-specified monotone scenarios", {
  design <- ppbr_design(1:3, 0.25, 3, 12)
  sim <- ppbr_simulate(design, c(0.10, 0.25, 0.45), nsim = 5, seed = 1)
  expect_s3_class(sim, "ppbr_simulation")
  expect_equal(nrow(sim$trials), 5)
  expect_equal(sum(sim$selection$probability), 1)
  expect_error(ppbr_simulate(design, c(0.2, 0.1, 0.3), 2), "nondecreasing")
})

test_that("simulation handles dose-boundary candidates without warnings", {
  design <- ppbr_design(1:6, 0.10, 3, 60)
  scenarios <- list(
    lower_target = c(0.05, 0.10, 0.18, 0.27, 0.38, 0.50),
    upper_target = c(0.01, 0.03, 0.05, 0.07, 0.09, 0.10)
  )

  result <- expect_no_warning(ppbr_simulate(design, scenarios, nsim = 3, seed = 11))
  expect_s3_class(result, "ppbr_simulation")
  expect_equal(nrow(summary(result)), 2L)
})
