# Regression/verification tests for the native C++ replacements of what used to be R
# callbacks to mvtnorm::pmvnorm(), MASS::mvrnorm(), and Matrix::nearPD() (see src/neuron.cpp's
# comments on bvn_cdf_uppertail(), dg_mvn_draw_cholesky(), and makePositiveDefinite()). mvtnorm,
# MASS, and Matrix are Suggests-only now (used here as reference implementations, not runtime
# dependencies), so every test below skips gracefully if the relevant one isn't installed.

test_that("bvn_cdf_uppertail (Owen's T) matches mvtnorm::pmvnorm", {
  skip_if_not_installed("mvtnorm")
  set.seed(1)
  for (i in 1:200) {
    h <- runif(1, -3, 3); if (abs(h) < 0.05) h <- 0.05
    k <- runif(1, -3, 3); if (abs(k) < 0.05) k <- 0.05
    rho <- runif(1, -0.9, 0.9)
    native <- bvn_cdf_uppertail_test(h, k, rho)
    ref <- mvtnorm::pmvnorm(
      lower = c(h, k), sigma = matrix(c(1, rho, rho, 1), 2, 2), keepAttr = FALSE
    )
    expect_equal(native, as.numeric(ref), tolerance = 1e-8)
  }
})

test_that("makePositiveDefinite (native nearPD) matches Matrix::nearPD", {
  skip_if_not_installed("Matrix")
  set.seed(2)
  for (trial in 1:20) {
    n <- sample(3:15, 1)
    M <- matrix(runif(n * n, -1, 1), n, n)
    M <- (M + t(M)) / 2
    diag(M) <- 1
    native <- makePositiveDefinite_test(M, TRUE, TRUE, TRUE)
    ref <- as.matrix(Matrix::nearPD(M, corr = TRUE, keepDiag = TRUE, do2eigen = TRUE)$mat)
    expect_equal(native, ref, tolerance = 1e-6, ignore_attr = TRUE)
  }
})

test_that("makePositiveDefinite preserves a unit diagonal on a Toeplitz matrix like dg_simulation() builds", {
  # Regression guard: the crude eigenvalue-clipping version this replaced never restored the
  # diagonal to 1 after clipping, silently biasing every downstream DG-simulated spike rate.
  skip_if_not_installed("Matrix")
  tp <- toeplitz(seq(1, 0.1, length.out = 20))
  native <- makePositiveDefinite_test(tp, TRUE, TRUE, TRUE)
  expect_equal(diag(native), rep(1, 20), tolerance = 1e-10)

  ref <- as.matrix(Matrix::nearPD(tp, corr = TRUE, keepDiag = TRUE, do2eigen = TRUE)$mat)
  expect_equal(native, ref, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("dg_mvn_draw_cholesky() draws match the target covariance (large-sample)", {
  Sigma <- matrix(c(1, 0.6, 0.3, 0.6, 1, 0.5, 0.3, 0.5, 1), 3, 3)
  L <- t(chol(Sigma)) # R's chol() returns upper-tri R s.t. Sigma = R'R; want lower-tri L s.t. Sigma = LL'
  expect_equal(L %*% t(L), Sigma, tolerance = 1e-10, ignore_attr = TRUE)

  draws <- dg_mvn_draw_cholesky_test(L, 200000, 42)
  expect_equal(dim(draws), c(200000, 3))
  expect_equal(colMeans(draws), c(0, 0, 0), tolerance = 0.02)
  expect_equal(cov(draws), Sigma, tolerance = 0.01, ignore_attr = TRUE)
})

test_that("dg_mvn_draw_cholesky() is deterministic given the same seed, and differs otherwise", {
  L <- matrix(c(1, 0.5, 0, 0.866), 2, 2) # lower-tri Cholesky factor of a 2x2 example Sigma
  d1 <- dg_mvn_draw_cholesky_test(L, 10, 7)
  d2 <- dg_mvn_draw_cholesky_test(L, 10, 7)
  d3 <- dg_mvn_draw_cholesky_test(L, 10, 8)
  expect_identical(d1, d2)
  expect_false(isTRUE(all.equal(d1, d3)))
})
