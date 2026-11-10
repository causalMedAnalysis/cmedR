# Small simulated data set with known mediation structure:
#   C ~ N(0, 1)
#   D ~ Bernoulli(0.5), independent of C
#   M = 0.5 * D + 0.3 * C + e_M
#   Y = 0.4 * M + 0.2 * D + 0.1 * C + e_Y
# Under this model (no exposure-mediator interaction):
#   NDE = 0.2, NIE = 0.5 * 0.4 = 0.2, ATE = 0.4
sim_mediation_data <- function(n = 3000, seed = 20260101) {
  withr_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  }
  on.exit({
    if (!is.null(withr_seed)) {
      assign(".Random.seed", withr_seed, envir = globalenv())
    }
  })
  set.seed(seed)
  C <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  M <- 0.5 * D + 0.3 * C + rnorm(n)
  Y <- 0.4 * M + 0.2 * D + 0.1 * C + rnorm(n)
  data.frame(Y = Y, D = D, M = M, C = C)
}
