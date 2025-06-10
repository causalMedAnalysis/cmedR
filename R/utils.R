utils::globalVariables(c(".data"))

demean <- function(x, w = rep(1, length(x))) x - weighted.mean(x, w, na.rm = TRUE)

`%||%` <- function(a, b) if (!is.na(a)) a else b

# trimming functions for use with IPW estimators
trim <- function(x, min = 0.01, max = 1) {
  x[x<min] <- min
  x[x>max] <- max
  x
}

trimQ <- function(x, low = 0.01, high = 0.99) {
  min <- quantile(x, low)
  max <- quantile(x, high)

  x[x<min] <- min
  x[x>max] <- max
  x
}


# function to combine lists of vectors with the same structure
# (used in the parallelized bootstraps)
comb_list_vec <- function(...) {
  mapply(c, ..., SIMPLIFY = FALSE)
}

