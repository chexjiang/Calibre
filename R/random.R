#' Create deterministic, pair-specific random seeds
#'
#' The returned seeds do not depend on worker count or scheduling order.
#' @param n Number of pair-specific seeds to return.
#' @param seed Master seed.
#' @return An integer vector of length `n`.
#' @export
epilink_pair_seeds <- function(n, seed = 123L) {
  if (length(n) != 1L || !is.finite(n) || n < 0 || n != as.integer(n)) {
    stop("n must be one non-negative integer", call. = FALSE)
  }
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be one finite number", call. = FALSE)
  }
  n <- as.integer(n)
  if (!n) return(integer())
  modulus <- .Machine$integer.max - 1
  as.integer(((as.double(seed) - 1 + seq_len(n) - 1) %% modulus) + 1)
}

with_epilink_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}
