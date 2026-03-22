#' Safe division
#'
#' @param x Numerator.
#' @param denom Denominator.
#' @param eps Small constant.
#'
#' @return Numeric vector.
#' @noRd
safe_div <- function(x, denom, eps = 1e-12) {
  x / (denom + eps)
}

#' Compute percentile-like rank
#'
#' Larger values get larger percentile scores.
#'
#' @param x Numeric vector.
#'
#' @return A numeric vector ranging from 0 to 1.
#' @noRd
percentile_rank_desc <- function(x) {
  if (length(x) <= 1) {
    return(rep(1, length(x)))
  }
  1 - (data.table::frank(-x, ties.method = "average") - 1) / (length(x) - 1)
}

#' Assign scRNA tier from q_sc
#'
#' @param q_sc Composite percentile.
#' @param high Threshold for high tier.
#' @param mid Threshold for mid tier.
#'
#' @return Character vector.
#' @noRd
assign_scrna_tier <- function(q_sc, high = 0.8, mid = 0.5) {
  data.table::fifelse(
    q_sc >= high, "high",
    data.table::fifelse(q_sc >= mid, "mid", "low")
  )
}

#' Jaccard similarity
#'
#' @param a Character vector.
#' @param b Character vector.
#'
#' @return Numeric scalar.
#' @export
jaccard_similarity <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(1)
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}
