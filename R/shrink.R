#' @title Shrink operator
#'
#' @description Performs the shrink operation.
#'
#' @param a Value or vector of number(s) to be shrunk.
#' @param b Amount to be shrunk by.
#'
#' @return Value or vector of shrunk number(s).

shrink <- function(a, b) {
  output <- (abs(a) <= b)*0 + (abs(a) > b)*(a - b*sign(a))
  output
}
