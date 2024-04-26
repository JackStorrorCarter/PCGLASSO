#' @title Forward backward algorithm
#'
#' @description Runs one step of the forward backward algorithm for updating xi.
#'
#' @param xi Previous value to be updated.
#' @param fb_param1 First parameter equal to 2 gamma_fb colSums(S Delta xi)
#' @param fb_param2 Second parameter equal to 4 c gamma_fb
#'
#' @return Updated value for xi.

FB <- function(xi, fb_param1, fb_param2) {
  xi <- 0.5 * ((xi - fb_param1) + sqrt((xi - fb_param1)^2 + fb_param2))
  xi
}
