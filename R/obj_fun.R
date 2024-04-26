#' @title PCGLASSO objective function
#'
#' @description Evaluates the PCGLASSO objective function (negative penalised likelihood).
#'
#' @param Delta Partial correlation matrix Delta.
#' @param xi Vector of square root diagonals xi.
#' @param S Sample covariance matrix.
#' @param rho Penalty parameter.
#' @param c Diagonal parameter (double usual parameter value).
#'
#' @return Value of objective function evaluated at (Delta, xi).

obj_fun <- function(Delta, xi, S, rho, c){
  logdet <- - log(det(Delta))
  pen <- rho * sum(abs(Delta - diag(diag(Delta))))
  tr <- sum(S * (diag(xi) %*% Delta %*% diag(xi)))
  diags <- - c * sum(log(xi))
  logdet + pen + tr + diags
}
