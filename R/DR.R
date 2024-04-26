#' @title Douglas Rachford algorithm
#'
#' @description Runs one step of the Douglas Rachford algorithm for updating Delta.
#'
#' @param Delta_1 Previous value to be updated.
#' @param Delta_2 Previous value to be updated.
#' @param Delta_3 Previous value to be updated.
#' @param S_aux Upper triangle of the matrix diag(xi) S diag(xi)
#' @param rho Penalty parameter.
#' @param UT Matrix indicating upper triangle of a pxp matrix.
#' @param LT Matrix indicating lower triangle of a pxp matrix.
#'
#' @return Updated value for Delta.

DR <- function(Delta_1, Delta_2, Delta_3, S_aux, rho, UT, LT) {
  Delta_3 <- Delta_3 + Delta_1 - Delta_2
  eig <- eigen(Delta_3, symmetric = TRUE)
  Delta_2 <- eig$vectors %*% diag(0.5 * (eig$values + sqrt(eig$values^2 + 4))) %*% t(eig$vectors)
  Delta_1[UT] <- shrink(((2 * Delta_2) - Delta_3)[UT] - S_aux,
                        rho
  )
  Delta_1[LT] <- t(Delta_1)[LT]
  list(Delta_1 = Delta_1, Delta_2 = Delta_2, Delta_3 = Delta_3)
}
