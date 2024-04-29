#' @title PCGLASSO
#'
#' @description Numerical solution to the PCGLASSO.
#'
#' @param S Sample covariance matrix.
#' @param rho Penalty parameter.
#' @param c Diagonal parameter. Default is the largest possible value for which
#'   the solution exists.
#' @param Theta_start Starting value of Theta. Default is the inverse of S
#'   (generalised inverse if S is not positive definite).
#' @param threshold Threshold for stopping rule of algorithm.
#' @param max_iter Maximum number of iterations.
#'
#' @return Numerical solution Theta to the PCGLASSO optimisation problem.
#' @export

pcglasso <- function(S, rho, c = NULL, Theta_start = NULL, threshold = 10^(-5), max_iter = 10000) {
  if (!is.matrix(S)) {
    stop("S is not a matrix")
  }
  if (!is.numeric(S)) {
    stop("S is not a numeric matrix")
  }
  if (!isSymmetric(S)) {
    stop("S is not a symmetric matrix")
  }
  S_evals <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (min(S_evals) < -1e-08) {
    stop("S is not positive semidefinite")
  }
  k <- length(which(S_evals < 1e-08))
  if (!is.null(c)) {
    if (length(c) != 1 || !is.numeric(c)) {
      stop("Not a valid penalty parameter c")
    }
    if (c <= 0) {
      stop("c must be greater than 0")
    }
    if (identical(k, 0)) {
      if (c > 1) {
        stop("c must be less than or equal to 1")
      }
    } else {
      if (c > 1 / (k + 1)) {
        stop("c is too large - no solution exists")
      }
    }
  } else {
    c <- 1 / (k + 1)
  }
  if (length(rho) != 1 || !is.numeric(rho)) {
    stop("Not a valid penalty parameter rho")
  }
  if (rho < 0) {
    stop("rho must be greater than or equal to 0")
  }
  if (length(threshold) != 1 || !is.numeric(threshold) || threshold < 0) {
    stop("Not a valid threshold")
  }
  if (length(max_iter) != 1 || !is.numeric(max_iter) || !identical(max_iter, round(max_iter))) {
    stop("Not a valid number of maximum iterations")
  }

  p <- dim(S)[1]
  c <- 2 * c

  if (!is.null(Theta_start)){
    if (!is.matrix(Theta_start)) {
      stop("Theta_start is not a matrix")
    }
    if (!is.numeric(Theta_start)) {
      stop("Theta_start is not a numeric matrix")
    }
    if (!isSymmetric(Theta_start)) {
      stop("Theta_start is not a symmetric matrix")
    }
    if (!isTRUE(all.equal(dim(Theta_start), c(p, p)))) {
      stop("Dimensions of S and Theta_start do not match")
    }
    if (min(eigen(Theta_start, symmetric = TRUE, only.values = TRUE)$values) < 1e-08) {
      stop("Theta_start is not positive definite")
    }
  } else {
    if (min(S_evals) > 1e-08){
      Theta_start <- solve(S)
    } else {
      Theta_start <- solve(S + diag(1e-08,p))
    }
  }

  Delta <- cov2cor(Theta_start)
  Delta_2 <- Delta
  Delta_3 <- Delta
  xi <- sqrt(diag(Theta_start))

  F_val <- rep(0, max_iter+1)
  F_val[1] <- obj_fun(Delta, xi, S, rho, c)
  ind <- FALSE
  niter <- 1
  UT <- upper.tri(Delta)
  LT <- lower.tri(Delta)

  while (ind == FALSE) {
    Delta_old <- Delta
    xi_old <- xi

    ind1 <- FALSE
    S_aux <- (diag(xi) %*% S %*% diag(xi))[UT]
    while(ind1 == FALSE){
      Delta_old2 <- Delta
      DR_out <- DR(Delta, Delta_2, Delta_3, S_aux, rho, UT, LT)
      Delta <- DR_out$Delta_1
      Delta_2 <- DR_out$Delta_2
      Delta_3 <- DR_out$Delta_3
      ind1 <- (norm(Delta - Delta_old2, type = '2') / norm(Delta_old2, type = '2') < threshold/10)
      if(ind1){
        ind1 <- min(eigen(Delta, symmetric = TRUE, only.values = TRUE)$values) > 1e-08
      }
    }

    gamma_fb <- 0.9 / max(abs(eigen(Delta * S, symmetric = TRUE, only.values = TRUE)$values))
    fb_param2 <- 4 * c * gamma_fb
    ind2 <- FALSE
    while(ind2 == FALSE){
      xi_old2 <- xi
      fb_param1 <- gamma_fb * 2 * colSums(S * Delta * xi)
      xi <- FB(xi = xi, fb_param1, fb_param2)
      ind2 <- (norm(xi - xi_old2, type = '2') / norm(xi_old2, type = '2') < threshold/10)
    }

    if (niter == max_iter) {
      ind <- TRUE
      warning("Maximum number of iterations reached")
    }
    niter <- niter + 1
    F_val[niter] <- obj_fun(Delta, xi, S, rho, c)
    if (stopping_rule(Delta, Delta_old, xi, xi_old, F_val[niter], F_val[niter - 1], threshold)) {
      ind <- TRUE
    }
  }
  Theta <- diag(xi) %*% Delta %*% diag(xi)
  Theta
}
