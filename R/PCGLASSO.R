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

pcglasso <- function(S, rho, c = NULL, Theta_start = NULL, threshold = 10^(-4), max_iter = 10000) {
  if (!is.matrix(S)) {
    stop("S is not a matrix")
  }
  if (!is.numeric(S)) {
    stop("S is not a numeric matrix")
  }
  if (!isSymmetric(S)) {
    stop("S is not a symmetric matrix")
  }
  p <- dim(S)[1]
  S_diags <- sqrt(diag(S))
  S <- cov2cor(S)
  S_eigen <- eigen(S, symmetric = TRUE)
  S_evals <- S_eigen$values
  S_evecs <- S_eigen$vectors
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
    if (identical(k, as.integer(0))) {
      if (c > 1) {
        stop("c must be less than or equal to 1")
      }
    } else {
      if (c >= 1 - k / p) {
        warning("c is too large - no solution exists")
      }
    }
  } else {
    if (identical(k, as.integer(0))){
      c <- 1
    } else{
      c <- 0.75 * (1 - k / p)
    }
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
    Theta_start <- S_diags * Theta_start * rep(S_diags, each = p)
    if (RSpectra::eigs_sym(Theta_start, k = 1, which = "SA", opts = list(retvec = FALSE))$values[1] < 1e-08) {
      stop("Theta_start is not positive definite")
    }
  } else {
    if (identical(k, as.integer(0))){
      Theta_start <- S_evecs %*% ( (1 / (S_evals)) * t(S_evecs) )
    } else {
      Theta_start <- S_evecs %*% ( (1 / (S_evals + 1 - min(S_evals))) * t(S_evecs) )
    }
  }

  Theta <- Theta_start
  Delta <- cov2cor(Theta_start)
  Delta_2 <- Delta
  Delta_3 <- Delta
  xi <- sqrt(diag(Theta_start))

  ind <- FALSE
  niter <- 1
  UT <- upper.tri(Delta)
  LT <- lower.tri(Delta)

  while (ind == FALSE) {
    Theta_old <- Theta
    Delta_old <- Delta
    xi_old <- xi

    gamma_fb <- 0.9 / RSpectra::eigs_sym(Delta * S, k = 1, which = "LA", opts = list(retvec = FALSE))$values[1]
    fb_param2 <- 4 * c * gamma_fb
    ind2 <- FALSE
    val <- sum(S * Theta) - c * sum(log(xi))
    threshold_FB <- max(10^(-3) * 0.9^(niter - 1), threshold * 0.1)
    while(ind2 == FALSE){
      xi_old1 <- xi
      fb_param1 <- gamma_fb * 2 * colSums(S * Delta * xi)
      xi <- FB(xi = xi, fb_param1, fb_param2)
      ind2 <- ( sum(abs(xi - xi_old1)) < sum(abs(xi_old1)) * threshold_FB )
      if(ind2){
        ind2 <- sum(S * (xi * Delta * rep(xi, each=p))) - c * sum(log(xi)) < val + 1e-08
      }
    }

    ind1 <- FALSE
    S_aux <- xi * S * rep(xi, each=p)
    val <- -log(det(Delta)) + sum(S_aux * Delta) + rho * (sum(abs(Delta)) - p)
    threshold_DR <- max(10^(-3) * 0.9^(niter - 1), threshold * 0.1)
    while(ind1 == FALSE){
      Delta_old1 <- Delta
      DR_out <- DR(Delta, Delta_2, Delta_3, S_aux[UT], rho, UT, LT)
      Delta <- DR_out$Delta_1
      Delta_2 <- DR_out$Delta_2
      Delta_3 <- DR_out$Delta_3
      ind1 <- ( sum(abs(Delta - Delta_old1)) <= (max(sum(abs(Delta_old1)) - p, 1e-08) * threshold_DR ) )
      if(ind1){
        ind1 <- RSpectra::eigs_sym(Delta, k = 1, which = "SA", opts = list(retvec = FALSE))$values[1] > 1e-08
      }
      if(ind1){
        ind1 <- -log(det(Delta)) + sum(S_aux * Delta) + rho * (sum(abs(Delta)) - p) < val + 1e-08
      }
    }

    if (niter == max_iter) {
      ind <- TRUE
      warning("Maximum number of iterations reached")
    } else{
      niter <- niter + 1
      Theta <- xi * Delta * rep(xi, each=p)
      ind <- sum(abs(Delta - Delta_old)) / max(sum(abs(Delta_old)) - p, 1e-08) + sum(abs(xi - xi_old)) / sum(abs(xi_old)) < threshold
    }
  }
  (1/S_diags) * Theta * rep(1/S_diags, each = p)
}
