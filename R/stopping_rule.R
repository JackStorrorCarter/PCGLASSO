#' @title Stopping rule
#'
#' @description Evaluates a stopping rule for a Delta, xi update.
#'
#' @param Delta New approximated Delta.
#' @param Delta_old Previous approximated Delta.
#' @param xi New approximated xi.
#' @param xi_old Previous approximated xi.
#' @param obj New objective function value.
#' @param obj_old Previous objective function value.
#' @param threshold Threshold for stopping rule to be achieved.
#'
#' @return Indicator of stopping rule.

stopping_rule <- function(Delta, Delta_old, xi, xi_old,
                          obj, obj_old, threshold) {
  if( obj > obj_old ){
    return(FALSE)
  }
  Delta_change <- norm(Delta - Delta_old, type = '2') / norm(Delta_old, type = '2')
  xi_change <- norm( xi - xi_old, type = '2') / norm(xi_old, type = '2')
  obj_change <- abs(obj - obj_old)
  Delta_change + xi_change + obj_change < threshold
}
