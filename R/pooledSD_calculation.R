#' Pooled Standard Deviation Calculation
#'
#' This function allows you to calculate a pooled standard deviation from a vector of standard deviations.
#' @param sdvec vector of standard deviation values.
#' @keywords pooledSD
#' @export
#' @examples
#' pooledSD()

pooledSD <- function(sdvec) {
  sum <- 0
  for (sd in sdvec){
    sqr <- sd^2
    sum <- sum + sqr
  }
  sum <- sum / length(sdvec)
  pooledsd <- sqrt(sum)
  return(pooledsd)
}
