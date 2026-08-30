#' Function for conditional sampling - sample the CDE
#'
#' @param cde.obj An hdrcde::cde object
#' @param Ytest Predicted value of testing data
#' @import dplyr
#' @export

CDE_sample <-
function(cde.obj, Ytest){
  cde_sample <- cde.obj$z %>% 
    apply(1, cumsum)
  # remove all 0 columns
  non0col <- which(colSums(cde_sample != 0) > 0)
  cde_sample<- cde_sample[, non0col] 
  cde_sample <- apply(cde_sample, 2, function(y){y/max(y)})
  cde_list <- list(x_grid = cde.obj$x[non0col], 
              y_grid = cde.obj$y,
              cde_sample = cde_sample)
  
  Ytest_close <- sapply(Ytest, 
                       function(i_y) which.min(abs(i_y - cde_list$x_grid)))
  u <- runif(n = length(Ytest))
  Y_sample <- cde_list$y_grid[
    sapply(seq_along(Ytest), function(i_obs){
      which(u[i_obs] < cde_list$cde_sample[, Ytest_close[i_obs]])[1]
    })]
  return(Y_sample)
}

#' Function for conditional sampling - create cde object
#'
#' @param Yhat True value of calibration data
#' @param resid Predicted value of calibration data
#' @param method_param The argument `method` for function `hdrcde::cde.bandwidths`
#' @param xden_param The argument `xden` for function `hdrcde::cde.bandwidths`
#' @param xmar Number of bins for conditional density estimation
#' @param ymar Number of bins for conditional density estimation
#' @param plot whether to generate the cde plot. Default is FALSE
#' @importFrom hdrcde cde cde.bandwidths
#' @import dplyr
#' @export

fit_cde <- function(Yhat, resid, method_param=2, xden_param="uniform", xmar=50, ymar=100, plot="FALSE"){
  # calculate a,b
  ab <- cde.bandwidths(
    x = Yhat, y = resid, method=method_param, xden=xden_param
  )
  a <- ab$a
  b <- ab$b
  # use hdr cde
  fit_cde <- hdrcde::cde(
    x = Yhat, y = resid,
    nxmargin = xmar,
    nymargin = ymar,
    a=a,
    b=b
  )
  if(plot){plot(fit_cde)}
  return(fit_cde)
}
