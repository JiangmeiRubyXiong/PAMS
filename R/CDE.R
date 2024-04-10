#' Function for conditional sampling
#'
#' @param Yhat True value of calibration data
#' @param resid Predicted value of calibration data
#' @param Ytest Predicted value of testing data
#' @param xmar Number of bins for conditional density estimation
#' @param ymar Number of bins for conditional density estimation
#' @importFrom hdrcde cde cde.bandwidths
#' @export

CDE <-
function(Yhat, resid, Ytest, xmar=50, ymar=100){
  # calculate a,b
  ab <- cde.bandwidths(
    x = Yhat, y = resid, method=2, xden="uniform"
  )
  a <- ab$a
  b <- ab$b
  # a=0.01
  # b=0.01
  # use hdr cde
  fit_cde <- hdrcde::cde(
    x = Yhat, y = resid,
    nxmargin = xmar,
    nymargin = ymar,
    a=a,
    b=b
  )
  cde_sample <- fit_cde$z %>% 
    apply(1, cumsum)
  # remove all 0 columns
  non0col <- which(colSums(cde_sample != 0) > 0)
  cde_sample<- cde_sample[, non0col] 
  cde_sample <- apply(cde_sample, 2, function(y){y/max(y)})
  cde_list <- list(x_grid = fit_cde$x[non0col], 
              y_grid = fit_cde$y,
              cde_sample = cde_sample)
  
  Ytest_close <- sapply(Ytest, 
                       function(i_y) order(abs(i_y - cde_list$x_grid))[1])
  u <- runif(n = length(Ytest))
  Y_sample <- cde_list$y_grid[
    sapply(seq_along(Ytest), function(i_obs){
      which(u[i_obs] < cde_list$cde_sample[, Ytest_close[i_obs]])[1]
    })]
  return(Y_sample)
}
