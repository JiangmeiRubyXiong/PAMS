#' Function for conditional sampling
#'
#' @param Ytrue True value of calibration data
#' @param Yhat Predicted value of calibration data
#' @param Ytest Predicted value of testing data
#' @param xmar Number of bins for conditional density estimation
#' @param ymar Number of bins for conditional density estimation
#' @export

CDE <-
function(Ytrue, Yhat, Ytest, xmar=50, ymar=100){
  fit_cde <- hdrcde::cde(
    x = Yhat, y = Ytrue,
    nxmargin = xmar,
    nymargin = ymar
  )
  cde_sample <- fit_cde$z %>% 
    apply(1, cumsum)
  cde_sample <- cde_sample / max(cde_sample)
  cde_list <- list(x_grid = fit_cde$x, 
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
