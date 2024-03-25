#' Refines image imputation with multiple imputation and conditional density estimation
#' 
#' @param yhat Predicted value for testing data.
#' @param calibYhat Predicted value for calibration data. For method `bySample`, it should be a vector of each training slide predicted value predicted by the other training slides.
#' @param calibY True value for calibration data.
#' @param cdeType Different types of conditional density estimation with calibration data, conditioning on the predicted value of predicted value of testing data.
#' @param sampleID (optional) Sample ID of each cell, only need to be provided for method `bySample`. It should be a list containing two vectors: the first for `calibY`, the second for `Yhat`.
#' @param n.mi (optional) Number of multiple imputations. Default to be 10.
#' @param x.mar (optional) Number of bins for conditional density estimation. Default to be 50.
#' @param y.mar (optional) umber of bins for conditional density estimation. Default to be 1000.
#' @importFrom hdrcde cde
#' @importFrom stats rnorm runif
#' @import dplyr
#' @return
#' \item{resmat}{The matrix of multiple imputations. The number of rows is equal to the number of cells, and the number of columns is the number of multiple imputation.}
#'
#' @examples
#'# create mock data for different types of conditional density estimation
#'## method `simple` and `local`
#'yhat_sl <- rnorm(1000, sd=20)
#'calibY_sl <- runif(50, min=-5, max=5)
#'calibYhat_sl <- runif(50, min=-5, max=5)
#'test.imp.fun <- imimp(yhat_sl, calibYhat_sl, calibY_sl, cdeType = "simple", x.mar = 50, y.mar=100)
#'test.imp.fun <- imimp(yhat_sl, calibYhat_sl, calibY_sl, cdeType = "local", x.mar = 50, y.mar=100)
#'## method `bySample`
#'yhat_b <- rnorm(600)
#'calibY_b <- runif(700, min=-5, max=5)
#'calibYhat_b <- runif(700, min=-5, max=5)
#'sample_b <- list(rep(state.name[7:16], each=70), rep(state.name[1:6], each=100))
#'test.imp.fun <- imimp(yhat_b, calibYhat_b, calibY_b, cdeType = "bySample", 
#'                      x.mar = 50, y.mar=100, sampleID=sample_b)
#'
#' @export

imimp <-
function(yhat, calibYhat, calibY, cdeType = c("simple", "local", "bySample"), sampleID=NULL, n.mi=10, x.mar=50, y.mar=1000){
  
  if(cdeType=="simple"){
    res.mat <- matrix(NA, nrow = length(yhat), ncol=n.mi)
    for(i in 1:n.mi){
      res.mat[,i]<- sample(calibYhat, length(yhat), replace = T) + yhat
    }
  }
  
  if(cdeType=="local"){
    res.mat <- matrix(NA, nrow = length(yhat), ncol=n.mi)
    for(i in 1:n.mi){
      res.mat[,i]<- CDE(calibY, calibYhat, yhat)
    }
  }
  
  if(cdeType=="bySample"){
    if(is.null(sampleID)){sampleID=rep("1", length(yhat))}
    res.mat <- matrix(NA, nrow = length(yhat), ncol=n.mi)
    test.slides <- unique(sampleID[[2]])
    train.slides <- unique(sampleID[[1]])
    for(i in 1:n.mi){
      for(test.i in test.slides){
        index.test.i <- which(sampleID[[2]] == test.i)
        train.slide.i <- sample(train.slides, 1)
        index.cali.i <- which(sampleID[[1]] == train.slide.i)
        res.mat[index.test.i, i]<- CDE(calibY[index.cali.i], calibYhat[index.cali.i], 
                                             yhat[index.test.i], xmar = x.mar, ymar=y.mar)
      }
    }
    res.mat <- as.data.frame(res.mat)
    res.mat$slides <- sampleID[[2]]
  }
  
  return(res.mat)
}
