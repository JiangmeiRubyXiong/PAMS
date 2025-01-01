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
#' @param plot whether to generate the cde plot. Default is FALSE
#' @param diag.plot whether to generate the diagnostic plot. Default is FALSE
#' @param path The file name and path to where the disgnostic plot is stored. Default to be "~/diag.pdf"
#' @importFrom hdrcde cde cde.bandwidths
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
  function(yhat, calibYhat, calibY, cdeType = c("simple", "local", "bySample"), 
           sampleID=NULL, n.mi=10, x.mar=50, y.mar=1000, plot=FALSE, 
           diag.plot=FALSE, path="~/diag.pdf"){
    
    error_calib <- calibY - calibYhat
    
    if(cdeType=="simple"){
      res.mat <- matrix(NA, nrow = length(yhat), ncol=n.mi)
      for(i in 1:n.mi){
        res.mat[,i]<- sample(error_calib, length(yhat), replace = T) + yhat
      }
    }
    
    if(cdeType=="local"){
      cde.obj <- fit_cde(calibYhat, error_calib, xmar = x.mar, ymar = y.mar, plot=plot)
      res.mat <- sapply(vector("list", length = n.mi), function(x){CDE_sample(cde.obj, yhat) + yhat},
                        simplify = TRUE)
    }
    
    if(cdeType=="bySample"){
      if(is.null(sampleID)){sampleID=rep("1", length(yhat))}
      res.mat <- matrix(NA, nrow = length(yhat), ncol=n.mi)
      test.slides <- unique(sampleID[[2]])
      train.slides <- unique(sampleID[[1]])
      train.calibYhat.list <- split(calibYhat, sampleID[[1]])
      train.error.list <- split(error_calib, sampleID[[1]])
      cde.train.list <- mapply(fit_cde, train.calibYhat.list, train.error.list, MoreArgs = list(xmar = x.mar, ymar = y.mar),SIMPLIFY = FALSE)
      for(i in 1:n.mi){
        for(test.i in test.slides){
          index.test.i <- which(sampleID[[2]] == test.i)
          train.slide.i <- sample(train.slides, 1)
          res.mat[index.test.i, i]<- CDE_sample(cde.train.list[[train.slide.i]], 
                                                yhat[index.test.i])+yhat[index.test.i]
        }
      }
    }
    
    return(res.mat)
  }

imp.diag.scatter <- function(imimp.obj, yhat, ytrue, SlideID=NULL, path="~/diag.pdf"){
  # plot a set of figures: yhat v.s. true error; yhat v.s. error sampled, 
  # AND error sampled for each imputation (?)
  # if SlideID is not specified, let it be a vector of "1"s
  if(is.null(SlideID)){SlideID <- rep("1", length(yhat))}
  slides <- unique(SlideID)
  pdf(file=path)
  for(slide.i in slides){
    par(mfrow=c(3,3))
    slide.idx <- which(SlideID==slide.i)
    x <- yhat[slide.idx]
    for(i in 1:n.mi){
      y <- (res.mat[,i]-yhat)[slide.idx]
      if(slide.idx > 100){
        idx=sample(1:length(x), 100)
        plot(x[idx],y[idx], xlab="Yhat", ylab="Error Sampled", main=paste(slide.i, i, sep=":Imp"))
      } else {
        plot(x,y, xlab="Yhat", ylab="Error Sampled", main=paste(slide.i, i, sep=":Imp"))
      }
    }
  }
  dev.off()
}
