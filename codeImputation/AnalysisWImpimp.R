###########################
########## Set up #########
###########################

library(parallel)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(knitr)
library(DT)
library(qreport)
library(hdrcde)
library(reshape2)


load('/media/disk2/atlas_mxif/data/cell_phenotyped.rds')
# cell = readRDS('/media/disk2/atlas_mxif/data/threshold_quantification_global_20230315.rds')
# cell$Slide_Region <- paste(cell$SlideID, cell$region, sep="_")

# Markers that have no missingness in all datasets + MUC5AC
commons <- c("BCATENIN", "CD11B", "CD20", "CD3D", "CD4", "CD45", "CD68", "CD8",      
             "CDX2", "CGA", "COLLAGEN", "ERBB2", "FOXP3", "HLAA", "LYSOZYME", "MUC2",     
             "MUC5AC", "NAKATPASE", "OLFM4", "PANCK", "PCNA", "PEGFR", "PSTAT3", "SMA",
             "SOX9", "VIMENTIN")

allMarkers = names(sc)[grepl("Median_Cell_",names(sc))]
nzNormedMarkers = paste0('nzNorm_',allMarkers)
sc = do.call(rbind, lapply(split(sc, sc$SlideID), function(df){df[,nzNormedMarkers] = log10(1+sweep(df[,allMarkers], 2, colMeans(replace(df[,allMarkers], df[,allMarkers]==0, NA), na.rm = TRUE ), FUN = '/' ) ); df }) )
slides <- unique(sc$SlideID)
sc$batch <- ifelse(sc$SlideID %in% slides[17:42], "MAP", "HTA")
sc_HTAfull <- sc[sc$batch=="HTA", ]
HTAslides <- slides[1:16]

################################
########## fold splits #########
################################

## 6 training slides, 6 testing slides # all from HTA
AD_slides <- HTAslides[table(sc_HTAfull$SlideID, sc_HTAfull$broadTissue)[,1]>0]
SSL_slides <- HTAslides[table(sc_HTAfull$SlideID, sc_HTAfull$broadTissue)[,2]>0]
## training 
set.seed(1)
sampleslide <- c(sample(AD_slides,6), sample(SSL_slides,6))
trainSlides <- sampleslide[c(1:3,7:9)]

## testing
testSlides <-  sampleslide[c(4:6,10:12)]

################################
########## fit model ###########
################################

## use lasso 

# Function for create data
marker="MUC5AC"
createFoldData <- function(dat=sc_HTAfull, slideTrain=slideTrain, marker="MUC5AC"){
  commons <- c("BCATENIN", "CD11B", "CD20", "CD3D", "CD4", "CD45", "CD68", "CD8",      
               "CDX2", "CGA", "COLLAGEN", "ERBB2", "FOXP3", "HLAA", "LYSOZYME", "MUC2",     
               "MUC5AC", "NAKATPASE", "OLFM4", "PANCK", "PCNA", "PEGFR", "PSTAT3", "SMA",
               "SOX9", "VIMENTIN")
  # normalized expression name in the dataset
  nzmarkername <- paste("nzNorm_Median_Cell",marker, sep="_")
  #define predictor and response variables in testing set. Use the common markers for all 3 batches as the prediction variables
  commons_regress <- paste("nzNorm_Median_Cell",commons, sep="_")
  commons_regress <- commons_regress[commons_regress!=nzmarkername]
  HTA_train1 <- as.matrix(dat[dat$SlideID%in%slideTrain,commons_regress])
  # Need to code tissue type and tumor (binary) as discrete numbers to run XGBoost
  HTA_train2 <- cbind(HTA_train1, apply(dat[dat$SlideID%in%slideTrain,c("broadTissue","Tumor")], 2, function(x)as.numeric(as.factor(x))))
  HTA_train_marker <- dat[dat$SlideID%in%slideTrain, nzmarkername]
  return(list(HTA_train1, HTA_train2, HTA_train_marker, dat_original = dat[dat$SlideID%in%slideTrain, ]))
}

## mean imputation
trainData <- createFoldData(dat=sc_HTAfull, 
                            slideTrain = trainSlides, 
                            marker = marker)
cv_model2 <- glmnet:::cv.glmnet(trainData[[2]], trainData[[3]], alpha = 1, nfolds = 5)
# find optimal lambda value that minimizes test MSE
best_lambda2 <- cv_model2$lambda.min
best_model2 <- glmnet::glmnet(trainData[[2]], trainData[[3]], alpha = 1, lambda = best_lambda2)
testData <- createFoldData(dat=sc_HTAfull, 
                          slideTrain = testSlides, 
                          marker = marker)
Yhat.test <- as.numeric(
  predict(best_model2,newx=testData[[2]],s= best_lambda2))

# slide-specific calibration
calibY <- c()
calibYhat <- c()
calibSlides <- c()
for(i in trainSlides){
  trainslidei <- setdiff(trainSlides,i)
  trainDatai <- createFoldData(dat=sc_HTAfull, 
                              slideTrain = trainslidei, 
                              marker = marker)
  
  # model with tumor info
  cv_model2 <- glmnet:::cv.glmnet(trainDatai[[2]], trainDatai[[3]], alpha = 1, nfolds = 5)
  # find optimal lambda value that minimizes test MSE
  best_lambda2 <- cv_model2$lambda.min
  best_model2 <- glmnet::glmnet(trainDatai[[2]], trainDatai[[3]], alpha = 1, lambda = best_lambda2)
  caliData <- createFoldData(dat=sc_HTAfull, 
                             slideTrain = i, 
                             marker = marker)
  Yhat <- as.numeric(
    predict(best_model2,newx=caliData[[2]],s= best_lambda2))
  
  calibY <- c(calibY, caliData[[3]])
  calibYhat <- c(calibYhat, Yhat)
  calibSlides <- c(calibSlides, rep(i, length(Yhat)))
}

################################
########## Use imimp ###########
################################

imp.simple <- imimp(Yhat.test, calibYhat , calibY, cdeType = "simple", x.mar = 50, y.mar=100)
imp.local<- imimp(Yhat.test, calibYhat, calibY, cdeType = "local")
imp.bysample<-imimp(Yhat.test, calibYhat , calibY, cdeType = "bySample",
                    x.mar = 50, y.mar=100, sampleID=list(calibSlides, sc$SlideID[sc$SlideID %in% testSlides]))
# alternatively; use one sample

imp.list <- list(imp.simple, imp.local, imp.bysample)
##############################
######### results ############
##############################

# plot: imputation quantile v.s. real quantile
slides.vec <- sc$SlideID[sc$SlideID %in% testSlides]
tumor.idx <- which(sc$Tumor[sc$SlideID %in% testSlides]==1)
# library(tidyverse)
# library(tidyr)
# library(tibble)
library(reshape2)
quantile.list <- lapply(imp.list, function(imp.mat){
  imp.df <- data.frame(imp.mat, slides=slides.vec)
  imp.df$tumorType <- sc$broadTissue[which(sc$SlideID %in% testSlides)]
  imp.df <- imp.df[tumor.idx,]
  imp.vec <- melt(imp.df, id.vars = c("slides", "tumorType"), variable.name = "im")
  quantile.df <- as.data.frame(imp.vec %>% 
                              group_by(slides, im, tumorType) %>%  
                              summarise(quantile = scales::percent(c(0.9, 0.95, 0.99)),
                                        value = quantile(value, c(0.9, 0.95, 0.99))))
  return(quantile.df)
})
quantile.list.df <- do.call(rbind, quantile.list)
quantile.list.df$methods <- rep(c("simple", "local", "bySample"), each=(nrow(quantile.list.df)/3))
quantile.imp.summary <- as.data.frame(quantile.list.df %>% group_by(quantile, im, tumorType, methods) %>% summarise(means=mean(value), .groups = "keep")
)
data.subset <- sc[which(sc$SlideID %in% testSlides),c(paste0("nzNorm_Median_Cell_", marker), "SlideID", "broadTissue")]
data.subset <- data.subset[tumor.idx,]
names(data.subset) <- c("value", "slides", "tumorType")
real.quantiles <- data.subset %>% 
  group_by(slides) %>%  
  summarise(tumorType=first(tumorType),
    quantile = scales::percent(c(0.9, 0.95, 0.99)),
            real_value = quantile(value, c(0.9, 0.95, 0.99)))
real.quantiles.summary <- as.data.frame(real.quantiles %>% group_by(quantile, tumorType) %>% summarise(means=mean(real_value), .groups = "keep")
)
ll <- list(quantile.list.df, real.quantiles)

save(ll, file = paste0("~/ImageImputation/plots/imimp/", marker,"quantileIMP.RData"))

## plot
library(ggplot2)

ggplot()+
  geom_boxplot(data=quantile.imp.summary , mapping = aes(x=means, y=methods, group=methods, fill=tumorType))+
  geom_vline(real.quantiles.summary, mapping = aes(xintercept=means, group=tumorType), linetype=5, linewidth=1, alpha=0.5)+
  facet_grid(tumorType~quantile, scales = "free_x")+ggtitle("Mean quantiles across imputations", marker)+theme_bw()
ggsave(paste0("~/ImageImputation/plots/imimp/all_", marker,".png"), width=9.85, height=7.79, dpi=300)

######################################
############Sanity Check##############
######################################

# check 1: for each testing slide, split 50-50
########## use the first half as calibration 
########## and the second half as testing
########## separate slides

## mean imputation
trainData <- createFoldData(dat=sc_HTAfull, 
                            slideTrain = trainSlides, 
                            marker = marker)
cv_model2 <- glmnet:::cv.glmnet(trainData[[2]], trainData[[3]], alpha = 1, nfolds = 5)
# find optimal lambda value that minimizes test MSE
best_lambda2 <- cv_model2$lambda.min
best_model2 <- glmnet::glmnet(trainData[[2]], trainData[[3]], alpha = 1, lambda = best_lambda2)

# results
test.check <- list()

for(i in testSlides){
  testData <- createFoldData(dat=sc_HTAfull, 
                             slideTrain = i, 
                             marker = marker)
  Yhat.test.check <- as.numeric(
    predict(best_model2,newx=testData[[2]],s= best_lambda2))
  idx.calib <- sample(1:length(testData[[3]]), round(length(testData[[3]])/2),replace = FALSE)
  
  
  imp.simple <- imimp(Yhat.test.check[-idx.calib], Yhat.test.check[idx.calib], testData[[3]][idx.calib], cdeType = "simple", x.mar = 50, y.mar=100, test.check = TRUE)
  imp.simple.hist <- imp.simple[[1]]
  imp.simple.mat <- imp.simple[[2]]
  colnames(imp.simple.mat) <- paste0("simple", 1:10)
  imp.simple <- data.frame(imp.simple.mat)
  imp.simple.df <- melt(imp.simple, variable.name = "impGroup")
  imp.local<- imimp(Yhat.test.check[-idx.calib], Yhat.test.check[idx.calib], testData[[3]][idx.calib], cdeType = "local", test.check = TRUE)
  imp.local.hist <- imp.local[[1]]
  imp.local.mat <- imp.local[[2]]
  colnames(imp.local.mat) <- paste0("local", 1:10)
  imp.local <- data.frame(imp.local.mat)
  imp.local.df <- melt(imp.local, variable.name = "impGroup")
  trueY.df <- data.frame(value=testData[[3]][-idx.calib], impGroup="trueY")
  mean.imp.df <- data.frame(value=Yhat.test.check[-idx.calib], impGroup="mean_imp")
  train.error <- data.frame(yhat=Yhat.test.check[idx.calib], error=testData[[3]][idx.calib]-Yhat.test.check[idx.calib])
  test.check.df <- data.frame(imp.simple.hist, imp.local.hist)
  output.df <- rbind(imp.simple.df, imp.local.df, trueY.df, mean.imp.df)
  
  # positive percentage and count
  posCount <- paste0("Number of postive cells: ", sum(testData[[3]]>0), "fraction: ", mean(testData[[3]]>0))
  
  test.check[[i]] <- list(test.check.df, output.df, train.error, posCount)
  }

save(test.check, file=paste0("~/ImageImputation/test.check.", marker, ".RData"))

# Check for error distribution: why local is not as good as simple?
i <- testSlides[6]
testData <- createFoldData(dat=sc_HTAfull, 
                           slideTrain = i, 
                           marker = marker)
Yhat.test.check <- as.numeric(
  predict(best_model2,newx=testData[[2]],s= best_lambda2))
idx.calib <- sample(1:length(testData[[3]]), round(length(testData[[3]])/2),replace = FALSE)

imp.simple <- imimp(Yhat.test.check[-idx.calib], Yhat.test.check[idx.calib], testData[[3]][idx.calib], cdeType = "simple", x.mar = 50, y.mar=100)
colnames(imp.simple) <- paste0("simple", 1:10)
imp.simple <- data.frame(imp.simple)
imp.simple.df <- melt(imp.simple, variable.name = "impGroup")

imp.local<- imimp(Yhat.test.check[-idx.calib], Yhat.test.check[idx.calib], testData[[3]][idx.calib], cdeType = "local", plot=TRUE)
colnames(imp.local) <- paste0("local", 1:10)
imp.local <- data.frame(imp.local)
imp.local.df <- melt(imp.local, variable.name = "impGroup")

trueY.df <- data.frame(value=testData[[3]][-idx.calib], impGroup="trueY")
mean.imp.df <- data.frame(value=Yhat.test.check[-idx.calib], impGroup="mean_imp")

load(paste0("~/ImageImputation/figure/list_check6", marker, ".RData"))
error_calib <- list_check[[1]]
sample01 <- list_check[[2]]
# "HTA11_8099_0000_02_01"
true.error <- data.frame(x=testData[[3]][-idx.calib] - Yhat.test.check[-idx.calib])
ggplot() +
  geom_histogram(error_calib, mapping=aes(x=error_calib,y=after_stat(density)),alpha=0.2, position="identity",bins=1000, fill="red")+
  geom_histogram(sample01, mapping=aes(x=x,y=after_stat(density)),alpha=0.2, position="identity",bins=1000,fill="blue")+
  geom_freqpoly(true.error, mapping=aes(x=x, y=after_stat(density)), bins=1000)+theme_bw()

# TRUE error distribution
fit_cde(Yhat.test.check[-idx.calib], true.error$x, plot = TRUE)
fit_cde(Yhat.test.check[idx.calib], testData[[3]][idx.calib] - Yhat.test.check[idx.calib], plot = TRUE)

# for(i in testSlides){
#    data.i <- test.check[[i]]
#    ggplot() +
#      stat_ecdf(data=data.i, aes(x=value, group=impGroup, color=impGroup),geom = "step", alpha=0.5)+
#      stat_ecdf(data=subset(data.i, impGroup="trueY") , aes(x=value), geom = "step")+theme_bw()
#  ggsave(paste0("~/ImageImputation/plots/imimp/sancheck_", marker,"_",i,".png"), width=8.85, height=7.79, dpi=300)
# }

##############################
######## Rubin's Rule ########
##############################
marker="MUC5AC"
load(paste0("~/ImageImputation/plots/imimp/", marker,"quantileIMP.RData")) #ll
impQuantile <- ll[[1]]
trueQuantile <- ll[[2]]
# within imputation: formula 9.2
within.var <- as.data.frame(impQuantile %>% group_by(methods, quantile, im)%>%dplyr::summarize(SEi = (var(value)*5/6), .groups = "keep"))
within.var.pool <- as.data.frame(within.var%>%group_by(methods, quantile)%>%dplyr::summarize(VARW = mean(SEi), .groups = "keep"))
# between imputation: formula 9.3
between.var <- as.data.frame(impQuantile %>% group_by(methods, quantile, slides)%>%dplyr::summarize(thetai = (mean(value)), .groups = "keep"))
between.var.pool <- as.data.frame(between.var%>%group_by(methods, quantile)%>%dplyr::summarize(VARB = var(thetai), means=mean(thetai), .groups = "keep"))
# pool variance: formula 9.4
all.vars <- join(within.var.pool, between.var.pool, by=c("methods", "quantile"))
all.vars$VARB <- all.vars$VARB*(11/10)
all.vars$vars <- all.vars$VARW+all.vars$VARB
# calculate df
source("~/ImageImputation/dfCalc.R")
all.vars.df <- dfCalc(all.vars$VARB, all.vars$vars, 10 ,6, 1)
all.vars$tval95 <- qt(0.975, all.vars.df)

# real quantile mean and variance
true.plot <- trueQuantile %>% group_by(quantile)%>%dplyr::summarise(means=mean(real_value), vars=var(real_value), methods="true") %>% as.data.frame()
true.plot$tval95 <- qt(0.975, 5)

df.plot <- rbind(true.plot, all.vars[,c("methods", "quantile", "means", "vars", "tval95")])
df.plot$sd <- sqrt(df.plot$vars)

ggplot() +
  geom_errorbar(data=df.plot,aes(xmin=means-tval95*sd, xmax=means+tval95*sd, y=methods), width=.2)+
  geom_point(data=df.plot,aes(x=means, y=methods))+facet_wrap(~quantile, scales = "free_x") +
  ggtitle(marker)+theme_bw()
ggsave(p, file=paste0("~/ImageImputation/figure/RubinPool", marker, ".png"), width=9.85, height=7.79, dpi=300)
