n.rep <- 500

#Bias: parameter value: the average difference between true and imputed (baseline: test+train; methods: test, imputation methods)
#Width of confidence interval: ratio between (test+train) and (imputation methods)
#Coverage: parameter value: how often on average, each method covers the true parameter (baseline: test+train; methods: test, imputation methods)

library(foreach)
library(doParallel)

library(GammaGateR)
library(data.table)
library(gridExtra)
library(ggpubr)
library(car)
library(RESI)
library(emmeans)
library(parallel)
library(hrbrthemes)
library(ggplot2)
library(fossil)
library(dplyr)
library(nlme)
library(glmnet)
library(modelr)
library(xgboost)
library(dplyr)
library(caret)
library(hdrcde)
library(reshape2)
library(purrr)
library(gee)
source("~/ImImp/R/CDE.R")
source("~/ImImp/R/imimp.R")
source("~/ImageImputation/dfCalc.R")
source("~/ImageImputation/Functions_imputation.R")
load("~/ImageImputation/nzNormedMarkers.RData")
# dataSchurch<- read.csv(file='/media/disk2/atlas_mxif/data/CRC_clusters_neighborhoods_markers.csv', header = TRUE)
# originalNames <- names(dataSchurch)
# markers <- originalNames[grep("\\.\\.\\.", originalNames)] |>
#   strsplit("\\.\\.\\.") |> lapply("[[",1) |> unlist()
# marker_datanames <- originalNames[grep("\\.\\.\\.", originalNames)]
# 
# allMarkers = markers
# nzNormedMarkers = paste0('nzNorm_',allMarkers)
# dataSchurch = do.call(rbind, lapply(split(dataSchurch, dataSchurch$TMA_AB), function(df){df[,nzNormedMarkers] = log10(1+sweep(df[,marker_datanames], 2, colMeans(replace(df[,marker_datanames], df[,marker_datanames]==0, NA), na.rm = TRUE ), FUN = '/' ) ); df }) )
# 
# dataSchurchLite <- dataSchurch[,c("File.Name", "Region", "groups", "patients", nzNormedMarkers)]
# write.csv(dataSchurchLite, file="~/ImageImputation/dataSchurchLite.csv", row.names = FALSE)
# save(nzNormedMarkers, file="~/ImageImputation/nzNormedMarkers.RData")

dataSchurch <- read.csv("~/ImageImputation/dataSchurchLite.csv")
patients <- unique(dataSchurch$patients)
patientGroup <- dataSchurch %>% group_by(patients, groups)%>%reframe()
patientCLR <- unlist(subset(patientGroup, groups==1, patients))
patientDII <- setdiff(patients, patientCLR)
dataSchurch.patientsplit <- split(dataSchurch, dataSchurch$patients)
#setup parallel backend to use many processors
cl <- 20
# registerDoParallel(cl)
sample.size=200
# foreach(rep.i = 1:500) %dopar% {
for (marker in c("CD20","CD38","CD138","CD4", "CD8", "CD68")){
for(rep.i in 1:500){
  print(c(marker, rep.i))
  set.seed(rep.i)
  # Bootstrap Dataset
  patients.sampled <- c(sample(patientCLR, (sample.size/2), replace=TRUE),
                        sample(patientDII, (sample.size/2), replace=TRUE))
  resampled_data <- dataSchurch.patientsplit[patients.sampled]
  names(resampled_data) <- as.character(1:sample.size)
  resampled_data <- do.call(rbind, resampled_data)
  new.id <- unlist(lapply(strsplit(rownames(resampled_data), "\\."), "[[", 1))
  subData <- resampled_data
  subData$patients <- new.id
  subData$File.Name <- paste0(subData$patients, subData$File.Name)
  PGtable <- table(subData$patients, subData$groups)
  CLR_patients <- rownames(PGtable)[PGtable[,1]>0]
  DII_patients <- rownames(PGtable)[PGtable[,2]>0]
  CLR_total <- length(CLR_patients)
  DII_total <- length(DII_patients)
  CLR_training <- round(CLR_total*2/3)
  DII_training <- round(DII_total*2/3)
  ## training
  trainPatients <- c(CLR_patients[sample(1:CLR_total, CLR_training)],
                     DII_patients[sample(1:DII_total, DII_training)])
  ## testing
  testPatients <-  setdiff(1:sample.size, trainPatients)
  ## slide - patient correspondence
  SlidePatients <- subData %>% group_by(patients, File.Name) %>% reframe()
  ## training slides
  trainSlide <- SlidePatients$File.Name[which(SlidePatients$patients %in% trainPatients)]
  testSlide <- SlidePatients$File.Name[which(SlidePatients$patients %in% testPatients)]
  # dir.create(Dir)
  # save(patients.sampled, file=paste0(Dir, "samplePatients.RData"))
  # imputation(marker=marker, trainSlides=trainSlide, testSlides=testSlide,
  # dataset=subData, DIR=Dir, phenotype = FALSE)
  Dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_", marker, "_phenotype/rep", rep.i, "/")
  # f.name <- paste0(Dir, "subDatafit", marker, ".RData")
  # imp.list.strat <- get(load(f.name))
  # summary_quantile(imp.result = imp.list.strat,marker=marker, testSlides=testSlide,
  #                  dataset=subData,DIR=Dir)
  f.name <- paste0(Dir, "dataSchurch", marker,"quantileIMP.RData")
  resultList <- get(load(f.name))

  # descriptive_quantiles(resultList, DIR=Dir)
  gee_metric(result.list = resultList, DIR=Dir, dataset=subData, tru.dir=paste0("~/ImageImputation/metricsAllMarker/metrics_",marker, "_phenotype/true.df.plot.geediff.RData"))

}}
 #stop cluster
stopCluster(cl)

marker="CD38"
# for(marker in c("CD38","CD138","CD4", "CD8", "CD68")){
dirs=paste0("~/ImageImputation/metricsAllMarker/metrics_", marker, "_phenotype/") 

# for(rep.i in 1:500){
  
foreach(rep.i = 1:500) %dopar% {
  Dir <- paste0(dirs,"rep", rep.i, "/")
  if(length(list.files(Dir))<5){
    print(rep.i)
    set.seed(rep.i)
    
    # Bootstrap Dataset
    patients.sampled <- c(sample(patientCLR, (sample.size/2), replace=TRUE),
                          sample(patientDII, (sample.size/2), replace=TRUE))
    
    resampled_data <- dataSchurch.patientsplit[patients.sampled]
    names(resampled_data) <- as.character(1:sample.size)
    resampled_data <- do.call(rbind, resampled_data)
    new.id <- unlist(lapply(strsplit(rownames(resampled_data), "\\."), "[[", 1))
    subData <- resampled_data
    subData$patients <- new.id
    subData$File.Name <- paste0(subData$patients, subData$File.Name)
    
    PGtable <- table(subData$patients, subData$groups)
    CLR_patients <- rownames(PGtable)[PGtable[,1]>0]
    DII_patients <- rownames(PGtable)[PGtable[,2]>0]
    
    CLR_total <- length(CLR_patients)
    DII_total <- length(DII_patients)
    CLR_training <- round(CLR_total*2/3)
    DII_training <- round(DII_total*2/3)
    ## training
    trainPatients <- c(CLR_patients[sample(1:CLR_total, CLR_training)], 
                       DII_patients[sample(1:DII_total, DII_training)])
    ## testing
    testPatients <-  setdiff(1:sample.size, trainPatients)
    
    ## slide - patient correspondence
    SlidePatients <- subData %>% group_by(patients, File.Name) %>% reframe()
    ## training slides
    trainSlide <- SlidePatients$File.Name[which(SlidePatients$patients %in% trainPatients)]
    testSlide <- SlidePatients$File.Name[which(SlidePatients$patients %in% testPatients)]
  
  print(rep.i)

    f.name <- paste0(Dir, "subDatafit", marker, ".RData")
    imp.list.strat <- get(load(f.name))
    summary_quantile(imp.result = imp.list.strat,marker=marker, testSlides=testSlide, 
                     dataset=subData,DIR=Dir)
    f.name <- paste0(Dir, "dataSchurch", marker,"quantileIMP.RData")
    resultList <- get(load(f.name))
    # descriptive_quantiles(resultList, DIR=Dir)
    gee_metric(result.list = resultList, DIR=Dir, dataset=subData, tru.dir = paste0(dirs,"true.df.plot.geediff.RData"))
  }
}
# generate true data


marker="CD20"
## put metrics together
metric.list <- metric.list.nophe <- list()
for(rep.i in 1:500){
  Dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker,"_phenotype/rep", rep.i, "/")
  load(paste0(Dir, "metric.df.RData"))
  metric.df$rep <- rep.i
  metric.list[[rep.i]] <- metric.df
  
  Dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker,"_nophe/rep", rep.i, "/")
  load(paste0(Dir, "metric.df.RData"))
  metric.df$rep <- rep.i
  metric.list.nophe[[rep.i]] <- metric.df
}
metric.df.all <- do.call(rbind,metric.list)
rm(list = "metric.list")
metric.df.all$methods <- factor(metric.df.all$methods, levels=rev(c("mean", "simple","local", "bySample")))
metric.df.all$phenotype <- "TRUE"

metric.df.nophe.all <- do.call(rbind,metric.list.nophe)
rm(list = "metric.list.nophe")
metric.df.nophe.all$methods <- factor(metric.df.all$methods, levels=rev(c("mean", "simple","local", "bySample")))
metric.df.nophe.all$phenotype <- "FALSE"

metric.df.all <- rbind(metric.df.all, metric.df.nophe.all)
rm(list = "metric.df.nophe.all")
# plot 1: bias

ggplot(metric.df.all)+
  geom_boxplot(aes(bias, methods, fill=phenotype), alpha=0.7)+facet_wrap(~quantile, scales = "free_x", nrow = 3)+
  geom_vline(xintercept=0)+theme_bw()+ggtitle(paste0("Bias in 500 reps:", marker))

# plot 2: ratio: imputed/true

ggplot(metric.df.all)+
  geom_boxplot(aes(ratio, methods, fill=phenotype), alpha=0.7)+facet_wrap(~quantile, scales = "free_x", nrow=3)+
  geom_vline(xintercept=1)+theme_bw()+ggtitle(paste0("Ratio of imputed/true in 500 reps", marker))

# plot 3: coverage
options(na.rm=TRUE)
metric.coverage <- metric.df.all %>% group_by(quantile, methods, phenotype) %>% summarise(perc=mean(coverage, na.rm=TRUE))
library(tidyr)
View(spread(metric.coverage, methods, perc))
# is the distribution really skewed?
ggplot(dataSchurch, aes(nzNorm_CD20, group=File.Name, after_stat(density))) +
  geom_freqpoly(binwidth = 0.05, alpha=0.1)+
  geom_vline(xintercept = c(0.1180571, 0.4817484 , 1.147745 ))+
  theme_bw() 

# Check: confidence interval
conv.truth.list <- list()
for(rep.i in 1:500){
  Dir <- paste0("~/ImageImputation/metrics/rep", rep.i, "/")
  load(paste0(Dir, "coverageTruth.RData"))
  conv.truth.list[[rep.i]] <- coverage.truth
}
coverages <- do.call(rbind, conv.truth.list)
colMeans(coverages)


# R square
model <- lm(nzNorm_CD38~.,data=dataSchurch[,nzNormedMarkers])
summary(model)
model <- lm(nzNorm_CD38~.,data=dataSchurch[,c(nzNormedMarkers, "groups")])
summary(model)
