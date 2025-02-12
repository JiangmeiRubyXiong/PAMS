library(RNifti)
library(neurobase)
library(hdrcde)
library(readr)
library(dplyr)
library(parallel)
library(EveTemplate)
library(papayaWidget)
source("~/ImImp/R/CDE.R")
source("~/ImImp/R/imimp.R")
source("~/ImageImputation/dfCalc.R")
source("~/ImageImputation/Functions_imputation.R")
# 4 steps
# 1. Bootstrap: 500 simulations, sample size 200 (out of 89)
# 2. Split: calibration and testing (half and half)
# 3. Generate calibration error from calibration dataset
# 4. Add error to dataset

# File locations:
# Training data: /media/disk2/beijing_dti/Neuroimaging/train_T1_FA/trainingSlices
# Testing data: /media/disk2/beijing_dti/processed/testingData
# Synthesized data: /media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain
# Resliced synthesized data (2mm): 

## In the scenario of application, this equals to 
## having a dataset with about 1/4 individual missing FA but T1 is available
## train the model with 1/2 of the dataset
## generate epirical error with 1/4 of the rest that are not missing

# Generate calibration error

testt_subject_data <- read.csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")
# get errors
FA_diff <- FA_synth_list <- FA_truth_list <- list() # errors for all calibration data
FA_atlas=readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_WMmask.nii.gz")
FA_atlas_dim=readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_WMmask.nii.gz")
data.dir <- "/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain_2mm"


# create a mask
FA_mask <- readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_WMmask.nii.gz")
# the mask is created with ALL SLIDES
for(i in 1:180){
  FA_truth=readNifti(paste0("/media/disk2/beijing_dti/enhanced/unpack/",sprintf("%03d", i),"/scalars_standard/DTI_FA_resliced.nii.gz"))
  print(range(FA_truth))
  # remove background
  FA_mask[FA_truth<0.2] <- 0
}
writeNifti(FA_mask, "/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask2.nii.gz")
FA_mask <- readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask2.nii.gz")
for(i in 0:88){
  id <- i
  id.old <- testt_subject_data$original_row[i+1]
  FA_synth=readNifti(paste0(data.dir, "/subject_", id, "_reassembled_brain.nii.gz"))
  FA_synth[which(is.na(FA_synth))] <- 0
  FA_truth=readNifti(paste0("/media/disk2/beijing_dti/enhanced/unpack/",sprintf("%03d", id.old),"/scalars_standard/DTI_FA_resliced.nii.gz"))

  FA_synth <-  FA_synth[-which(FA_mask==0)];FA_truth <-  FA_truth[-which(FA_mask==0)]
  # print(range(FA_truth))
  FA_synth_list[[i+1]] <- FA_synth
  FA_truth_list[[i+1]] <- FA_truth
}
saveRDS(FA_synth_list, file="/media/disk2/beijing_dti/Neuroimaging/FA_synth_list.rds")
saveRDS(FA_truth_list, file="/media/disk2/beijing_dti/Neuroimaging/FA_truth_list.rds")
# method_param=2; xden_param="uniform"; xmar=50; ymar=100
# ab <- cde.bandwidths(
#   x = FA_synth, y = FA_synth-FA_truth, method=method_param, xden=xden_param
# )
# a <- ab$a
# b <- ab$b
# # use hdr cde
# fit_cde <- hdrcde::cde(
#   x = FA_synth, y = FA_synth-FA_truth,
#   nxmargin = xmar,
#   nymargin = ymar,
#   a=a,
#   b=b
# )
# plot(fit_cde)


# one simulation: 25min for one simulation
btsp_simulation <- function(seed.i, data.dir="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain_2mm",
                            phenotype=FALSE, n.imp=10){
  set.seed(seed.i)
  DIR="/media/disk2/beijing_dti/Neuroimaging/imp_res"
  # bootstrap 200 out of 89
  # read data
  testt_subject_data <- read_csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")
  btsp.idx <- sample(1:89, 200, replace=TRUE)
  # split 100 as calibration, 100 as testing
  btsp.idx.calib <- btsp.idx[1:100]
  
  btsp.idx.testing <- btsp.idx[101:200]
  
  
  #######################
  ###### use imimp ######
  #######################
  FA_synth_list <- readRDS("/media/disk2/beijing_dti/Neuroimaging/FA_synth_list.rds")
  FA_truth_list <- readRDS("/media/disk2/beijing_dti/Neuroimaging/FA_truth_list.rds")
  
  Yhat.test <- do.call(c, FA_synth_list[btsp.idx.testing])
  calibYhat <- do.call(c, FA_synth_list[btsp.idx.calib])
  calibY <- do.call(c, FA_truth_list[btsp.idx.calib])
  
  imp.simple <- imimp(Yhat.test, calibYhat , calibY, cdeType = "simple", x.mar = 50, y.mar=100, n.mi=n.imp)
  
  imp.local<- imimp(Yhat.test, calibYhat, calibY, cdeType = "local", x.mar = 50, y.mar=100, n.mi=n.imp)
  samples.id.calib <- mapply(function(a, x){rep(a, length(x))}, a=1:100, x=FA_synth_list[btsp.idx.calib])
  samples.id.testing <- mapply(function(a, x){rep(a, length(x))}, a=101:200, x=FA_synth_list[btsp.idx.testing])
  if(phenotype){
    calib.group <- testt_subject_data$Sex[btsp.idx.calib]
    calib.F <- btsp.idx.calib[calib.group=="F"]
    calib.M <- btsp.idx.calib[calib.group=="M"]
    testing.group <- testt_subject_data$Sex[btsp.idx.testing,]
    testing.F <- btsp.idx.testing[testing.group=="F"]
    testing.M <- btsp.idx.testing[testing.group=="M"]
    
    Yhat.test.F <- do.call(rbind, FA_synth_list[[btsp.idx.testing.F]])
    calibYhat.F <- do.call(rbind, FA_synth_list[[btsp.idx.calib.F]])
    calibY.F <- do.call(rbind, FA_truth_list[[btsp.idx.calib.F]])
    
    Yhat.test.M <- do.call(rbind, FA_synth_list[[btsp.idx.testing.M]])
    calibYhat.M <- do.call(rbind, FA_synth_list[[btsp.idx.calib.M]])
    calibY.M <- do.call(rbind, FA_truth_list[[btsp.idx.calib.M]])
    
    testSlidesAll <- dataset$File.Name[dataset$File.Name %in% testSlides]
    
    imp.bySample <- as.list(seq (1,100,1))
    imp.bySample[[which(testing.group=="F")]] <- imimp(Yhat.test.F, calibYhat.F, calibY.F, cdeType = "bySample",
                                           x.mar = 50, y.mar=100,
                                           sampleID=list(unlist(samples.id.calib[[calib.group=="F"]]),
                                                         unlist(samples.id.testing[[testing.group=="F"]])),
                                           n.mi=n.imp)
    imp.bySample[[which(testing.group=="M")]] <- imimp(Yhat.test.M, calibYhat.M, calibY.M, cdeType = "bySample",
                                                       x.mar = 50, y.mar=100,
                                                       sampleID=list(unlist(samples.id.calib[[calib.group=="M"]]),
                                                                     unlist(samples.id.testing[[testing.group=="M"]])),
                                                       n.mi=n.imp)
  }else{
    imp.bySample<-imimp(Yhat.test, calibYhat , calibY, cdeType = "bySample",
                        x.mar = 50, y.mar=100, sampleID=list(unlist(samples.id.calib),unlist(samples.id.testing)))
    }
  imp.list <- list(Yhat.test, imp.simple, imp.local, imp.bySample)
  saveRDS(imp.list, file=paste0(DIR, "/imputation", seed.i, ".rds"))
}

seq.list <- as.list(seq (1,200,1)) 
mclapply(seq.list, btsp_simulation, mc.cores=10, mc.preschedule = FALSE)


norm.datavec <- function(datavec){
  datavec[datavec>1] <- 1
  datavec[datavec<0] <- 0
  data.out <- datavec
  return(data.out)
}

# reshape the result
FA_atlas_zeros= readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask2.nii.gz")
imp_to_nifti <- function(data.vec, path, mask=FA_atlas_zeros){
  data.vec <- norm.datavec(data.vec)
  mask[mask!=0] <- data.vec
  # Save as a NIfTI file
  writeNifti(mask, file = paste0(path, ".nii.gz"))
}
#
# split result to individual vectors
res_to_vec <- function(i){
  reslist.path <- paste0("/media/disk2/beijing_dti/Neuroimaging/imp_res/imputation", i,".rds")
  res.list <- readRDS(reslist.path)
  res.folder <- paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_brains/simulation", i,"/")
  dir.create(res.folder)
  
  patient.split <- rep(1:100, each = 27408)
  yhat <- split(res.list[[1]], patient.split)
  res.folder1 <- paste0(res.folder,"meanImp/")
  dir.create(res.folder1)
  for(yhat.i in 1:length(yhat)){
    pathi=paste0(res.folder1, "meanImp",yhat.i)
    imp_to_nifti(yhat[[yhat.i]], pathi)
  }
  # create multiple imputation folders
  multi.folders <- c(paste0(res.folder,"simpleImp/"),
                     paste0(res.folder,"localImp/"),
                     paste0(res.folder,"hierImp/"))
  for(foldi in multi.folders){
    dir.create(foldi)
    for(imp.i in 1:10){
      dir.create(paste0(foldi, "Imputation", imp.i, '/'))
    }
  }
  
  
  for(listi in 2:4){
    multImp=res.list[[listi]]
    for(imp.i in 1:10){
      multImp.i <- split(multImp[,imp.i], patient.split)
      for(ii in 1:length(multImp.i)){
        pathi=paste0(multi.folders[(listi-1)], "Imputation", imp.i, "/meanImp", ii)
        imp_to_nifti(multImp.i[[ii]], pathi)
      }
    }
  }
}

mclapply(seq.list, res_to_vec, mc.cores=10, mc.preschedule = FALSE)


# create list of files for atlas
files_brain <- list.files("/media/disk2/beijing_dti/Neuroimaging/imputation_brains/", recursive = TRUE, full.names = FALSE)
files_brain0 <- paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_brains/",files_brain)
files_brain_out <- gsub("\\.nii\\.gz$", "", files_brain)
files_brain_out <- paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_csv/",files_brain_out,".csv")
write.table(files_brain0, file = "/media/disk2/beijing_dti/Neuroimaging/listBrainFiles.txt", row.names = FALSE, col.names = FALSE)
write.table(files_brain_out, file = "/media/disk2/beijing_dti/Neuroimaging/listBraincsv.txt", row.names = FALSE, col.names = FALSE)

# Downstream Analysis
# Combine all the csvs
# row_LR_pair <- data.frame(left=c( 8,  2, 39, 10, 5, 15, 19,
#                                   11,  4, 16, 20,  7,   6, 27,  3, 14, 1, 18, 59),
#                           right=c(70, 13, 99, 72,  67, 76, 80,
#                                   73, 66, 77, 81, 69, 68, 55, 65, 75, 64, 79, 60), 
#                           names=c("ANGULAR", "CAUDATE_NUCLEUS", 'Cingulum_cingulate_gyrus', 'CUNEUS',
#                                   'INFERIOR_FRONTAL', 'INFERIOR_OCCIPITAL',
#                                   'INFERIOR_TEMPORA', 'LINGUAL', 'MIDDLE_FRONTAL', 'MIDDLE_OCCIPITAL',
#                                   'MIDDLE_TEMPORA', 'POSTCENTRAL', 'PRECENTRAL',  'PUTAMEN',
#                                   'SUPERIOR_FRONTAL', 'SUPERIOR_OCCIPITAL', 'SUPERIOR_PARIETAL',
#                                   'SUPERIOR_TEMPORAL', 'THALAMUS'))
# row_LR_single <- data.frame(singles = c(12, 52, 40, 51, 53, 71),
#                            names= c('Fusiform_right',"Body_of_corpus_callosum", "Cingulum_hippocampus_left",
#                                     "Genu_of_corpus_callosum_left", "Splenium_of_corpus_callosum_left", 'Precuneus_left'))

# keep only white matter


row_LR_pair <- data.frame(left=c( 8, 39, 10,
                                  5, 15, 19,
                                  11,  4, 16, 20,
                                  7,   6, 
                                  3, 14, 1, 18),
                          right=c(70,  99, 72,
                                  67, 76, 80,
                                  73, 66, 77, 81, 
                                  69, 68, 
                                  65, 75, 64, 79), 
                          names=c("ANGULAR",  'Cingulum_cingulate_gyrus', 'CUNEUS',
                                  'INFERIOR_FRONTAL', 'INFERIOR_OCCIPITAL', 'INFERIOR_TEMPORA',
                                  'LINGUAL', 'MIDDLE_FRONTAL', 'MIDDLE_OCCIPITAL', 'MIDDLE_TEMPORA',
                                  'POSTCENTRAL', 'PRECENTRAL',
                                  'SUPERIOR_FRONTAL', 'SUPERIOR_OCCIPITAL', 'SUPERIOR_PARIETAL', 'SUPERIOR_TEMPORAL'))
row_LR_single <- data.frame(singles = c(74, 9,
                                        52, 40, 51, 53),
                            names= c('Fusiform_right', 'Precuneus_left',
                                     "Body_of_corpus_callosum", "Cingulum_hippocampus_left",
                                     "Genu_of_corpus_callosum_left", "Splenium_of_corpus_callosum_left"))

library(stringr)

# Function to convert to proper case
format_camel_case <- function(input_string) {
  parts <- strsplit(input_string, "_")[[1]]  # Split the string by underscores
  parts <- tolower(parts)                   # Convert all parts to lowercase
  parts[1] <- paste0(toupper(substr(parts[1], 1, 1)), substr(parts[1], 2, nchar(parts[1]))) # Capitalize first word
  paste(parts, collapse = "_")              # Combine parts back with underscores
}

row_LR_pair$names <- unlist(lapply(row_LR_pair$names, format_camel_case))

row_wm <- list(row_LR_pair, row_LR_single)

csv_filter <- function(input_path, rows=row_wm, input_base ="/media/disk2/beijing_dti/Neuroimaging/imputation_csv",
                       output_base = "/media/disk2/beijing_dti/Neuroimaging/imputation_rds"){
  
  df_csv <- read.delim(input_path, sep="", header = FALSE, stringsAsFactors = FALSE)
  df_csv_mean <- df_csv[-(1:2), 1:2]
  names(df_csv_mean) <- c("LabelID", "Mean")
  
  # combine the region names
  row_left <- rows[[1]]$left;row_right <- rows[[1]]$right; 
  df_csv_mean_LR <- (as.numeric(unlist(subset(df_csv_mean, LabelID %in% row_left, Mean)))+
                       as.numeric(unlist(subset(df_csv_mean, LabelID %in% row_right, Mean))))/2
  df_csv_mean_S <- as.numeric(unlist(subset(df_csv_mean, LabelID %in% rows[[2]]$singles, Mean)))
  df_csv_out <- data.frame(names=c(rows[[1]]$names, rows[[2]]$names),
                           value=c(df_csv_mean_LR, df_csv_mean_S))
  
  # Derive the output file path by replacing the base directory and extension
  relative_path <- sub(input_base, "", input_path) # Get relative path from the base directory
  output_path <- paste0(output_base, dirname(relative_path), "/",
                           sub("\\.csv$", ".rds", basename(relative_path)))
  
  # Create the output directory if it doesn't exist
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  saveRDS(df_csv_out, output_path)
}
files_brain_out <- read.table("/media/disk2/beijing_dti/Neuroimaging/listBraincsv.txt")[[1]]
# Apply the function to all CSV files
mclapply(files_brain_out, csv_filter, mc.cores=30, mc.preschedule = FALSE)


##### sex variable in all simulations
testt_subject_data <- read_csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")
btsp_simulation_sex <- function(seed.i, data.dir="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain_2mm",
                            phenotype=FALSE, n.imp=10, sub_data = testt_subject_data){
  set.seed(seed.i)
  DIR="/media/disk2/beijing_dti/Neuroimaging/imp_res"
  # bootstrap 200 out of 89
  # read data

  btsp.idx <- sample(1:89, 200, replace=TRUE)
  # split 100 as calibration, 100 as testing
  btsp.idx.calib <- btsp.idx[1:100]
  
  btsp.idx.testing <- btsp.idx[101:200]
  return(sub_data$Sex[btsp.idx.testing])
  }

sex.var <- lapply(1:200, btsp_simulation_sex)

# Combine the data results
singleImpData <- function(path.data){
  files_imp <- list.files(path = path.data, full.names = TRUE)
  list.df <- lapply(files_imp, readRDS)
  names_brain <- list.df[[1]][[1]]
  list.df <- data.frame(t(do.call(cbind, lapply(list.df, "[[", 2))))
  names(list.df) <- names_brain
  return(list.df)
}
# GLM
regres <- function(dataset){
  modeli <- glm(sex~., family = "binomial",dataset)
  sum_modeli <- summary(modeli)
  return(sum_modeli$coefficients[,1:2])
}

# logistic regression Function

# this function is applied to one simulation one type of MULTIPLE imputation
gee_results_neuro <- function(sim.idx, type="hier",sex.vars=sex.var){

  sex.var.i <- sex.var[[sim.idx]]
  DIR="/media/disk2/beijing_dti/Neuroimaging/imputation_rds/simulation1/hierImp/Imputation1/"
  DIR = paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_rds/simulation", sim.idx, "/", type, "Imp")
  files.sim <- list.dirs(DIR, full.names = TRUE)[-1]
  imp.dataset <- lapply(files.sim, singleImpData)
  imp.dataset <- lapply(imp.dataset, cbind, sex=sex.var.i)
  # combine DIR# combine all the individuals in the 
  if(!type %in% c("truth, mean")){
    regres.out <- lapply(imp.dataset, regres)
    ### use Rubin's rule to combine
    df.plot.geediff$means[idx] <- mean(m0.param$mean)
    VARB <- var(m0.param$mean)
    # vars <- sqrt(var(m0.param$mean)*(n.imp+1)/n.imp+mean((m0.param$SE)^2))
    # df.plot.geediff$sd[idx] <- vars
    # df.plot.geediff$stats[idx] <- qt(0.975, dfCalc(VARB, vars, 10 , 200, 2))
  }
   
  
  
  
  
  for(i in quantile3){
    ### true data for simulation subset
    m0 <- gee(real_value ~ groups, id=patients, family='gaussian',
              corstr = "exchangeable", data=subset(trueData, quantile==i))
    idx <- df.plot.geediff$quantile==i&df.plot.geediff$methods=="true"
    df.plot.geediff$means[idx] <- m0$coefficients[2]
    df.plot.geediff$sd[idx] <- sqrt(m0$robust.variance[2,2])
    df.plot.geediff$stats[idx] = qnorm(0.975)
    ### mean imputation
    m0 <- gee(real_value ~ groups, id=patients, family='gaussian',
              corstr = "exchangeable", data=subset(meanData, quantile==i))
    idx <- df.plot.geediff$quantile==i&df.plot.geediff$methods=="mean"
    df.plot.geediff$means[idx] <- m0$coefficients[2]
    df.plot.geediff$sd[idx] <- sqrt(m0$robust.variance[2,2])
    df.plot.geediff$stats[idx] = qnorm(0.975)
    
    for(m in c("simple", "local", "bySample")){
      n.imp <- length(levels(quantile.list.df.strat$im))
      m0.param <- data.frame(mean=rep(NA, n.imp), SE=rep(NA, n.imp))
      names(trueQuantiletrain)[which(names(trueQuantiletrain)=="real_value")] <- "value"
      for(ii in 1:n.imp){
        dataImp <- subset(quantile.list.df.strat, (quantile==i&methods==m&im==paste0("X",ii)))
        trueQuantiletrainSub <- subset(trueQuantiletrain, quantile==i)
        dataImp <- rbind(trueQuantiletrainSub, dataImp[, names(trueQuantiletrain)])
        m0 <- gee(value ~ groups, id=patients, family='gaussian', corstr = "exchangeable",
                  data=dataImp)
        m0.param[ii,"mean"] <- m0$coefficients[2]
        m0.param[ii,"SE"] <- summary(m0)$coefficients[2,4]
      }
      idx <- df.plot.geediff$quantile==i&df.plot.geediff$methods==m

    }
  }
  save(df.plot.geediff, file = paste0(DIR, "df.plot.geediff.RData"))
  
  # calculate metrics:
  true.df.plot.geediff <- get(load(tru.dir))
  #Bias: parameter value: the average difference between true and imputed (baseline: test+train; methods: test, imputation methods)
  
  # Width of confidence interval: ratio between (test+train) and (imputation methods)
  df.plot.geediff$cinf.upper <- df.plot.geediff$means + df.plot.geediff$sd * df.plot.geediff$stats
  df.plot.geediff$cinf.lower <- df.plot.geediff$means - df.plot.geediff$sd * df.plot.geediff$stats
  df.plot.geediff$int.length <- df.plot.geediff$sd * df.plot.geediff$stats
  
  df.plot.geediff.subsetTrue <- subset(df.plot.geediff, methods=="true")
  df.plot.geediff <- subset(df.plot.geediff, methods!="true")
  bias <- df.plot.geediff[, "means"] - true.df.plot.geediff[, "means"]
  
  ratio <- df.plot.geediff[, "int.length"] / true.df.plot.geediff[, "int.length"]
  #Coverage: parameter value: how often on average, each method covers the true parameter (baseline: test+train; methods: test, imputation methods)
  coverage <- (df.plot.geediff[, "cinf.lower"] <= true.df.plot.geediff[, "means"]) &
    (true.df.plot.geediff[, "means"] <= df.plot.geediff[, "cinf.upper"])
  metric.df <- data.frame(quantile = df.plot.geediff[, "quantile"],
                          methods = df.plot.geediff[, "methods"], means=df.plot.geediff[, "means"],
                          btspmean = df.plot.geediff.subsetTrue[, "means"],
                          bias=bias, ratio=ratio, coverage=coverage, 
                          coverage.check=df.plot.geediff$int.length>abs(bias),
                          sd=df.plot.geediff$sd,
                          stats=df.plot.geediff$stats)
  
  save(metric.df, file = paste0(DIR, "metric.df.RData"))
  
  #sanity check: how does the subset truth CI covers the true parameter
  coverage.truth <- (df.plot.geediff.subsetTrue[, "cinf.lower"] <= true.df.plot.geediff[, "means"]) &
    (true.df.plot.geediff[, "means"] <= df.plot.geediff.subsetTrue[, "cinf.upper"])
  
  save(coverage.truth, file = paste0(DIR, "coverageTruth.RData"))
  
  # return gee results
}

#### DO the same thing for the TRUTH



####### Combine the calibration data
