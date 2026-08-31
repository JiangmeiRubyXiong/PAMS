library(RNifti)
library(neurobase)
library(hdrcde)
library(readr)
library(dplyr)
library(parallel)
library(sandwich)
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


## recycle mean imputation result
FA_atlas_zeros= readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask3.nii.gz")
row_wm <- readRDS(file="/media/disk2/beijing_dti/Neuroimaging/brain_ROI.rds")
### remove regions with too little observations (< 100)
small.region <- which(table(c(FA_atlas_zeros))[as.character(row_wm$label)]<100, arr.ind = TRUE)
row_wm_100 <- row_wm[-small.region, ]
format_string <- function(input_str) {
  # Replace "_wm" and "_" with space
  cleaned <- gsub("_wm|_", " ", input_str)
  
  # Split and capitalize (except "of")
  words <- strsplit(cleaned, " ")[[1]]
  words <- sapply(words, function(word) {
    if (tolower(word) == "of") {
      "of"
    } else {
      paste0(toupper(substring(word, 1, 1)), tolower(substring(word, 2)))
    }
  })
  
  # Add \n after "Cingulum" or "of"
  result <- c()
  for (word in words) {
    result <- c(result, word)
    if (word == "Cingulum" || word == "of") {
      result <- c(result, "\n")
    }
  }
  
  # Collapse into one string, remove any extra space before/after \n
  output <- gsub(" \n", "\n", paste(result, collapse = " "))
  output <- gsub("\n ", "\n", output)
  
  return(output)
}

for(i in 1:nrow(row_wm_100)){
  row_wm_100$names[i] <- format_string(row_wm_100$names[i])
}

# saveRDS(row_wm_100, file="/media/disk2/beijing_dti/Neuroimaging/brain_ROI_100.rds")

## Generate region level calibration error for all subjects 
FA_synth_list <- readRDS(file="/media/disk2/beijing_dti/Neuroimaging/FA_synth_list.rds")
FA_truth_list <- readRDS(file="/media/disk2/beijing_dti/Neuroimaging/FA_truth_list.rds")
  ## take region average
norm.datavec <- function(datavec){
  datavec[datavec>1] <- 1
  datavec[datavec<0] <- 0
  data.out <- datavec
  return(data.out)
}
imp_to_csv <- function(data.vec, path, mask=FA_atlas_zeros, row.match = row_wm_100){
  data.vec <- norm.datavec(data.vec)
  brain_mask <- mask
  brain_mask[brain_mask!=0] <- data.vec
  # take regional average
  brain_csv = sapply(split(brain_mask, mask), mean)[-1]
  brain_df <- data.frame(values = brain_csv, label = as.integer(names(brain_csv)))
  row.match <- left_join(row.match, brain_df, by="label")
  return(row.match)
}
FA_synth_region <- lapply(FA_synth_list, imp_to_csv)
FA_truth_region <- lapply(FA_truth_list, imp_to_csv)
saveRDS(FA_synth_region, file="/media/disk2/beijing_dti/Neuroimaging/FA_synth_region.rds")
saveRDS(FA_truth_region, file="/media/disk2/beijing_dti/Neuroimaging/FA_truth_region.rds")

calibration_error <- mapply(function(a,b)return(a$values-b$values), FA_truth_region, FA_synth_region, SIMPLIFY = FALSE)
saveRDS(calibration_error, file="/media/disk2/beijing_dti/Neuroimaging/calibration_error.rds")


FA_synth_region <- readRDS("/media/disk2/beijing_dti/Neuroimaging/FA_synth_region.rds")
FA_truth_region <- readRDS("/media/disk2/beijing_dti/Neuroimaging/FA_truth_region.rds")
calibration_error <- readRDS("/media/disk2/beijing_dti/Neuroimaging/calibration_error.rds")
testt_subject_data <- read.csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")

# store calibration and testing splits for MAR
# testing: 3,7 split for male and female 
set.seed(1)
imputation_seed_generation <- lapply(1:500, 
                                     function(i){
                                       set.seed(i)
                                       test_sampled <- sample(1:nrow(testt_subject_data), 500, 
                                                              replace = TRUE)
                                       sex_sampled <- testt_subject_data[test_sampled,"Sex"]
                                       test_sampled_M <- test_sampled[sex_sampled=="M"]
                                       test_sampled_F <- test_sampled[sex_sampled=="F"]
                                       
                                       n_test_M <- round(length(test_sampled_M)*0.7)
                                       test_M <- test_sampled_M[1:n_test_M]
                                       calib_M <- test_sampled_M[-(1:n_test_M)]
                                       n_test_F <- round(length(test_sampled_F)*0.3)
                                       test_F <- test_sampled_F[1:n_test_F]
                                       calib_F <- test_sampled_F[-(1:n_test_F)]
                                       
                                       btsp.idx.calib <- c(calib_M, calib_F)
                                       btsp.idx.testing <- c(test_M, test_F)
                                       
                                       return(list(btsp.idx.calib, btsp.idx.testing))
                                       })
# saveRDS(imputation_seed_generation, file="/media/disk2/beijing_dti/Neuroimaging/imputation_seed_generation.rds")

imputation_seed_generation <- readRDS("/media/disk2/beijing_dti/Neuroimaging/imputation_seed_generation.rds")

# apply rubin's rule to multiple imputations
region.list <- function(roi, data.list, n.imp=100){
  # reorder the nested list in the multiple imputation lists
  lists <- lapply(data.list, "[[", roi)
  # combine the imputations results
  list.df <- do.call(rbind,lists)
  # rubin's rule calculation
  means <- list.df[,1]
  ses <- list.df[,2]
  VARB <- var(means)
  vars <- sqrt(var(means)*(n.imp+1)/n.imp+mean(ses^2))
  stats <- qt(0.975, dfCalc(VARB, vars, 10 , 500, 2))
  return(data.frame(mean=mean(means), sd=vars, stats=stats))
}
row_wm_100 <- readRDS("/media/disk2/beijing_dti/Neuroimaging/brain_ROI_100.rds")
btsp_imp_df <- function(seed.i, n.imp=100, sub_data = testt_subject_data, 
                          sample.idx = imputation_seed_generation, FA_synth_regions=FA_synth_region,
                        FA_truth_regions = FA_truth_region,calibration_errors= calibration_error,
                          roi.df = row_wm_100){
  set.seed(seed.i)
  # DIR=paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_region/simulation", seed.i)
  # dir.create(DIR)
  # bootstrap 200 out of 89
  # read data
  btsp.idx <- sample.idx[[seed.i]]
  # split 100 as calibration, 100 as testing
  # btsp.idx.calib <- btsp.idx[1:250]
  # btsp.idx.testing <- btsp.idx[251:500]
  
  btsp.idx.calib <- btsp.idx[[1]]
  btsp.idx.testing <- btsp.idx[[2]]
  
  testing.data <- FA_synth_regions[btsp.idx.testing]
  calibration.error.original <- calibration_errors[btsp.idx.calib]
  regions <- testing.data[[1]]$names
  
  # create calibration data
  
  calib.data <- FA_truth_regions[btsp.idx.calib]

  calib.i <- data.frame(do.call(rbind,lapply(calib.data, "[[", "values")))
  names(calib.i) <- regions
  calib.i$Sex <- sub_data$Sex[btsp.idx.calib]
  calib.i$Age <- sub_data$Age[btsp.idx.calib]
  
  # complete case data
  list.res.copmlete <- list()
  for(ROI in regions){
    formulai <- as.formula(paste0(ROI," ~ Sex"))
    m0 <- lm(formulai, data=calib.i)
    estimate <- m0$coefficients[2]
    stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
    list.res.copmlete[[ROI]] <- data.frame(Estimate = estimate, StdError = stdError)
  }

  # Mean imputation
  mean.imp.all <- data.frame(do.call(rbind,lapply(testing.data, "[[", "values")))
  names(mean.imp.all) <- regions
  mean.imp.all$Sex <- sub_data$Sex[btsp.idx.testing]
  mean.imp.all$Age <- sub_data$Age[btsp.idx.testing]
  mean.imp.all <- rbind(mean.imp.all, calib.i)
  
  list.res.mean <- list()
  for(ROI in regions){
    formulai <- as.formula(paste0(ROI," ~ Sex"))
    m0 <- lm(formulai, data=mean.imp.all)
    estimate <- m0$coefficients[2]
    stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
    list.res.mean[[ROI]] <- data.frame(Estimate = estimate, StdError = stdError)
  }
  # list.res <- do.call(rbind, list.res)
  # saveRDS(list.res, paste0(DIR, "/meanImp.rds"))

  # multiple imputation
  
  # 1. NOT stratify by sex

  imputation.list <- list()
  for(i.imp in 1:n.imp){
    testing.data.frame <- roi.df
    ## resample calibration error values
    len.calib.error <- length(calibration.error.original)
    calibration.error <- calibration.error.original[sample(1:len.calib.error, len.calib.error, replace = TRUE)]
    for(ind in 1:length(testing.data)){
      out.testing <- testing.data[[ind]]
      out.testing$values <- out.testing$values + unlist(calibration.error[sample(1:length(calibration.error),1)])
      out.testing <- out.testing[, c("names", "values")]
      names(out.testing) <- c("names", as.character(ind))
      testing.data.frame <- left_join(testing.data.frame, out.testing, by="names")
    }
    regions <- testing.data.frame$names
    testing.data.frame <- subset(testing.data.frame, select = -c(label, names))
    testing.data.frame <- data.frame(t(testing.data.frame ))
    names(testing.data.frame) <- regions
    testing.data.frame$Sex <- sub_data$Sex[btsp.idx.testing]
    testing.data.frame$Age <- sub_data$Age[btsp.idx.testing]
    testing.data.frame <- rbind(testing.data.frame, calib.i)
    
    list.res <- list()
    for(ROI in regions){
      formulai <- as.formula(paste0(ROI," ~ Sex"))
      m0 <- lm(formulai, data=testing.data.frame)
      estimate <- m0$coefficients[2]
      stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
      list.res[[ROI]] <- data.frame(Estimate = estimate, StdError = stdError)
    }
    imputation.list[[i.imp]] <- list.res
  }  
  imputation_result <- lapply(regions, region.list ,data.list = imputation.list)
  
  # imputation_result <- do.call(rbind, imputation_result)
  # saveRDS(imputation_result, paste0(DIR, "/multipleImp.rds"))
  
  # 2. Stratify by sex
  imputation.sex.list <- list()
  calibration.error.female.original <- calibration.error.original[which(sub_data$Sex[btsp.idx.calib]=="F")]
  
  calibration.error.male.original <- calibration.error.original[which(sub_data$Sex[btsp.idx.calib]=="M")]

  
  ind.sex <- sub_data$Sex[btsp.idx.testing]
  
  for(i.imp in 1:n.imp){
    testing.data.frame <- roi.df
    ## resample calibration error values
    calibration.error.female <- calibration.error.female.original[
      sample(1:length(calibration.error.female.original), 
             length(calibration.error.female.original), replace=TRUE)]
    
    calibration.error.male <- calibration.error.male.original[
      sample(1:length(calibration.error.male.original), 
             length(calibration.error.male.original), replace=TRUE)]
    
    for(ind in 1:length(testing.data)){
      
      out.testing <- testing.data[[ind]]
      if(ind.sex[ind] == "F"){
        out.testing$values <- out.testing$values + unlist(calibration.error.female[sample(1:length(calibration.error.female),1)])
      } else {
        out.testing$values <- out.testing$values + unlist(calibration.error.male[sample(1:length(calibration.error.male),1)])
      }
      out.testing <- out.testing[, c("names", "values")]
      names(out.testing) <- c("names", as.character(ind))
      testing.data.frame <- left_join(testing.data.frame, out.testing, by="names")
    }
    regions <- testing.data.frame$names
    testing.data.frame <- subset(testing.data.frame, select = -c(label, names))
    testing.data.frame <- data.frame(t(testing.data.frame ))
    names(testing.data.frame) <- regions
    testing.data.frame$Sex <- sub_data$Sex[btsp.idx.testing]
    testing.data.frame$Age <- sub_data$Age[btsp.idx.testing]
    testing.data.frame <- rbind(testing.data.frame, calib.i)
    
    list.res.robust <- list.res <- list()
    for(ROI in regions){
      formulai <- as.formula(paste0(ROI," ~ Sex"))
      m0 <- lm(formulai, data=testing.data.frame)
      estimate <- m0$coefficients[2]
      stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
      list.res[[ROI]] <- data.frame(Estimate = estimate, StdError = stdError)
    }
    imputation.sex.list[[i.imp]] <- list.res
  }  
  imputation_result_sex <- lapply(regions, region.list ,data.list = imputation.sex.list)

  # imputation_result_sex <- do.call(rbind, imputation_result_sex)
  # 
  # saveRDS(imputation_result_sex, paste0(DIR, "/multipleImpSex.rds"))
  all.res <- list(complete = list.res.copmlete, meanImp = list.res.mean, multiImp = imputation_result, multiImpStrat = imputation_result_sex)
  saveRDS(all.res, paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_region_mar/M7F3/simulation", seed.i, ".rds"))
}

mclapply(1:500, btsp_imp_df, mc.cores=20, mc.preschedule = FALSE)

## function: calculate metrics: bias, se, etc for all methods

## function: input: dataset, output: plot data frame

background <- readNifti("/home/xiongj3/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz")
mask.wm <- readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask3.nii.gz")
mask.wm[!mask.wm %in% row_wm$label] <- 0

plot.brain.dfcreate <- function(param.data, mask = mask.wm, bg = background,
                                region_match_df=row_wm, slices=c(35 ,40, 45, 50, 55,60)){
  names(region_match_df) <- c("value", "regions")
  param.data <- data.frame(param.data)
  colnames(param.data) <- "Estimate"
  
  new_data <- data.frame(param.data, region_match_df)
  FA_atlas_new <- mask
  for(i in 1:nrow(new_data)){
      FA_atlas_new[which(mask==new_data$value[i])] <- new_data$Estimate[i]
  }
  
  bg.list <- value.list <- list()
  for(i in 1:length(slices)){
    slicei <- slices[i]
    # Slice white matter background
    bg_df <- reshape2::melt(bg[,,slicei])  # Converts to long format with Var1, Var2, and value
    bg_df <- bg_df[bg_df$value>0,]
    colnames(bg_df) <- c("x", "y", "intensity")
    bg_df$z <- slicei
    bg.list[[i]] <- bg_df
    
    # Convert matrices to data frames with correct coordinates
    # Take slice
    value_df <- reshape2::melt(FA_atlas_new[,,slicei])
    mask_df <- reshape2::melt(mask[,,slicei])
    # carve out non-zero parts
    value_df <- value_df[mask_df$value>0,]
    mask_df <- mask_df[mask_df$value>0,]
    # add region name to dataset
    mask_df <- left_join(mask_df, region_match_df, by="value")
    value_df$region <- mask_df$region
    # Rename columns for clarity
    colnames(value_df)[1:3] <- c("x", "y", "intensity")
    value_df$z <- slicei
    value.list[[i]] <- value_df
  }
  bg.df <- do.call(rbind, bg.list)
  value.df <- do.call(rbind, value.list)
  return(list(bg.df, value.df))
}



