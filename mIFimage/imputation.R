library(tools)
library(raster)
library(terra)
library(readr)
library(dplyr)
library(parallel)
library(sandwich)
library(nlme)
source("~/ImageImputation/dfCalc.R")

# normalization methods
#  normalized_array = np.log10(1 + (img_array / mean_intensity))

# scale the narmalized image back to synthesized images

# read data

synth_wd <- "/media/disk2/GCA_norm_inference_baseline"
normed_wd <- "/media/disk2/GCA_norm"
######################################################################################
# construct data frame of patient and image information
test.IDs <- list.files(synth_wd)
test.ID.id <- lapply(strsplit(test.IDs, "_"),"[[", 1)
test.IDs.tissue <- lapply(strsplit(test.IDs, "_"), function(x){if(length(x)==2)return(x[[2]])else{return(NA)}})
test.patient <- substr(test.IDs, 1, 6)
test.tissue <- substr(test.IDs, 7, 8)
test.tissue.full <- substr(test.IDs, 7, 9)
test.ID.df <- data.frame(ID = test.IDs,testID = unlist(test.ID.id), tissue = unlist(test.IDs.tissue),
                         patient = test.patient, tissue = test.tissue,
                         tissue.full = test.tissue.full)

markers <- c('CD11B', 'CD20', 'CD3D', 'CD4', 'CD45', 'CD68', 'CD8', 'CGA', 
             'ERBB2', 'FOXP3', 'LYSOZYME', 'OLFM4', 'PCNA', 'SOX9')
set_number <- c(1,2,2,3,3,3,3,6,6,6)

status <- c('inactive', 'active', 'active', 'inactive', 'active', 'normal', 
            'active', 'active', 'active', 'active')

patient <- c('GCA004', 'GCA045', 'GCA059', 'GCA069', 'GCA062', 'GCA072', 'GCA077', 
             'GCA113', 'GCA118', 'GCA132')
synth_ids <- data.frame(patient, status, set_number)
test.ID.df <- left_join(test.ID.df, synth_ids, by="patient")


## combine with test.ID.df 
# calib.df <- data.frame(patient = c(calib.patient, test.patient),
#                        testSplit = c(mapply(rep, c("calib", "test"),
#                                             c(sample.size, length(patients_unique)-sample.size))))
# test.ID.df <- left_join(test.ID.df, calib.df, by="patient")
test.ID.list <- split(test.ID.df, test.ID.df$patient)
saveRDS(test.ID.df, file="~/ImImp/codeImputation/pixn2n/test.ID.df.rds")
saveRDS(test.ID.list, file="~/ImImp/codeImputation/pixn2n/test.ID.list.rds")
test.ID.df$status <- ifelse(test.ID.df$status=="normal", "inactive", test.ID.df$status)
######################################################################################
markers <- c('CD11B', 'CD20', 'CD3D', 'CD4', 'CD45', 'CD68', 'CD8', 'CGA', 
             'ERBB2', 'FOXP3', 'LYSOZYME', 'OLFM4', 'PCNA', 'SOX9')
test.ID.df <- readRDS("~/ImImp/codeImputation/pixn2n/test.ID.df.rds")
test.ID.list <- readRDS("~/ImImp/codeImputation/pixn2n/test.ID.list.rds")
# linear transformation from [c,d] to [a,b] in the transformed image: a+[(pixel-c)(b-a)/(d-c)]
linear_transform <- function(input_vec, ref_vec){
  ref_range <- range(ref_vec)
  a <- ref_range[1]; b <- ref_range[2]
  input_range <- range(input_vec)
  c <- 0; d <- 255
  
  output_vec <- a + (input_vec-c)*(b-a)/(d-c)
  
  perc <- c(0.90, 0.95, 0.99)
  perc_out <- c(mean(output_vec), quantile(output_vec, perc))
  perc_ref <- c(mean(ref_vec), quantile(ref_vec, perc))
  return(list(input = perc_out, ref = perc_ref))
}
######################################################################################
# Level 1: Apply to all markers (each list is one marker)
# level 2: Apply to all individuals (tissue)
### TIB and ACB tissues: belonging to each patient, use as different covariates

## read in all values, scale & take 90% perc
##### function: construct file path for Marker A
##### Start from synthsized image, get the tissue number (if any)
##### and go back to the normalized path
path_gen <- function(id.idx, ID.df = test.ID.df, marker){
  IDi <- ID.df[id.idx, "testID"]
  Tissuei <- ID.df[id.idx, "tissue"]
  # synth path
  if(is.na(Tissuei)){
    syn.path <- paste0("/media/disk2/GCA_norm_inference_baseline/", IDi, 
                       "/RECONSTRUCTED/GCA_norm_v2_N_to_N_lose_5_marker_instance_norm_1024p_pix2pixHD2/",
                       IDi, "_", marker, "_reconstructed.tif"
    )
  }else{
    syn.path <- paste0("/media/disk2/GCA_norm_inference_baseline/", IDi,"_",Tissuei, 
                       "/RECONSTRUCTED/GCA_norm_v2_N_to_N_lose_5_marker_instance_norm_1024p_pix2pixHD2/",
                       IDi,"_",Tissuei,"_", marker, "_reconstructed.tif"
    )
  }
  # norm path
  setid <- ID.df[id.idx, "set_number"]
  if(is.na(Tissuei)){
    norm.path <- paste0("/media/disk2/GCA_norm/Set0", setid, "/", IDi, "/AFRemoved/",
                        IDi, "_", marker, ".tif"
    )
  }else{
    norm.path <- paste0("/media/disk2/GCA_norm/Set0", setid, "/", IDi, "/AFRemoved/",
                        IDi, "_",Tissuei,"_", marker, ".tif"
    )
  }
  
  # mask path
  if(is.na(Tissuei)){
    mask.path <- paste0("/media/disk2/GCA_norm_inference_baseline/", IDi, 
                        "/MASKS/", IDi, "_TISSUE_MASK.tif"
    )
  }else{
    mask.path <- paste0("/media/disk2/GCA_norm_inference_baseline/", IDi,"_",Tissuei, 
                        "/MASKS/", IDi,"_",Tissuei,"_TISSUE_MASK.tif"
    )
  }
  return(list(synth = syn.path, norm = norm.path, mask = mask.path))
}
## read data (and remove 0s)
extract_pixel_intensity_raster <- function(file_path, mask_path){
  # Load the .tif file
  img <- rast(file_path)
  img <- as.vector(img)
  
  mask <- rast(mask_path)
  mask <- as.vector(mask)
  
  img <- img[which(mask>0)]
  
  return(img)
}

## extract quantile from both truth and synth
scale_quantile <- function(data.idx, marker){
  path_3 <- path_gen(data.idx, marker=marker)
  path_img <- path_3[1:2]
  path_mask <- path_3[[3]]
  synth_norm_vec_list <- lapply(path_img, extract_pixel_intensity_raster, mask_path = path_mask)
  perc_0 <- linear_transform(synth_norm_vec_list[[1]], synth_norm_vec_list[[2]])
  return(perc_0)
}

## Run for one marker, one percentile
analysis_marker <- function(marker){
  quantilessynth <- lapply(1:nrow(test.ID.df),scale_quantile, marker=marker)
  saveRDS(quantilessynth, file=paste0("/media/disk2/GCA_norm/tissue_quantiles/",marker, ".rds"))
}


## apply to all marker and save the result
# run in screen
# mclapply(markers, analysis_marker, mc.cores=5, mc.preschedule = FALSE)

## multiple imputation
# For each marker
# 1. split data into calibration and testing
# 2. extract empirical errors (per tissue)
# 3. sample errors and add to each testing tissue
# 4. calculate rubin's rule estimation results (GEE)



marker_imputation <- function(marker, seed.i=1, n.imp = 20, type="strat"){
  ## split data into calibration and testing
  # read calibration and testing split
  test.ID.df <- readRDS("~/ImImp/codeImputation/pixn2n/test.ID.df.rds")
  # read data
  data.marker <- readRDS(paste0("/media/disk2/GCA_norm/tissue_quantiles/", marker, ".rds"))
  ## bootstrap
  ## calibration/test split: sample tissue
  sample.size <- 200
  set.seed(seed.i)
  all.patients <- sample(1:length(data.marker), sample.size, replace = TRUE)
  calib.patient <- all.patients[1:(sample.size/2)]
  test.patient <- all.patients[(sample.size/2+1): sample.size]
  
  ## assemble calibration and testing data
  calib_data_list <- data.marker[calib.patient]
  test_data_list <- data.marker[test.patient]
  
  # extract: truth & synth, calibration 
  calib_data_truth <- do.call(rbind, lapply(calib_data_list, "[[", "ref"))
  calib_data_synth <- do.call(rbind, lapply(calib_data_list, "[[", "input"))
  # extract: synth, calibration
  test_data_synth <- do.call(rbind, lapply(test_data_list, "[[", "input"))
  # take diff
  empirical_error <- calib_data_truth - calib_data_synth
  
  colnames(calib_data_truth) <- c("mean", "q90", "q95", "q99")
  colnames(calib_data_synth) <- c("mean", "q90", "q95", "q99")
  
  # check
  calib_compare <- data.frame(truth90 = calib_data_truth[,"q90"],
                              synth90 = calib_data_synth[,"q90"],
                              test.ID.df[calib.patient, c("status", "tissue.1", "patient")])
  # sample calibration error
  imputed_list <- imputed_list_robust <- list()
  for (i in 1:n.imp){
    
    if(type!="simple"){
      tissues <- test.ID.df$status[calib.patient]
      tissues_test <- test.ID.df$status[test.patient]
      empirical.TI.idx.o <- which(tissues == "active")
      empirical.AC.idx.o <- which(tissues !="active")
      # resample idx
      empirical.TI.idx <- sample(empirical.TI.idx.o, length(empirical.TI.idx.o), replace = TRUE)
      empirical.AC.idx <- sample(empirical.AC.idx.o, length(empirical.AC.idx.o), replace = TRUE)
      
      sample_idx <- rep(NA, nrow(test_data_synth))
      sample_idx[which(tissues_test=="active")] <- sample(empirical.TI.idx, 
                                                      sum(tissues_test=="active"), 
                                                      replace = TRUE)
      
      sample_idx[which(tissues_test!="active")] <- sample(empirical.AC.idx, 
                                                      sum(tissues_test!="active"),
                                                      replace = TRUE)
      
    }else{
      # resample idx
      empirical_error_idx <- sample(1:nrow(empirical_error), nrow(empirical_error), replace=TRUE)
      
      sample_idx <- sample(empirical_error_idx, nrow(test_data_synth), replace = TRUE)
    }
    error_sampled <- empirical_error[sample_idx,]
    imputed <- data.frame(test_data_synth + error_sampled)
    names(imputed) <- colnames(calib_data_truth)
    
    # use linear regression model
    res_list <- list()
    for(type.i in names(imputed)){
      ## calibration data frame and test data frame rbind()
      imputed_i_df <- data.frame(values = c(calib_data_truth[,type.i], imputed[[type.i]]),
                                 test.ID.df[all.patients, c("status", "tissue.1", "patient")])
      
      # combine status
      imputed_i_df$status <- ifelse(imputed_i_df$status=="normal", "inactive", imputed_i_df$status)
      imputed_i_df$status <- factor(imputed_i_df$status, levels=c("inactive", "active"))
      
        m0 <- lm(values ~ status, imputed_i_df)
        estimate <- m0$coefficients[2]
        stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
        res_list[[type.i]] <- list(coef = estimate, sd = stdError, type = type.i)

      
    }
    imputed_list[[i]] <- res_list
  }
  
  ## Rubin's rule
  # extract patient outcome
  type_output <- list()
  for(i in names(imputed_list[[1]])){
    imputed_i_list <- lapply(imputed_list, "[[", i)
    imputed_i_res_df <- data.frame(lapply(names(imputed_i_list[[1]]),
                                          function(nm){unlist(lapply(imputed_i_list, "[[", nm))}))
    names(imputed_i_res_df) <- names(imputed_i_list[[1]])
    type_output[[i]] <- rubin_imp_calc(imputed_i_res_df)
  }
  type_output_df <- do.call(rbind, type_output)
  return(type_output_df)
}

## Rubin's rule function (for one type tissue, and for one summary data type)
rubin_imp_calc <- function(list.df, sample.size = 200){
  n.imp <- nrow(list.df)
  # rubin's rule calculation
  means <- list.df[["coef"]]
  ses <- list.df[["sd"]]
  VARB <- var(means)
  vars <- sqrt(var(means)*(n.imp+1)/n.imp+mean(ses^2))
  stats <- qt(0.975, dfCalc(VARB, vars, n.imp , sample.size, 2))
  return(data.frame(mean=mean(means), sd=vars, stats=stats, type = list.df$type[1]))
}



## Mean imputation & truth at every simulation ( and complete )
mean_imputation <- function(marker, seed.i=1){
  ## split data into calibration and testing
  # read calibration and testing split
  test.ID.df <- readRDS("~/ImImp/codeImputation/pixn2n/test.ID.df.rds")
  # read data
  data.marker <- readRDS(paste0("/media/disk2/GCA_norm/tissue_quantiles/", marker, ".rds"))
  ## bootstrap
  ## calibration/test split: sample tissue
  sample.size <- 200
  set.seed(seed.i)
  all.patients <- sample(1:length(data.marker), sample.size, replace = TRUE)
  calib.patient <- all.patients[1:(sample.size/2)]
  test.patient <- all.patients[(sample.size/2+1): sample.size]
  
  ## assemble calibration and testing data
  calib_data_list <- data.marker[calib.patient]
  test_data_list <- data.marker[test.patient]
  
  # extract: truth, calibration 
  calib_data_truth <- do.call(rbind, lapply(calib_data_list, "[[", "ref"))
  # extract: synth, test
  test_data_synth <- do.call(rbind, lapply(test_data_list, "[[", "input"))
  # extract: truth, test
  test_data_truth <- do.call(rbind, lapply(test_data_list, "[[", "ref"))
  
  all.simulation <- rbind(calib_data_truth, test_data_synth)
  colnames(all.simulation) <- c("mean", "q90", "q95", "q99")
  
  truth.simulation <- rbind(calib_data_truth, test_data_truth)
  colnames(truth.simulation) <- colnames(calib_data_truth) <- c("mean", "q90", "q95", "q99")
  
  # use linear regression model
  truth_list <- res_list <- complete_list <- list()
  for(type.i in colnames(all.simulation)){
    ## calibration data frame and test data frame rbind()
    imputed_i_df <- data.frame(values = all.simulation[,type.i],
                               test.ID.df[all.patients, c("status", "tissue.1", "patient")])
    # combine status
    imputed_i_df$status <- ifelse(imputed_i_df$status=="normal", "inactive", imputed_i_df$status)
    imputed_i_df$status <- factor(imputed_i_df$status, levels=c("inactive", "active"))
    m0 <- lm(values ~ status, imputed_i_df)
    estimate <- m0$coefficients[2]
    stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
    res_list[[type.i]] <- data.frame(mean = estimate, 
                               sd = stdError, stats=1.96, type = type.i)
    
    imputed_i_df$values <- truth.simulation[,type.i]
    m0 <- lm(values ~ status, imputed_i_df)
    estimate <- m0$coefficients[2]
    stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
    truth_list[[type.i]] <- data.frame(data.frame(mean = estimate, 
                                                  sd = stdError, stats=1.96, type = type.i))
    # complete data
    complete_i_df <- data.frame(values = calib_data_truth[,type.i],
                               test.ID.df[calib.patient, c("status", "tissue.1", "patient")])
    complete_i_df$status <- ifelse(complete_i_df $status=="normal", "inactive", complete_i_df$status)
    complete_i_df$status <- factor(complete_i_df $status, levels=c("inactive", "active"))
    m0 <- lm(values ~ status, complete_i_df)
    estimate <- m0$coefficients[2]
    stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
    complete_list[[type.i]] <- data.frame(mean = estimate, 
                                     sd = stdError, stats=1.96, type = type.i)
    
  }
  res_list_df <- do.call(rbind, res_list)
  truth_list_df <- do.call(rbind, truth_list)
  complete_list_df <- do.call(rbind, complete_list)
  out.list <- list(meanImp= res_list_df, truthImp = truth_list_df, complete = complete_list_df)
  return(out.list)
}

# repeat the same mean_imputation for 200 simulations
## format result data
# Define the function to run for each seed
run_imputation_mkr <- function(seed_i) {
  markers <- c('CD11B', 'CD20', 'CD3D', 'CD4', 'CD45', 'CD68', 'CD8', 'CGA', 
               'ERBB2', 'FOXP3', 'LYSOZYME', 'OLFM4', 'PCNA', 'SOX9')
  mkr.truth <- mkr.mean <- mkr.multi <- mkr.strat <- mkr.complete <- list()
  for(mkr in markers){
    res.mean <-  mean_imputation(marker = mkr, seed.i = seed_i)
    mkr.truth[[mkr]] <- res.mean[["truthImp"]]
    mkr.mean[[mkr]] <-  res.mean[["meanImp"]]
    mkr.complete[[mkr]] <-  res.mean[["complete"]]
    mkr.multi[[mkr]] <-  marker_imputation(marker = mkr, seed.i = seed_i, type="simple", n.imp = 50)
    mkr.strat[[mkr]] <-  marker_imputation(marker = mkr, seed.i = seed_i, type="strat", n.imp = 50)
  }
  out.list <- list(mkr.truth, mkr.mean, mkr.complete, mkr.multi, mkr.strat)
  return(out.list)
}

## Run this
res.list.marker <- mclapply(1:500, run_imputation_mkr, mc.cores = 20, mc.preschedule = FALSE)  # Adjust based on total cores
saveRDS(res.list.marker, file="/media/disk2/GCA_norm/sim_result_updated.rds")

res.list.marker <- mclapply(1:500, run_imputation_mkr, mc.cores = 20, mc.preschedule = FALSE)  # Adjust based on total cores
saveRDS(res.list.marker, file="/media/disk2/GCA_norm/sim_result_m50.rds")
