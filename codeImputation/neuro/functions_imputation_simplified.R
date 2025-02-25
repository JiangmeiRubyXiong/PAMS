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

testt_subject_data <- read.csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv") # generated from python codes
testt_subject_data$NewID <- as.character(testt_subject_data$NewID) # the identifier of each object
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
  
  # identifier: the file name for each individual at /media/disk2/beijing_dti/enhanced/unpack, a three digit character
  new.id <- testt_subject_data$NewID[i+1]
  FA_synth_list[[new.id]] <- FA_synth
  FA_truth_list[[new.id]] <- FA_truth
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
# Step 1: list all bilateral regions
# Step 2: go back to the mask, and check which regions are still in the mask
# Step 3: recode the mask to remove bilateral differences
# Step 4: run FAWM with the new mask

# Function to convert to proper case
format_camel_case <- function(input_string) {
  parts <- strsplit(input_string, "_")[[1]]  # Split the string by underscores
  parts <- tolower(parts)                   # Convert all parts to lowercase
  parts[1] <- paste0(toupper(substr(parts[1], 1, 1)), substr(parts[1], 2, nchar(parts[1]))) # Capitalize first word
  paste(parts, collapse = "_")              # Combine parts back with underscores
}

regions_all <- read.table("/media/disk2/beijing_dti/Neuroimaging/Eve_Atlas-master/JHU_MNI_SS_WMPM_Type-III_SlicerLUT.txt")
names(regions_all) <- c("label", "name", "c1", "c2", "c3", "c4")
regions_all$name <- unlist(lapply(regions_all$name, format_camel_case))
list_names <- strsplit(regions_all$name, split="_")
regions_all$LR <- unlist(lapply(list_names, last))

regions_all$LR <- ifelse(regions_all$LR %in% c("left", "right"), regions_all$LR, "All")

region_pair <- unlist(lapply(list_names, function(x){
  if(last(x) %in% c("left", "right"))x=x[-length(x)]
if(any(x==""))x=x[-which(x=="")]
paste0(x, collapse="_")}))

regions_all$name <- region_pair

regions <- regions_all
regions_WM <- c(1,3,4,5,6,7,8,9,10,11,14,15,16,18,19,20,64,65,66,67,68,69,70,72,73,74,75,76,77,79,80,81,39,40,51,52,53,99)
regions <- regions[regions$label %in% regions_WM,]

# sort through the data frame to pair the left and right parts
pairs <- table(regions$name)
singles <- names(pairs)[pairs==1]
singles_code <- regions$label[which(regions$name %in% singles)]
# add the pair back to the lists
for(i in singles_code){
  idx <- which (regions_all$label==i)
  namei <- regions_all$name[idx]
  sidei <- regions_all$LR[idx]
  rowidx <- which(regions_all$name == namei & regions_all$LR != sidei)
  regions <- rbind (regions, regions_all[rowidx, ])
}

# In mask, assign all left side label to right side
# While at the same time, collect all regions that are present

FA_mask <- readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask2.nii.gz")
regions_names <- unique(regions$name)
mask_regions <- c()
for(namei in regions_names){
  left_num <- regions$label[regions$name==namei & regions$LR=="left"]
  right_num <- regions$label[regions$name==namei & regions$LR=="right"]
  if(any(FA_mask %in% c(left_num, right_num))){
    print(c(left_num, right_num))
    mask_regions <- c(mask_regions, namei)
    if(any(FA_mask == left_num)){FA_mask[FA_mask == left_num] <- right_num}
  }
}
writeNifti(FA_mask,"/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask3.nii.gz")
FA_mask <- readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask3.nii.gz")

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
# 
# 
# row_LR_pair <- data.frame(left=c( 8, 39, 10,
#                                   5, 15, 19,
#                                   11,  4, 16, 20,
#                                   7,   6, 
#                                   3, 14, 1, 18),
#                           right=c(70,  99, 72,
#                                   67, 76, 80,
#                                   73, 66, 77, 81, 
#                                   69, 68, 
#                                   65, 75, 64, 79), 
#                           names=c("ANGULAR",  'Cingulum_cingulate_gyrus', 'CUNEUS',
#                                   'INFERIOR_FRONTAL', 'INFERIOR_OCCIPITAL', 'INFERIOR_TEMPORA',
#                                   'LINGUAL', 'MIDDLE_FRONTAL', 'MIDDLE_OCCIPITAL', 'MIDDLE_TEMPORA',
#                                   'POSTCENTRAL', 'PRECENTRAL',
#                                   'SUPERIOR_FRONTAL', 'SUPERIOR_OCCIPITAL', 'SUPERIOR_PARIETAL', 'SUPERIOR_TEMPORAL'))
# row_LR_single <- data.frame(singles = c(74, 9,
#                                         52, 40, 51, 53),
#                             names= c('Fusiform_right', 'Precuneus_left',
#                                      "Body_of_corpus_callosum", "Cingulum_hippocampus_left",
#                                      "Genu_of_corpus_callosum_left", "Splenium_of_corpus_callosum_left"))
# 
# library(stringr)
# 
# 
# 
# row_wm <- list(row_LR_pair, row_LR_single)
# 

row_wm <- subset(regions, LR=="right"&name%in%mask_regions, c("label", "names"))
row_wm$names <- lapply(row_wm$names, function(string){gsub("[()]","", string)}) %>% unlist()
row_wm$names <- lapply(row_wm$names, function(string){gsub("-","_", string)}) %>% unlist()
saveRDS(row_wm, file="/media/disk2/beijing_dti/Neuroimaging/brain_ROI.rds")

#### Here redid FAWM with the new FA atlas
### wont need csv_filter since it is processed already

# csv_filter <- function(input_path, rows=row_wm, input_base ="/media/disk2/beijing_dti/Neuroimaging/imputation_csv",
#                        output_base = "/media/disk2/beijing_dti/Neuroimaging/imputation_rds"){
#   
#   df_csv <- read.delim(input_path, sep="", header = FALSE, stringsAsFactors = FALSE)
#   df_csv_mean <- df_csv[-(1:2), 1:2]
#   names(df_csv_mean) <- c("LabelID", "Mean")
#   
#   # combine the region names
#   # row_left <- rows[[1]]$left;row_right <- rows[[1]]$right; 
#   # df_csv_mean_LR <- (as.numeric(unlist(subset(df_csv_mean, LabelID %in% row_left, Mean)))+
#   #                      as.numeric(unlist(subset(df_csv_mean, LabelID %in% row_right, Mean))))/2
#   # df_csv_mean_S <- as.numeric(unlist(subset(df_csv_mean, LabelID %in% rows[[2]]$singles, Mean)))
#   # df_csv_out <- data.frame(names=c(rows[[1]]$names, rows[[2]]$names),
#   #                          value=c(df_csv_mean_LR, df_csv_mean_S))
#   
#   # Derive the output file path by replacing the base directory and extension
#   relative_path <- sub(input_base, "", input_path) # Get relative path from the base directory
#   output_path <- paste0(output_base, dirname(relative_path), "/",
#                            sub("\\.csv$", ".rds", basename(relative_path)))
#   
#   # Create the output directory if it doesn't exist
#   output_dir <- dirname(output_path)
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }
# 
#   saveRDS(df_csv_out, output_path)
# }
# files_brain_out <- read.table("/media/disk2/beijing_dti/Neuroimaging/listBraincsv.txt")[[1]]
# # Apply the function to all CSV files
# mclapply(files_brain_out, csv_filter, mc.cores=30, mc.preschedule = FALSE)


##### sex and age variable in all simulations
testt_subject_data <- read_csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")
btsp_simulation_sex <- function(seed.i, data.dir="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain_2mm",
                            phenotype=FALSE, n.imp=10, sub_data = testt_subject_data){
  set.seed(seed.i)
  DIR="/media/disk2/beijing_dti/Neuroimaging/imp_res"
  # bootstrap 200 out of 89
  # read data
  testt_subject_data <- read_csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")
  btsp.idx <- sample(1:89, 200, replace=TRUE)
  # split 100 as calibration, 100 as testing
  btsp.idx.calib <- btsp.idx[1:100]
  
  btsp.idx.testing <- btsp.idx[101:200]
  return(sub_data[btsp.idx.testing, c("Sex", "Age")])
  }

sex.var <- lapply(1:200, btsp_simulation_sex)



# This function is applied to one simulation one type of MULTIPLE imputation
# Imput: simulation number, input type, 
# Combine testing, and calibration
# Regression for all regions
# returns coefficients and variance


row_wm <- readRDS("/media/disk2/beijing_dti/Neuroimaging/brain_ROI.rds")
## keep only region of interest
read_brain_csv <- function(input_path, roi_df = row_wm){
    # read in data
    df_csv <- read.delim(input_path, sep="", header = FALSE, stringsAsFactors = FALSE)
    df_csv_mean <- df_csv[-(1:2), 1:2]
    df_csv_mean <- data.frame(lapply(df_csv_mean, as.numeric))
    names(df_csv_mean) <- c("label", "value")
    # filter out regions not of interest
    df_csv_mean <- left_join(row_wm, df_csv_mean, by="label")
    return(df_csv_mean[,c("names", "value")])
}

singleImpData <- function(path.data, sex.vec=sex.var, sim.i){
  # read all individual data in the same simulation
  files_imp <- list.files(path = path.data, full.names = TRUE)
  list.df <- lapply(files_imp, read_brain_csv)
  list.df <- data.frame(t(do.call(cbind, lapply(list.df, "[[", 2))))

  # read calibration data
  calibi <- readRDS(paste0("/media/disk2/beijing_dti/enhanced/Calib_df/calib_sim_",sprintf("%03d", sim.i),".rds"))
  list.df <- data.frame(list.df, sex.vec[[sim.i]])
  names(list.df) <- names(calibi)
  list.df <- rbind(list.df, calibi)
  names_brain <- names(list.df)[1:(ncol(list.df)-2)]
  
  list.res <- list()
  for(ROI in names_brain){
    formulai <- as.formula(paste0(ROI," ~ Sex + Age"))
    list.res[[ROI]] <- summary(lm(formulai, data=list.df))[["coefficients"]][2,1:2]
  }
  # regression analysis
  return(list.res)
}

# function: reorder the nested list in the multiple imputation lists,
# then combine the imputations results
# and do rubin's rule calculation
source("~/ImageImputation/dfCalc.R")
region.list <- function(roi, data.list, n.imp=10){
  # reorder the nested list in the multiple imputation lists
  lists <- lapply(data.list, "[[", roi)
  # combine the imputations results
  list.df <- do.call(rbind,lists)
  # rubin's rule calculation
  means <- list.df[,1]
  ses <- list.df[,2]
  VARB <- var(means)
  vars <- sqrt(var(means)*(n.imp+1)/n.imp+mean(ses^2))
  stats <- qt(0.975, dfCalc(VARB, vars, 10 , 200, 2))
  return(data.frame(mean=mean(means), sd=vars, stats=stats))
  }


# this function deals with all simulations in the analysis
lm_results_neuro_simulation <- function(sim.idx){
  # find the location of simulation folder of data frames
  DIR=paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_csv/simulation", sim.idx)
  # do the same analysis for a set of data
 
  hier_list_results <- lapply(list.dirs(DIR)[3:12], singleImpData, sim.i = sim.idx)
  local_list_results <- lapply(list.dirs(DIR)[14:23], singleImpData, sim.i = sim.idx)
  simple_list_results <- lapply(list.dirs(DIR)[26:35], singleImpData, sim.i  = sim.idx)
  mean_list_results <- singleImpData(list.dirs(DIR)[24], sim.i = sim.idx)
  
  # Rubin's rule for multiple imputation results
  brain_roi <- names(hier_list_results[[1]])
  # original results list: imputation, regions ; reorder to regions, imputation
  hier_list_results <- lapply(brain_roi, region.list ,data.list = hier_list_results)
  local_list_results <- lapply(brain_roi, region.list ,data.list = local_list_results)
  simple_list_results <- lapply(brain_roi, region.list ,data.list = simple_list_results)
  
  res.list <- list(meanImp = mean_list_results, simpleImp = simple_list_results,
                   localImp = local_list_results, hierImp = hier_list_results)
  saveRDS(res.list, file=paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/simulation", sim.idx, ".rds"))
}


#### function: put data to nifti plots

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

## function:  generate and combine dataframe, and plot wrap 6 slices and 4 methods
brain.plot <- function(){

}




