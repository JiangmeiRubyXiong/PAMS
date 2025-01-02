library(oro.nifti)
library(neurobase)
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


# get errors
FA_diff <- list() # errors for all calibration data
FA_atlas=c(readnii("/home/xiongj3/fsl/data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz"))
for(i in 0:88){
  id <- i
  id.old <- testt_subject_data$original_row[i+1]
  FA_synth=c(readnii(paste0(data.dir, "/subject_", id, "_reassembled_brain.nii.gz")))
  FA_synth[is.na(FA_synth)] <- 0
  FA_truth=c(readnii(paste0("/media/disk2/beijing_dti/enhanced/unpack/",sprintf("%03d", id.old),"/scalars_standard/DTI_FA_resliced.nii.gz")))
  # remove background
  FA_synth <-  FA_synth[-which(FA_atlas==0)];FA_truth <-  FA_truth[-which(FA_atlas==0)]
  FA_diff[[i+1]] <- FA_synth-FA_truth
}
saveRDS(FA_diff, file="/media/disk2/beijing_dti/Neuroimaging/FA_diff.rds")
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


# one simulation
btsp_simulation <- function(seed.i=1, data.dir="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain_2mm"){
  set.seed(seed.i)
  # bootstrap 200 out of 89
  # read data
  testt_subject_data <- read_csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")
  btsp.idx <- sample(1:89, 200, replace=TRUE)
  # split 100 as calibration, 100 as testing
  btsp.idx.calib <- sample(btsp.idx, 100, replace=TRUE)
  btsp.idx.testing <- setdiff(btsp.idx, btsp.idx.calib)
  
  ###### use imimp
  
  imp.simple <- imimp(Yhat.test, calibYhat , calibY, cdeType = "simple", x.mar = 50, y.mar=100, n.mi=n.imp)
  imp.local<- imimp(Yhat.test, calibYhat, calibY, cdeType = "local", n.mi=n.imp)

}


# Downstream Analysis