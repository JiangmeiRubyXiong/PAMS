################################################################################
##                  ---------       ----------       -------------
## Simulation & -> |  Calib  |     | Training |     | Synthesized | <- Multiple 
## Overall Truth   |  Image  |     |   Image  |     |    Image    |  imputations
##                  ---------       ----------       -------------
##                     |                 |                |         
##                     |                 |                |         
##  Data          ------------------------------------------------
##  Processing   | truncate 0-1; Summary by region; select region |
##                ------------------------------------------------
##                   |                |                  |         
##                   |                |                  |    
##                ---------       ----------        -------------
##               |  Calib  |     | Training |      | Synthesized | 
##               |   df    |     |     df   |      |      df     | 
##                ---------       ----------        -------------
##                   \             /      \              /  
##                    \          /         \           /  
##                     \        /           \        /  
##                    ------------       ------------
##                   |   Linear   |     |   Linear   |
##                   | Regression |     | Regression |
##                    ------------       ------------
################################################################################

## Data Processing: Calibration data
## Instead of sample first, and process later (in multiple imputation),
## Here we process first, and sample later

## Processing for true data in testing set (testing + calibration)
## Truncate
FA_truth_list <- readRDS(file="/media/disk2/beijing_dti/Neuroimaging/FA_truth_list.rds")
FA_truth_trunc <- lapply(FA_truth_list, norm.datavec)
# Summary by region and select region
FA_mask <- readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask2.nii.gz")
testt_subject_data <- read.csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")

mapply(function(dataset, i, mask=FA_mask, test.data=testt_subject_data){
  ## file names for the truncated files
  mask[which(mask!=0)] <- dataset
  id <- i
  id.old <- test.data$original_row[i+1]
  file.name <- 
    paste0("/media/disk2/beijing_dti/enhanced/unpack/",sprintf("%03d", id.old),"/scalars_standard/DTI_FA_resliced_trunc.nii.gz")
  ## write Nifti for the truncated files
  writeNifti(mask, file = file.name)
}, FA_truth_trunc, 0:88)
### FAWMed them
CalibBrain <- CalibBraincsv <- c()
for(i in 0:88){
  id <- i
  id.old <- testt_subject_data$original_row[i+1]
  CalibBrain <- c(CalibBrain,
    paste0("/media/disk2/beijing_dti/enhanced/unpack/",sprintf("%03d", id.old),"/scalars_standard/DTI_FA_resliced_trunc.nii.gz"))
  CalibBraincsv <- c(CalibBraincsv, paste0("/media/disk2/beijing_dti/enhanced/Calib_trunc/trunc",sprintf("%03d", id),".csv"))
}
write.table(CalibBrain, file = "/media/disk2/beijing_dti/Neuroimaging/listCalibFiles.txt", row.names = FALSE, col.names = FALSE)
write.table(CalibBraincsv, file = "/media/disk2/beijing_dti/Neuroimaging/listCalibcsv.txt", row.names = FALSE, col.names = FALSE)



btsp_simulation_calib <- function(seed.i, data.dir="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain_2mm",
                            phenotype=FALSE, n.imp=10){
  set.seed(seed.i)
  DIR="/media/disk2/beijing_dti/Neuroimaging/imp_res"
  # bootstrap 200 out of 89
  # read data
  testt_subject_data <- read_csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")
  btsp.idx <- sample(1:89, 200, replace=TRUE)
  # split 100 as calibration, 100 as testing
  btsp.idx.calib <- btsp.idx[1:100]
  
}

