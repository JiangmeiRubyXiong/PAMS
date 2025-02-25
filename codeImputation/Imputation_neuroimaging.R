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
################################################################################

################################################################################
## Data Processing: Calibration data (and testing true data) ###################
## Instead of sample first, and process later (in multiple imputation),  #######
## Here we process first, and sample later. ####################################
################################################################################

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
############################### FAWM Process ###################################
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
## Do not run: pick out only region of interest
# mclapply(CalibBraincsv, csv_filter, 
#          input_base ="/media/disk2/beijing_dti/enhanced/Calib_trunc",
#          output_base = "/media/disk2/beijing_dti/enhanced/Calib_trunc_rds",
#          mc.cores=30, mc.preschedule = FALSE)

## function to combine the selection of calibration data in each simulation
calib_data <- function(data.vec){
  data.vec <- data.vec-1
  data_files <- paste0("/media/disk2/beijing_dti/enhanced/Calib_trunc/trunc",sprintf("%03d", data.vec),".csv")
  data_all <- lapply(data_files,read_brain_csv)
  names_brain <- data_all[[1]][[1]]
  data_all <- data.frame(t(do.call(cbind, lapply(data_all, "[[", 2))))
  names(data_all) <- names_brain
  return(data_all)
}

btsp_calib_df <- function(seed.i, data.dir="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain_2mm",
                                phenotype=FALSE, n.imp=10, sub_data = testt_subject_data){
  set.seed(seed.i)
  DIR="/media/disk2/beijing_dti/Neuroimaging/imp_res"
  # bootstrap 200 out of 89
  # read data
  
  btsp.idx <- sample(1:89, 200, replace=TRUE)
  # split 100 as calibration, 100 as testing
  btsp.idx.calib <- btsp.idx[1:100]
  
  btsp.idx.testing <- btsp.idx[101:200]
  calib_data_i <- calib_data(btsp.idx.calib)
  calib_data_i$Sex <- sub_data$Sex[btsp.idx.calib]
  calib_data_i$Age <- sub_data$Age[btsp.idx.calib]
  saveRDS(calib_data_i,  file=paste0("/media/disk2/beijing_dti/enhanced/Calib_df/calib_sim_",sprintf("%03d", seed.i),".rds"))
}

lapply(1:200, btsp_calib_df)

calib_data_i <- calib_data((1:89))
calib_data_i$Sex <- testt_subject_data$Sex
calib_data_i$Age <- testt_subject_data$Age
saveRDS(calib_data_i,  paste0("/media/disk2/beijing_dti/enhanced/Calib_df/calib_all.rds"))

################################################################################
################################################################################
################################################################################

################################################################################
## Data Processing: Synthesized Data. ##########################################
## Sample first, and process later #############################################
################################################################################
## Multiple imputation
seq.list <- as.list(seq (1,200,1)) 
mclapply(seq.list, btsp_simulation, mc.cores=10, mc.preschedule = FALSE)
## Imputation processing: truncate imputation outcome,
## then save as nifiti file (for FAWM)
mclapply(seq.list, res_to_vec, mc.cores=10, mc.preschedule = FALSE)

############################### FAWM Process ###################################
# create list of files for atlas
files_brain <- list.files("/media/disk2/beijing_dti/Neuroimaging/imputation_brains/", recursive = TRUE, full.names = FALSE)
files_brain0 <- paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_brains/",files_brain)
files_brain_out <- gsub("\\.nii\\.gz$", "", files_brain)
files_brain_out <- paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_csv/",files_brain_out,".csv")
write.table(files_brain0, file = "/media/disk2/beijing_dti/Neuroimaging/listBrainFiles.txt", row.names = FALSE, col.names = FALSE)
write.table(files_brain_out, file = "/media/disk2/beijing_dti/Neuroimaging/listBraincsv.txt", row.names = FALSE, col.names = FALSE)
############################### FAWM Process ###################################
# Do not run
# ## select only region of interest,
# ## and average left&right brain region pairs
# files_brain_out <- read.table("/media/disk2/beijing_dti/Neuroimaging/listBraincsv.txt")[[1]]
# # Apply the function to all CSV files
# row_wm <- readRDS(file="/media/disk2/beijing_dti/Neuroimaging/brain_ROI.rds")
# mclapply(files_brain_out, csv_filter, mc.cores=30, mc.preschedule = FALSE)

## sex variable in all simulations
testt_subject_data <- read_csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")
sex.var <- lapply(1:200, btsp_simulation_sex)

################################################################################
## Linear Regression ###########################################################
## Piece together testing + calibration ########################################
################################################################################

mclapply(1:200, lm_results_neuro_simulation, mc.cores=20, mc.preschedule = FALSE)

## true data at every simulation
btsp_regression_df <- function(seed.i, data.dir="/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/synthetic_brain_2mm",
                          phenotype=FALSE, n.imp=10, sub_data = testt_subject_data){
  set.seed(seed.i)
  DIR="/media/disk2/beijing_dti/Neuroimaging/imp_res"
  # bootstrap 200 out of 89
  # read data
  
  btsp.idx <- sample(1:89, 200, replace=TRUE)
  # calculate true result for each simulation, and save at the same place as simulations
  testing_data_i <- calib_data(btsp.idx)
  rois <- names(testing_data_i)
  testing_data_i$Sex <- sub_data$Sex[btsp.idx]
  testing_data_i$Age <- sub_data$Age[btsp.idx]
  list.res <- list()
  for(ROI in rois){
    formulai <- as.formula(paste0(ROI," ~ Sex + Age"))
    list.res[[ROI]] <- summary(lm(formulai, data=testing_data_i))[["coefficients"]][2,1:2]
  }
  list.res <- do.call(rbind, list.res)
  saveRDS(list.res,  file=paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/truedata",seed.i,".rds"))
}
lapply(1:200, btsp_regression_df)

## True data no simulation
calib_data_i <- calib_data((1:89))
rois <- names(calib_data_i)
calib_data_i$Sex <- testt_subject_data$Sex
calib_data_i$Age <- testt_subject_data$Age
list.res <- list()
for(ROI in rois){
  formulai <- as.formula(paste0(ROI," ~ Sex + Age"))
  list.res[[ROI]] <- summary(lm(formulai, data=calib_data_i))[["coefficients"]][3,1:4]
}
list.res <- do.call(rbind, list.res)
saveRDS(list.res,  file="/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/truedata0.rds")
list.res <- as.data.frame(list.res)
list.res$regions <- as.character(rownames(list.res))
################################################################################
## Make Plot ###################################################################
################################################################################

## create atlas for plotting
FA_atlas_plot=readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_WMmask.nii.gz")
roi <- readRDS(file="/media/disk2/beijing_dti/Neuroimaging/brain_ROI.rds")
# Do not run # set non-ROI regions to be 0
# vals <- c(roi[[1]]$left, roi[[1]]$right, roi[[2]]$singles)
# FA_atlas_plot[which(!FA_atlas_plot %in% vals)] <- 0
# for(i in 1:nrow(roi[[1]])){
#   # set left label to be same as right
#   left.val <- roi[[1]]$left[i]; right.val <- roi[[1]]$right[i]
#   FA_atlas_plot[which(FA_atlas_plot==left.val)] <- right.val
# }
# writeNifti(FA_atlas_plot,"/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_plotMask.nii.gz" )
# region_match <- data.frame(value=c(roi[[1]]$right,roi[[2]]$singles), regions=c(roi[[1]]$names,roi[[2]]$names))
# library(dplyr)
# new_data <- left_join(list.res, region_match, by="regions")
# for(i in 1:nrow(new_data)){
#   FA_atlas_new <- FA_atlas_plot
#   FA_atlas_new[which(FA_atlas_plot==new_data$value[i])] <- new_data$Estimate[i]
# }
# writeNifti(FA_atlas_new, file = "/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/temp.nii.gz")

library(RNifti)
library(ggplot2)
library(raster)
library(sf)

## calculate metrics: bias, se, etc for all methods
all.data <- lapply(1:200, function(i){readRDS(file=paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/simulation",i,".rds"))})
all.data <- lapply(1:4, function(i){lapply(all.data, "[[", i)})
all.data.means <- lapply(all.data, function(lists){vector <- lapply(lists, do.call, what=rbind)
                                                   vector <- lapply(vector, function(matr){matr[,1]})
                                                   do.call(cbind, vector)})
all.data.mean <- lapply(all.data.means, rowMeans)
truedata0 <- readRDS("/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/truedata0.rds")
all.data.bias.mean <- lapply(all.data.mean, function(x, truemean=truedata0){x-truemean[,1]})

all.data.se <- lapply(all.data, function(lists){vector <- lapply(lists, do.call, what=rbind)
                                                vector <- lapply(vector, function(matr){matr[,2]})
                                                do.call(cbind, vector)})
# Mean SE with model
all.data.se.mean <- lapply(all.data.se, rowMeans)
# SE of estimator
row.se <- function(x){
  vec.sd <- list()
  for(i in 1:nrow(x)){
    vec.sd[[i]] <- sd(x[i,])
  }
  return(unlist(vec.sd))
}
all.data.est.se <- lapply(all.data.means, row.se)

all.data.ratio.2 <- mapply(function(upper, lower){return(upper/lower)},
                            all.data.se.mean, all.data.est.se, SIMPLIFY = FALSE)

# MSE
## Numerator
all.data.means <- lapply(all.data, function(lists){vector <- lapply(lists, do.call, what=rbind)
vector <- lapply(vector, function(matr){matr[,1]})
do.call(cbind, vector)})
mse.bias <-  function(x, truedata=truedata0){
  vec.mse <- list()
  x <- x - matrix(truedata[,1], nrow=nrow(truedata), ncol=ncol(x), byrow = FALSE)
  for(i in 1:nrow(x)){
    vec.mse[[i]] <- mean(x[i,]^2)
  }
  return(unlist(vec.mse))
}
all.data.mse.est <- lapply(all.data.means, mse.bias)
## Denominator
sim.res <- lapply(1:200, function(x){
  readRDS(paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/truedata",x,".rds"))[,1]})
sim.res <- do.call(cbind, sim.res)
sim.res.bias <- mse.bias(sim.res)
## ratio
all.data.ratio.3 <- lapply(all.data.mse.est,function(upper, lower=sim.res.bias){return(upper/lower)})

# Coverage
all.data.means <- lapply(all.data, function(lists){vector <- lapply(lists, do.call, what=rbind)
vector <- lapply(vector, function(matr){matr[,1]})
do.call(cbind, vector)})
all.data.bias.abs <- lapply(all.data.means, function(x, truedata=truedata0){
  abs(x - matrix(truedata[,1], nrow=nrow(truedata), ncol=ncol(x), byrow = FALSE))
})


all.data.stats <- lapply(all.data[2:4], function(lists){vector <- lapply(lists, do.call, what=rbind)
vector <- lapply(vector, function(matr){matr[,3]})
do.call(cbind, vector)})
all.data.stats <- append(list(matrix(1.96, nrow=21, ncol=200)), all.data.stats)

all.data.coverage <- mapply(function(stats, se, abs.bias){return(stats * se >= abs.bias)}, 
                            all.data.stats, all.data.se, all.data.bias.abs, SIMPLIFY = FALSE)
all.data.coverage.perc <- lapply(all.data.coverage, rowMeans)

## combine all 6 axis, and 4 different methods

plot.df.list <- lapply(all.data.bias, plot.brain.dfcreate)
bg.plot.df <- rbind(data.frame(plot.df.list[[1]][[1]], method="mean"),
                      data.frame(plot.df.list[[2]][[1]], method="simple"),
                      data.frame(plot.df.list[[3]][[1]], method="local"),
                      data.frame(plot.df.list[[4]][[1]], method="hier"))
bg.plot.df$method <- factor(bg.plot.df$method, levels=c("mean",
                                                        "simple",
                                                        "local",
                                                        "hier"))

mean.plot.df <- rbind(data.frame(plot.df.list[[1]][[2]], method="mean"),
                      data.frame(plot.df.list[[2]][[2]], method="simple"),
                      data.frame(plot.df.list[[3]][[2]], method="local"),
                      data.frame(plot.df.list[[4]][[2]], method="hier"))
mean.plot.df$method <- factor(mean.plot.df$method, levels=c("mean",
                                                        "simple",
                                                        "local",
                                                        "hier"))

p1 <- ggplot() +
  geom_tile(data = bg.plot.df, aes(x = x, y = y), fill = "grey90",alpha = 0.8) +
  geom_tile(data = mean.plot.df, aes(x = x, y = y, fill = intensity)) +
  scale_fill_gradient2()+ theme_blank()+theme(axis.text.x=element_blank(), 
                                                 axis.ticks.x=element_blank(), 
                                                 axis.text.y=element_blank(), 
                                                 axis.ticks.y=element_blank())+
  facet_grid(vars(method), vars(z))+ggtitle("Bias")
p1
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/bias1brain.jpg", height = 5, width=8)

############# Figure 2

plot.df.list.2 <- lapply(all.data.ratio.2, plot.brain.dfcreate)
mean.plot.df.2 <- rbind(data.frame(plot.df.list.2[[1]][[2]], method="mean"),
                      data.frame(plot.df.list.2[[2]][[2]], method="simple"),
                      data.frame(plot.df.list.2[[3]][[2]], method="local"),
                      data.frame(plot.df.list.2[[4]][[2]], method="hier"))
mean.plot.df.2$method <- factor(mean.plot.df.2$method, levels=c("mean",
                                                            "simple",
                                                            "local",
                                                            "hier"))
p2 <- ggplot() +
  geom_tile(data = bg.plot.df, aes(x = x, y = y), fill = "grey90",alpha = 0.8) +
  geom_tile(data = mean.plot.df.2, aes(x = x, y = y, fill = intensity)) +
  scale_fill_viridis_c()+ theme_blank()+theme(axis.text.x=element_blank(), 
                                              axis.ticks.x=element_blank(), 
                                              axis.text.y=element_blank(), 
                                              axis.ticks.y=element_blank())+
  facet_grid(vars(method), vars(z))+ggtitle("Accuracy of standard error estimation")
p2
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f2brain.jpg", height = 5, width=8)


############### Figure 3

plot.df.list.3 <- lapply(all.data.ratio.3, plot.brain.dfcreate)
mean.plot.df.3 <- rbind(data.frame(plot.df.list.3[[1]][[2]], method="mean"),
                        data.frame(plot.df.list.3[[2]][[2]], method="simple"),
                        data.frame(plot.df.list.3[[3]][[2]], method="local"),
                        data.frame(plot.df.list.3[[4]][[2]], method="hier"))
mean.plot.df.3$method <- factor(mean.plot.df.3$method, levels=c("mean",
                                                              "simple",
                                                              "local",
                                                              "hier"))
p3 <- ggplot() +
  geom_tile(data = bg.plot.df, aes(x = x, y = y), fill = "grey90",alpha = 0.8) +
  geom_tile(data = mean.plot.df.3, aes(x = x, y = y, fill = intensity)) +
  scale_fill_viridis_c()+ theme_blank()+theme(axis.text.x=element_blank(), 
                                              axis.ticks.x=element_blank(), 
                                              axis.text.y=element_blank(), 
                                              axis.ticks.y=element_blank())+
  facet_grid(vars(method), vars(z))+ggtitle("relative efficiency")
p3
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f3brain.jpg", height = 5, width=8)

############### Figure 4

plot.df.list.4 <- lapply(all.data.coverage.perc, plot.brain.dfcreate)
mean.plot.df.4 <- rbind(data.frame(plot.df.list.4[[1]][[2]], method="mean"),
                        data.frame(plot.df.list.4[[2]][[2]], method="simple"),
                        data.frame(plot.df.list.4[[3]][[2]], method="local"),
                        data.frame(plot.df.list.4[[4]][[2]], method="hier"))
mean.plot.df.4$method <- factor(mean.plot.df.4$method, levels=c("mean",
                                                                "simple",
                                                                "local",
                                                                "hier"))
p4 <- ggplot() +
  geom_tile(data = bg.plot.df, aes(x = x, y = y), fill = "grey90",alpha = 0.8) +
  geom_tile(data = mean.plot.df.4, aes(x = x, y = y, fill = intensity)) +
  scale_fill_viridis_c()+ theme_blank()+theme(axis.text.x=element_blank(), 
                                              axis.ticks.x=element_blank(), 
                                              axis.text.y=element_blank(), 
                                              axis.ticks.y=element_blank())+
  facet_grid(vars(method), vars(z))+ggtitle("Coverage of 95% CI")
p4
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f4brain.jpg", height = 5, width=8)


# Regular plots

## figure 1
brain.df.all.bias <- data.frame(bias.mean= unlist(all.data.bias.mean),
                                bias.se = unlist(all.data.se.mean),
                                methods =rep(c("mean",
                                                "simple",
                                                "local",
                                                "hier"), each=21),
                                regions=rep(row_wm$names, 4))
p1_old <- ggplot(brain.df.all.bias)+
  geom_errorbar(aes(xmin=bias.mean-1.96*bias.se, xmax=bias.mean+1.96*bias.se, y=methods), width=.2, position=position_dodge(.9)) +
  geom_point(aes(x=bias.mean, y=methods), alpha=0.7)+
  geom_vline(xintercept=0, linetype=3, alpha=0.6)+theme_minimal()+ggtitle(paste0("Mean Bias in 500 reps"))+xlab("Bias")+
  facet_wrap(~regions, nrow=5)
p1_old
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f1regular.jpg", height = 5, width=10)
######### figure 2
brain.df.all.bias <- data.frame(sd.ratio = unlist(all.data.ratio.2),
                                methods = rep(c("mean",
                                                 "simple",
                                                 "local",
                                                 "hier"), each=21),
                                regions=rep(row_wm$names, 4))
p2_old <- ggplot(brain.df.all.bias)+
  geom_bar(stat="identity",aes(sd.ratio, methods), alpha=0.7, position=position_dodge(), width=0.5)+
  geom_vline(xintercept = 1, linetype=3, alpha=0.6)+
  theme_minimal()+ggtitle(paste0("Standard error estimation accuracy"))+
  facet_wrap(~regions, nrow=5)
p2_old
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f2regular.jpg", height = 5, width=10)
######## figure 3
brain.df.ratio <- data.frame(sd.ratio = unlist(all.data.ratio.3),
                                methods = rep(c("mean",
                                                "simple",
                                                "local",
                                                "hier"), each=21),
                                regions=rep(row_wm$names, 4))
p3_old <- ggplot(brain.df.ratio)+
  geom_bar(stat="identity",aes(sd.ratio, methods), alpha=0.7, position=position_dodge(), width=0.5)+
  geom_vline(xintercept = 1, linetype=3, alpha=0.6)+
  theme_minimal()+ggtitle("Relative MSE (v.s. full data estimator)")+
  facet_wrap(~regions, nrow=5)
p3_old
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f3regular.jpg", height = 5, width=10)

######## figure 4
brain.df.ratio <- data.frame(sd.ratio = unlist(all.data.coverage.perc),
                             methods = rep(c("mean",
                                             "simple",
                                             "local",
                                             "hier"), each=21),
                             regions=rep(row_wm$names, 4))
p4_old <- ggplot(brain.df.ratio)+
  geom_bar(stat="identity",aes(sd.ratio, methods), alpha=0.7, position=position_dodge(), width=0.5)+
  geom_vline(xintercept = 0.95, linetype=3, alpha=0.6)+
  theme_minimal()+ggtitle("Coverage of 95% CI")+
  facet_wrap(~regions, nrow=5)
p4_old
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f4regular.jpg", height = 5, width=10)


