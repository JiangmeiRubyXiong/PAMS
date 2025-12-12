

library(ggplot2)
library(sandwich)
library(lmtest)

# calibration data metrics
FA_synth_region <- readRDS("/media/disk2/beijing_dti/Neuroimaging/FA_synth_region.rds")
FA_truth_region <- readRDS("/media/disk2/beijing_dti/Neuroimaging/FA_truth_region.rds")
calibration_error <- readRDS("/media/disk2/beijing_dti/Neuroimaging/calibration_error.rds")
testt_subject_data <- read.csv("/media/disk2/beijing_dti/Neuroimaging/testt_subject_data.csv")

imputation_seed_generation <- readRDS("/media/disk2/beijing_dti/Neuroimaging/imputation_seed_generation.rds")
btsp_calib_region_df <- function(seed.i, n.imp=100, sub_data = testt_subject_data, 
                        sample.idx = imputation_seed_generation, FA_synth_regions=FA_synth_region,
                        FA_truth_regions=FA_truth_region,
                        roi.df = row_wm_100){
  set.seed(seed.i)
  btsp.idx <- sample.idx[[seed.i]]

  calib.data <- FA_truth_regions[btsp.idx]
  regions <- calib.data[[1]]$names
  
  calib.data.all <- data.frame(do.call(rbind,lapply(calib.data, "[[", "values")))
  names(calib.data.all) <- regions
  calib.data.all$Sex <- sub_data$Sex[btsp.idx]
  calib.data.all$Age <- sub_data$Age[btsp.idx]

  list.res.calib <- list.res.complete <- list()
  for(ROI in regions){
    formulai <- as.formula(paste0(ROI," ~ Sex"))
    m0 <- lm(formulai, data=calib.data.all)
    estimate <- m0$coefficients[2]
    stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
    list.res.calib[[ROI]] <- data.frame(Estimate = estimate, StdError = stdError)
    
    # complete data only
    calib.data.complete <- calib.data.all[1:250,]
    m0 <- lm(formulai, data=calib.data.complete)
    estimate <- m0$coefficients[2]
    stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
    list.res.complete[[ROI]] <- data.frame(Estimate = estimate, StdError = stdError)
  }
  list.res.calib <- do.call(rbind, list.res.calib)
  rownames(list.res.calib) <- regions
  names(list.res.complete) <- regions
  saveRDS(list.res.calib, paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_region_100/truedata", seed.i, ".rds"))
  saveRDS(list.res.complete, paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_region_100/complete", seed.i, ".rds"))
}

mclapply(1:500, btsp_calib_region_df, mc.cores=20, mc.preschedule = FALSE)

# calculate true data
true.data <- FA_truth_region
regions <- true.data[[1]]$names
true.data.all <- data.frame(do.call(rbind,lapply(true.data, "[[", "values")))
names(true.data.all) <- regions
true.data.all$Sex <- testt_subject_data$Sex
true.data.all$Age <- testt_subject_data$Age
list.res.true <- list()
for(ROI in regions){
  formulai <- as.formula(paste0(ROI," ~ Sex"))
  m0 <- lm(formulai, data=true.data.all)
  estimate <- m0$coefficients[2]
  stdError <- sqrt(vcovHC(m0, type = "HC3")[2,2])
  test.res <- coeftest(m0, vcov = vcovHC(m0, type="HC3"))
  list.res.true[[ROI]] <- data.frame(Estimate = estimate, StdError = stdError,
                                     tval=test.res[2,3], pval=test.res[2,4])
}
list.res.true <- do.call(rbind, list.res.true)
rownames(list.res.true) <- regions
list.res.true <- data.frame(list.res.true)
saveRDS(list.res.true, "/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/truedata_100.rds")
truedata_100 <- readRDS("/media/disk2/beijing_dti/Neuroimaging/imputation_reg_res/truedata_100.rds")

## calculate metrics: bias, se, etc for all methods
all.data <- lapply(1:200, function(i){readRDS(file=paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_region/simulation",i,".rds"))})
list.orders <- names(all.data[[1]])
all.data <- lapply(list.orders, function(i){lapply(all.data, "[[", i)}) 
names(all.data) <- list.orders
all.data.means <- lapply(all.data, function(lists){vector <- lapply(lists, do.call, what=rbind)
                                                   vector <- lapply(vector, function(matr){matr[,1]})
                                                   do.call(cbind, vector)})
all.data.mean <- lapply(all.data.means, rowMeans)


all.data.bias.mean <- lapply(all.data.mean, function(x, truemean=truedata_100){x-truemean[,1]})

all.data.se <- lapply(all.data, function(lists){vector <- lapply(lists, do.call, what=rbind)
                                                vector <- lapply(vector, function(matr){matr[,2]})
                                                do.call(cbind, vector)})

# simulation truth calculation
true.data.sim <- lapply(1:200, function(i){readRDS(file=paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_region/truedata",i,".rds"))})
true.data.means <- do.call(cbind,lapply(true.data.sim, function(x){x[,1]}))
true.data.mean <- rowMeans(true.data.means)

row.se <- function(x){
  vec.sd <- list()
  for(i in 1:nrow(x)){
    vec.sd[[i]] <- sd(x[i,])
  }
  return(unlist(vec.sd))
}
true.data.est.se <- row.se(true.data.means)

# Mean SE with model
all.data.se.mean <- lapply(all.data.se, rowMeans)
# SE of estimator

all.data.est.se <- lapply(all.data.means, row.se)

all.data.ratio.2 <- mapply(function(upper, lower){return(upper/lower)},
                            all.data.se.mean, all.data.est.se, SIMPLIFY = FALSE)

# MSE
## Numerator
all.data.means <- lapply(all.data, function(lists){vector <- lapply(lists, do.call, what=rbind)
                                                    vector <- lapply(vector, function(matr){matr[,1]})
                                                    do.call(cbind, vector)})
mse.bias <-  function(x, truedata=truedata_100){
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
  readRDS(paste0("/media/disk2/beijing_dti/Neuroimaging/imputation_region/truedata",x,".rds"))[,1]})
sim.res <- do.call(cbind, sim.res)
sim.res.bias <- mse.bias(sim.res)
## ratio
all.data.ratio.3 <- lapply(all.data.mse.est,function(upper, lower=sim.res.bias){return(upper/lower)})

# Coverage
all.data.means <- lapply(all.data, function(lists){vector <- lapply(lists, do.call, what=rbind)
                                                    vector <- lapply(vector, function(matr){matr[,1]})
                                                    do.call(cbind, vector)})
all.data.bias.abs <- lapply(all.data.means, function(x, truedata=truedata_100){
  abs(x - matrix(truedata[,1], nrow=nrow(truedata), ncol=ncol(x), byrow = FALSE))
})


all.data.stats <- lapply(all.data[2:3], function(lists){vector <- lapply(lists, do.call, what=rbind)
                                                    vector <- lapply(vector, function(matr){matr[,3]})
                                                    do.call(cbind, vector)})
all.data.stats <- append(list(matrix(1.96, nrow=nrow(row_wm_100), ncol=200)), all.data.stats)

all.data.coverage <- mapply(function(stats, se, abs.bias){return(stats * se >= abs.bias)}, 
                            all.data.stats, all.data.se, all.data.bias.abs, SIMPLIFY = FALSE)
all.data.coverage.perc <- lapply(all.data.coverage, rowMeans)

## figure 1
truedata_100$regions <- rownames(truedata_100)

brain.df.all.bias <- data.frame(bias.mean= c(truedata_100$Estimate,true.data.mean,unlist(all.data.mean)),
                                bias.se = (c(truedata_100$StdError,true.data.est.se,unlist(all.data.est.se)))*sqrt(499/(500^2)),
                                methods =rep(c("True","BTSP",list.orders), each=nrow(row_wm_100)),
                                regions=c(rownames(truedata_100),rownames(true.data.sim[[1]]),rep(row_wm_100$names, length(list.orders))),
                                significant = rep( factor(truedata_100$pval<0.2), (length(list.orders)+2)))
brain.df.all.bias$methods <- factor(brain.df.all.bias$methods, levels=c(rev(list.orders), "BTSP", "True"))
p1_old <- ggplot(brain.df.all.bias)+
  geom_errorbar(aes(xmin=bias.mean-1.96*bias.se, xmax=bias.mean+1.96*bias.se, y=methods, color = significant), width=.2, position=position_dodge(.9)) +
  geom_point(aes(x=bias.mean, y=methods, color = significant), alpha=0.7)+ scale_x_continuous(n.breaks = 3)+ 
  geom_vline(aes(xintercept=Estimate), data=truedata_100, linetype=3, alpha=0.6)+theme_minimal()+ggtitle(paste0("Mean in 200 reps: n=500, with calibration"))+xlab("Bias")+
  facet_wrap(~regions, nrow=3, scales = "free_x")+ labs(colour="p.val < 0.2")
p1_old

ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f1regular.jpg", height = 5, width=10)

######### figure 2
brain.df.all.bias.2 <- data.frame(sd.ratio = unlist(all.data.ratio.2),
                                methods = rep(list.orders, each=17),
                                regions=rep(row_wm_100$names, length(list.orders)),
                                significant = rep( factor(truedata_100$pval<0.2), length(list.orders)))
brain.df.all.bias.2$methods <- factor(brain.df.all.bias.2$methods, levels=rev(list.orders))
p2_old <- ggplot(brain.df.all.bias.2)+
  geom_bar(stat="identity",aes(sd.ratio, methods, fill=significant), alpha=0.7, position=position_dodge(), width=0.5)+
  geom_vline(xintercept = 1, linetype=3, alpha=0.6)+
  theme_minimal()+ggtitle(paste0("Standard error estimation accuracy"))+
  facet_wrap(~regions, nrow=5)+ labs(fill="p.val < 0.2")
p2_old
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f2regular.jpg", height = 5, width=10)
######## figure 3
brain.df.ratio <- data.frame(sd.ratio = unlist(all.data.ratio.3),
                                methods = rep(list.orders, each=21),
                                regions=rep(row_wm$names, length(list.orders)),
                             significant = rep( factor(truedata_100$`Pr(>|t|)`<0.2), length(list.orders)))
brain.df.ratio$methods <- factor(brain.df.ratio$methods, levels=rev(list.orders))
p3_old <- ggplot(brain.df.ratio)+
  geom_bar(stat="identity",aes(sd.ratio, methods, fill=significant), alpha=0.7, position=position_dodge(), width=0.5)+
  geom_vline(xintercept = 1, linetype=3, alpha=0.6)+
  theme_minimal()+ggtitle("Relative MSE (v.s. full data estimator)")+
  facet_wrap(~regions, nrow=5, scales = "free_x")+ labs(fill="p.val < 0.2")
p3_old
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f3regular.jpg", height = 5, width=10)

######## figure 4
brain.df.ratio.4 <- data.frame(sd.ratio = unlist(all.data.coverage.perc),
                             methods = rep(list.orders, each=17),
                             regions=rep(row_wm_100$names, length(list.orders)),
                             significant = rep( factor(truedata_100$`pval`<0.2), length(list.orders)))
brain.df.ratio.4$methods <- factor(brain.df.ratio.4$methods, levels=rev(list.orders))
p4_old <- ggplot(brain.df.ratio.4)+
  geom_bar(stat="identity",aes(sd.ratio, methods, fill=significant), alpha=0.7, position=position_dodge(), width=0.5)+
  geom_vline(xintercept = 0.95, linetype=3, alpha=0.6)+
  theme_minimal()+ggtitle("Coverage of 95% CI")+
  facet_wrap(~regions, nrow=3)+labs(fill="p.val < 0.2")
p4_old
ggsave(filename = "/media/disk2/beijing_dti/Neuroimaging/f4regular.jpg", height = 5, width=10)

## checking correlation: synthesized data & calibration error
regions <- 1:17
FA_synth_patient <- list()
for(i in regions){
  FA_synth_patient[[i]] <- unlist(lapply(FA_synth_region, function(df, idx=i){df$value[idx]}))
}
FA_truth_patient <- list()
for(i in regions){
  FA_truth_patient[[i]] <- unlist(lapply(FA_truth_region, function(df, idx=i){df$value[idx]}))
}
calib_error_patient <- list()
for(i in regions){
  calib_error_patient[[i]] <- unlist(lapply(calibration_error, function(df, idx=i){df[idx]}))
}

# corr_check <- mapply(function(x,y){cor(x$value, y)}, x=FA_synth_region, y=calibration_error)
# hist(corr_check, main = "Correlation between Synthesized and Calibration Error")
# 
# corr_check2 <- mapply(function(x,y){cor(x$value, y)}, x=FA_truth_region, y=calibration_error)
# hist(corr_check2, main = "Correlation between True Data and Calibration Error")

corr_check_patient <- mapply(function(x,y){cor(x, y)}, x=FA_synth_patient, y=calib_error_patient)
corr_check2_patient <- mapply(function(x,y){cor(x, y)}, x=FA_truth_patient, y=calib_error_patient)
data.frame(syth_corr = corr_check_patient, truth_corr = corr_check2_patient, row.names = unlist(row_wm_100$names))

# plot: Y - Yhat (y axis) against Yhat (x axis)
calib_error_patient_df <- as.data.frame(calib_error_patient, col.names = unlist(row_wm_100$names))
calib_error_patient_df$ID <- 1:nrow(calib_error_patient_df)
FA_synth_patient_df<- as.data.frame(FA_synth_patient, col.names = unlist(row_wm_100$names))
FA_synth_patient_df$ID <- 1:nrow(FA_synth_patient_df)
FA_truth_patient_df<- as.data.frame(FA_truth_patient, col.names = unlist(row_wm_100$names))
FA_truth_patient_df$ID <- 1:nrow(FA_truth_patient_df)

calib_error_patient_df_plot <- reshape2::melt(calib_error_patient_df, id.vars="ID")
FA_synth_patient_df_plot <- reshape2::melt(FA_synth_patient_df, id.vars="ID")
FA_truth_patient_df_plot <- reshape2::melt(FA_truth_patient_df, id.vars="ID")

plot.df.1 <- dplyr::left_join(FA_synth_patient_df_plot, calib_error_patient_df_plot, by=join_by("ID", "variable"))
library(ggplot2)
ggplot(plot.df.1)+
  geom_point(aes(x=value.x, y=value.y), alpha=0.3, color="forestgreen")+
  xlab("Yhat")+ylab("Y - Yhat")+
  ggtitle("Y - Yhat against Yhat")+
  facet_wrap(~variable, nrow=4, scales = "free")+theme_minimal()

# plot: Y (y axis) against Yhat (x axis)

plot.df.2 <- dplyr::left_join(FA_synth_patient_df_plot, FA_truth_patient_df_plot, by=join_by("ID", "variable"))
p2 <- ggplot(plot.df.2)+
  geom_point(aes(x=value.x, y=value.y), alpha=0.3, color="blue")+
  xlab("Yhat")+ylab("Y")+
  ggtitle("Y against Yhat")+
  facet_wrap(~variable, nrow=4, scales = "free")+theme_minimal()

# plot: histogram of calibration errors
ggplot(calib_error_patient_df_plot, aes(x=value, y = stat(density))) + ggtitle("Histogram of calibration errors")+xlab("Y - Yhat")+
  geom_histogram(fill="orange" ,alpha=0.7)+facet_wrap(~variable, nrow=4, scales = "free")+theme_minimal()

# patient ID 68, pixel comparison across all regions
FA_atlas_zeros= c(readNifti("/media/disk2/beijing_dti/Neuroimaging/train_T1_FA/FA_2mm_FAmask3.nii.gz"))
FA_atlas_zeros = FA_atlas_zeros[FA_atlas_zeros>0]
FA_atlas_zeros_region <- which(FA_atlas_zeros %in% unlist(row_wm_100$label))
FA_atlas_zeros_region_label <- FA_atlas_zeros[FA_atlas_zeros_region]
FA_synth_list <- readRDS(file="/media/disk2/beijing_dti/Neuroimaging/FA_synth_list.rds")
FA_truth_list <- readRDS(file="/media/disk2/beijing_dti/Neuroimaging/FA_truth_list.rds")

ID68_synth_region <- FA_synth_list[[68]][FA_atlas_zeros_region]
ID68_truth_region <- FA_truth_list[[68]][FA_atlas_zeros_region]

plot.68.df <- data.frame(Yhat=ID68_synth_region, Y=ID68_truth_region, label = FA_atlas_zeros_region_label)
plot.68.df <- dplyr::left_join(plot.68.df, row_wm_100, by="label")
ggplot(plot.68.df) + geom_point(aes(Yhat, Y), alpha=0.2, color="orchid1")+
  ggtitle("ID:138")+ geom_abline(slope = 1, intercept = 0, linetype=3)+
  facet_wrap(~names)+theme_minimal()

# correlation
regions <- unlist(row_wm_100$names)
corr.68 <- list()
for(i in regions){
  corr.68[[i]] <- cor(plot.68.df$Yhat[plot.68.df$names==i], plot.68.df$Y[plot.68.df$names==i])
}

## for patient ID 11
ID11_synth_region <- FA_synth_list[[11]][FA_atlas_zeros_region]
ID11_truth_region <- FA_truth_list[[11]][FA_atlas_zeros_region]

plot.11.df <- data.frame(Yhat=ID11_synth_region, Y=ID11_truth_region, label = FA_atlas_zeros_region_label)
plot.11.df <- dplyr::left_join(plot.11.df, row_wm_100, by="label")
ggplot(plot.11.df) + geom_point(aes(Yhat, Y), alpha=0.2, color="orchid1", stroke=0)+
  ggtitle("ID:34")+ geom_abline(slope = 1, intercept = 0, linetype=3)+
  facet_wrap(~names)+theme_minimal()

corr.11 <- list()
for(i in regions){
  corr.11[[i]] <- cor(plot.11.df$Yhat[plot.11.df$names==i], plot.11.df$Y[plot.11.df$names==i])
}

## for patient ID 88
ID88_synth_region <- FA_synth_list[[88]][FA_atlas_zeros_region]
ID88_truth_region <- FA_truth_list[[88]][FA_atlas_zeros_region]

plot.88.df <- data.frame(Yhat=ID88_synth_region, Y=ID88_truth_region, label = FA_atlas_zeros_region_label)
plot.88.df <- dplyr::left_join(plot.88.df, row_wm_100, by="label")
ggplot(plot.88.df) + geom_point(aes(Yhat, Y), alpha=0.2, color="orchid1", stroke=0)+
  ggtitle("ID:179")+ geom_abline(slope = 1, intercept = 0, linetype=3)+
  facet_wrap(~names)+theme_minimal()

corr.88 <- list()
for(i in regions){
  corr.88[[i]] <- cor(plot.88.df$Yhat[plot.88.df$names==i], plot.88.df$Y[plot.88.df$names==i])
}

## for patient ID 83
ID83_synth_region <- FA_synth_list[[83]][FA_atlas_zeros_region]
ID83_truth_region <- FA_truth_list[[83]][FA_atlas_zeros_region]

plot.83.df <- data.frame(Yhat=ID83_synth_region, Y=ID83_truth_region, label = FA_atlas_zeros_region_label)
plot.83.df <- dplyr::left_join(plot.83.df, row_wm_100, by="label")
ggplot(plot.83.df) + geom_point(aes(Yhat, Y), alpha=0.2, color="orchid1", stroke=0)+
  ggtitle("ID:169")+ geom_abline(slope = 1, intercept = 0, linetype=3)+
  facet_wrap(~names)+theme_minimal()

corr.83 <- list()
for(i in regions){
  corr.83[[i]] <- cor(plot.83.df$Yhat[plot.83.df$names==i], plot.83.df$Y[plot.83.df$names==i])
}
### region mean correlation and comparison
corr.region <- list()
for(i in regions){
  corr.region[[i]] <- cor(plot.df.2$value.x[plot.df.2$variable==i], plot.df.2$value.y[plot.df.2$variable==i])
}
compare.plot.df <- data.frame(regionMean = unlist(corr.region), patient138=unlist(corr.68), patient34=unlist(corr.11), patient179=unlist(corr.88), patient169=unlist(corr.83),region=regions)
compare.plot.df <- reshape2::melt(compare.plot.df, id.vars="region")
ggplot(compare.plot.df, aes(x=value, y=region, color=variable))+geom_point(size=3, alpha=0.6)+
  scale_color_manual(values=c("red", "grey80", "grey70","grey60","grey50"))+
  ggtitle("Compare region mean and selected patient \n Y/Yhat correlation")+theme_minimal()

