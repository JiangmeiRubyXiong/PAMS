
sample.size=200

marker="CD138"
  home.dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker,"_nophe/")
  load("~/ImageImputation/nzNormedMarkers.RData")
  #create truth
  true.quantile.df <- as.data.frame(dataSchurch %>%
                                      group_by(File.Name, groups, patients) %>%
                                      reframe(quantile = scales::percent(c(0.9, 0.95, 0.99)),
                                              value = quantile(nzNorm_CD38, c(0.9, 0.95, 0.99))))
  true.quantile.param <- data.frame(quantile=c("90%", "95%", "99%"),
                                    means=NA, sd=NA, stats=NA,
                                    methods="True")
  for(i in 1:3){
    quantiles=c("90%", "95%", "99%")[i]
    m0 <- gee(value ~ groups, id=patients, family='gaussian', corstr = "exchangeable",
              data=subset(true.quantile.df, quantile==quantiles))
    true.quantile.param[i,"means"] <- m0$coefficients[2]
    true.quantile.param[i,"sd"] <- summary(m0)$coefficients[2,4]
    true.quantile.param[i,"stats"] <- qnorm(0.975)
  }
  
  true.quantile.param$cinf.upper <- true.quantile.param$means + true.quantile.param$sd * true.quantile.param$stats
  true.quantile.param$cinf.lower <- true.quantile.param$means - true.quantile.param$sd * true.quantile.param$stats
  true.quantile.param$int.length <- true.quantile.param$sd * true.quantile.param$stats*2
  
  save(true.quantile.param, file=paste0(home.dir,"true.df.plot.geediff.RData"))
  
  cl <- 25
  registerDoParallel(cl)
  foreach(rep.i = 1:500) %dopar% {
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
    Dir <- paste0(home.dir, "rep", rep.i, "/")
    
    dir.create(Dir)
    save(patients.sampled, file=paste0(Dir, "samplePatients.RData"))
    imputation(marker=marker, trainSlides=trainSlide, testSlides=testSlide,
               dataset=subData, DIR=Dir, phenotype = FALSE)
  }
  #stop cluster
  stopCluster(cl)

markers1 <- c("CD20","CD38","CD138","CD4", "CD8", "CD68")
markers2 <- c("CD20","CD38","CD138","CD4", "CD8", "CD68")
seq.list <- as.list(seq (1,500,1)) 

directory_pathaaa <- "~/ImageImputation/metricsAllMarker/metrics_CD4_nophe"

# Get the folder names (excluding the parent directory)
folder_namesaaa <- list.dirs(directory_pathaaa, full.names = FALSE, recursive = FALSE)
folder_namesaaa <- folder_namesaaa[folder_namesaaa != "."]  # Remove the parent directory itself

# Extract numbers from folder names like 'rep123'
numbers <- as.numeric(sub("rep", "", folder_namesaaa))
setdiff(1:500, numbers)


seq.list <- as.list(c(438, 402, 419, 437, 434, 420, 435, 404, 436, 433, setdiff(1:500, numbers)))

for(marker.i in markers2){
  imputation_truth(marker=marker.i, phenot=FALSE)
  imputation_truth(marker=marker.i, phenot=TRUE)
}

for(marker.i in markers2){
  imputation_truth(marker=marker.i, phenot=FALSE)
  home.dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker.i,"_nophe/")
  mclapply(X=seq.list, FUN=imputation_boostrap, home.dir=home.dir, sample.size = 200,marker=marker.i, phenot=FALSE, mc.cores=10)
}

# imputation_boostrap(seed.path = 1, home.dir = paste0("~/ImageImputation/metricsAllMarker/metrics_","CD138","_phenotype/"), marker="CD138", phenot=TRUE)
for(marker.i in markers1){
  imputation_truth(marker=marker.i, phenot=TRUE)
  home.dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker.i,"_phenotype/")
  mclapply(X=seq.list, FUN=imputation_boostrap, home.dir=home.dir, sample.size = 200,marker=marker.i, phenot=TRUE, mc.cores=10)
}

for(marker.i in markers1){
  imputation_truth(marker=marker.i, phenot=TRUE)
  home.dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker.i,"_phenotype/")
  mclapply(X=seq.list, FUN=imputation_boostrap, home.dir=home.dir, sample.size = 200,marker=marker.i,
           phenot=TRUE, mc.cores=10, mc.preschedule = FALSE)
  
  imputation_truth(marker=marker.i, phenot=FALSE)
  home.dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker.i,"_nophe/")
  mclapply(X=seq.list, FUN=imputation_boostrap, home.dir=home.dir, sample.size = 200,marker=marker.i, 
           phenot=FALSE, mc.cores=10, mc.preschedule = FALSE)
}

# no need to run, now added to the function gee_metric
# for(marker in markers2){
#   for(rep.i in 1:500){
#     Dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker,"_phenotype/rep", rep.i, "/")
#     load(paste0(Dir, "metric.df.RData"))
#     load(paste0(Dir, "df.plot.geediff.RData"))
#     metric.df$sd <- df.plot.geediff$sd[-(1:3)]
#     save(metric.df, file=paste0(Dir, "metric.df.RData"))
#     
#     Dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker,"_nophe/rep", rep.i, "/")
#     load(paste0(Dir, "metric.df.RData"))
#     load(paste0(Dir, "df.plot.geediff.RData"))
#     metric.df$sd <- df.plot.geediff$sd[-(1:3)]
#     save(metric.df, file=paste0(Dir, "metric.df.RData"))
#     }
# }

start.time <- Sys.time()
imputation(marker=marker, trainSlides=trainSlide, testSlides=testSlide,
           dataset=subData, DIR=Dir, phenotype = phenot)
Sys.time() - start.time


