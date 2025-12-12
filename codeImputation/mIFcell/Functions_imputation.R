createFoldData <- function(dat=dataSchurch, slideTrain=slideTrain, marker="MUC5AC"){
  # normalized expression name in the dataset
  nzmarkername <- paste("nzNorm", marker, sep="_")
  #define predictor and response variables in testing set. Use the common markers for all 3 batches as the prediction variables
  commons_regress <- nzNormedMarkers[nzNormedMarkers!=nzmarkername]
  HTA_train1 <- as.matrix(dat[dat$File.Name%in%slideTrain,commons_regress])
  # Need to code tissue type and tumor (binary) as discrete numbers to run XGBoost
  HTA_train2 <- cbind(HTA_train1, dat[dat$File.Name%in%slideTrain,"groups"])
  HTA_train_marker <- dat[dat$File.Name%in%slideTrain, nzmarkername]
  return(list(HTA_train1, HTA_train2, HTA_train_marker, dat_original = dat[dat$File.Name%in%slideTrain, ]))
}

imputation <- function(marker="CD20", trainSlides=trainSlide, testSlides=testSlide,
  dataset=dataSchurch, n.imp=10, DIR="~/ImageImputation/", phenotype=TRUE){
  trainData <- createFoldData(dat=dataset,
  slideTrain = trainSlides,
  marker = marker)
  cv_model2 <- glmnet:::cv.glmnet(
    ifelse(phenotype, trainData[[2]], trainData[[1]]), trainData[[3]], alpha = 1, nfolds = 5)
  
  best_lambda2 <- cv_model2$lambda.min
  best_model2 <- glmnet::glmnet(
    ifelse(phenotype,trainData[[2]], trainData[[1]]), trainData[[3]], alpha = 1, lambda = best_lambda2)
 
  testData <- createFoldData(dat=dataset,
                    slideTrain = testSlides,
                    marker = marker)

  Yhat.test <- as.numeric(
    predict(best_model2, newx = if (phenotype) testData[[2]] else testData[[1]], s = best_lambda2)
  )

  # slide-specific calibration
  calibY <- list()
  calibYhat <- list()
  calibSlides <- list()

  # split train slide into 10 folds
  modFold <- length(trainSlides)%%10
  if(modFold!=0){
    folds <- (length(trainSlides)-modFold)/10
    split.vec <- rep(1:10, each=folds)
    split.vec <- c(split.vec, 1:(length(trainSlides)-length(split.vec)))
  } else {
    folds <- round(length(trainSlides)/10)
    split.vec <- rep(1:10, each=folds)
  }

  ### add the rest rounded into split vector
  split.trains <- split(trainSlides, split.vec)

  for(i in 1:10){
    fold.slides <- split.trains[[i]]
    trainslidei <- setdiff(trainSlides,fold.slides)
    trainDatai <- createFoldData(dat=dataset,
    slideTrain = trainslidei,
    marker = marker)
    cv_model2 <- glmnet:::cv.glmnet(ifesle (phenotype, trainDatai[[2]] ,trainDatai[[1]]), trainDatai[[3]], alpha = 1, nfolds = 5)
    
    # find optimal lambda value that minimizes test MSE
    best_lambda2 <- cv_model2$lambda.min
    # best_model2 <- glmnet::glmnet(trainDataxi, trainDatai[[3]], alpha = 1, lambda = best_lambda2)
    best_model2 <- glmnet::glmnet(ifelse (phenotype, trainDatai[[2]] ,trainDatai[[1]]), trainDatai[[3]], alpha = 1, lambda = best_lambda2)
    
    caliData <- createFoldData(dat=dataset,
      slideTrain = fold.slides,
      marker = marker)
    Yhat.i <- as.numeric(
      predict(best_model2, newx = ifelse(phenotype, caliData[[2]] ,caliData[[1]]), s = best_lambda2)
    )
    
    
    rm(trainDatai, caliData, Yhat.i, best_lambda2, best_model2)
    gc()
  }
  ################################
  ########## Use imimp ###########
  ################################

  calibY <- unlist(calibY)
  calibYhat <- unlist(calibYhat)
  calibSlides <- unlist(calibSlides)

  imp.simple <- imimp(Yhat.test, calibYhat , calibY, cdeType = "simple", x.mar = 50, y.mar=100, n.mi=n.imp)
  imp.local<- imimp(Yhat.test, calibYhat, calibY, cdeType = "local", n.mi=n.imp)

  if(phenotype){
    # Subset `dataset` for each group only once
    SlideGroup_CLR <- subset(dataset, groups == 1, select = File.Name)
    SlideGroup_DII <- subset(dataset, groups == 2, select = File.Name)
    
    # Convert to vectors for faster lookups
    SlideGroup_CLR <- SlideGroup_CLR$File.Name
    SlideGroup_DII <- SlideGroup_DII$File.Name
    
    # Find indexes directly using `SlideGroup_CLR` and `SlideGroup_DII`
    calibSlides_CLR <- which(calibSlides %in% SlideGroup_CLR)
    calibSlides_DII <- which(calibSlides %in% SlideGroup_DII)
    
    testSlidesAll <- dataset$File.Name[dataset$File.Name %in% testSlides]
    testSlides_CLR <- which(testSlidesAll %in% SlideGroup_CLR)
    testSlides_DII <- which(testSlidesAll %in% SlideGroup_DII)
    
    imp.bySample <- matrix(NA, nrow=length(testSlidesAll), ncol=n.imp)
    imp.bySample[testSlides_CLR,] <- imimp(Yhat.test[testSlides_CLR], calibYhat[calibSlides_CLR],
                                           calibY[calibSlides_CLR], cdeType = "bySample",
                                           x.mar = 50, y.mar=100,
                                           sampleID=list(calibSlides[calibSlides_CLR], testSlidesAll[testSlides_CLR]), n.mi=n.imp)
    imp.bySample[testSlides_DII,] <-imimp(Yhat.test[testSlides_DII], calibYhat[calibSlides_DII],
                                          calibY[calibSlides_DII], cdeType = "bySample",
                                          x.mar = 50, y.mar=100,
                                          sampleID=list(calibSlides[calibSlides_DII], testSlidesAll[testSlides_DII]), n.mi=n.imp)
  }else{
    imp.bySample<-imimp(Yhat.test, calibYhat , calibY, cdeType = "bySample",
                        x.mar = 50, y.mar=100, sampleID=list(as.character(calibSlides), as.character(dataset$File.Name[dataset$File.Name %in% testSlides])))
  }

  imp.list <- list(Yhat.test, imp.simple, imp.local, imp.bySample)
  save(imp.list, file=paste0(DIR,paste(substitute(dataset)),"fit", marker, ".RData"))
}


summary_quantile <- function(imp.result=imp.list.strat, marker=marker, dataset=dataSchurch, testSlides=testSlide, DIR="~/ImageImputation/"){
  yhat.test <- imp.result[[1]]
  imp.result=imp.result[2:4]
  slides.vec <- dataset$File.Name[dataset$File.Name %in% testSlides]
  quantile.list.strat <- lapply(imp.result, function(imp.mat){
    imp.df <- data.frame(imp.mat, slides=slides.vec)
    imp.df$groups <- dataset$groups[which(dataset$File.Name %in% testSlides)]
    imp.vec <- melt(as.data.table(imp.df), id.vars = c("slides","groups"), variable.name = "im")
    quantile.df <- as.data.frame(imp.vec %>%
    group_by(slides, im, groups) %>%
    reframe(quantile = scales::percent(c(0.9, 0.95, 0.99)),
    value = quantile(value, c(0.9, 0.95, 0.99))))
    return(quantile.df)
  })
  quantile.list.df.strat <- do.call(rbind, quantile.list.strat)
  quantile.list.df.strat$methods <- rep(c("simple", "local", "bySample"), each=(nrow(quantile.list.df.strat)/3))
  quantile.imp.summary.strat  <- as.data.frame(quantile.list.df.strat  %>% group_by(quantile, im, groups, methods) %>% reframe(means=mean(value)))
  # test data true
  real.quantiles.test <- single_data_process(dataset[which(dataset$File.Name %in% testSlides),c(paste0("nzNorm_", marker), "File.Name", "groups")])
  # training data true
  real.quantiles.train <- single_data_process(dataset[which(!(dataset$File.Name %in% testSlides)),c(paste0("nzNorm_", marker), "File.Name", "groups")])
  # yhat
  quantiles.yhat <- single_data_process(data.frame(yhat.test, dataset[which(dataset$File.Name %in% testSlides),c("File.Name", "groups")]))
  ll <- list(quantile.list.df.strat, real.quantiles.test, real.quantiles.train, quantiles.yhat)
  save(ll, file = paste0(DIR,"dataSchurch", marker,"quantileIMP.RData"))
}

single_data_process <- function(data.subset.test){
  names(data.subset.test) <- c("value", "slides", "groups")
  real.quantiles.test <- data.subset.test %>%
    group_by(slides) %>%
    reframe(groups=first(groups),
    quantile = scales::percent(c(0.9, 0.95, 0.99)),
    real_value = quantile(value, c(0.9, 0.95, 0.99)))
  return(real.quantiles.test)
}

var.calc <- function(qunatiledf){
  # within imputation: formula 9.2
  n.slides <- unique(qunatiledf$slides)
  within.var <- as.data.frame(qunatiledf %>% group_by(methods, quantile, im)%>%dplyr::summarize(SEi = se.sq(value), .groups = "keep"))
  within.var.pool <- as.data.frame(within.var%>%group_by(methods, quantile)%>%dplyr::summarize(VARW = mean(SEi), .groups = "keep"))
  # between imputation: formula 9.3
  between.var <- as.data.frame(qunatiledf %>% group_by(methods, quantile, im)%>%dplyr::reframe(meanImp = mean(value), .groups = "keep"))
  between.var.pool <- as.data.frame(between.var%>%group_by(methods, quantile)%>%dplyr::summarize(VARB = var(meanImp), means=mean(meanImp), .groups = "keep"))
  # pool variance: formula 9.4
  library(dplyr)
  all.vars <- left_join(within.var.pool, between.var.pool, by=c("methods", "quantile"))
  all.vars$VARB <- all.vars$VARB
  all.vars$vars <- all.vars$VARW+all.vars$VARB*((n.imp+1)/n.imp)
  # calculate df
  source("~/ImageImputation/dfCalc.R")
  all.vars.df <- dfCalc(all.vars$VARB, all.vars$vars, n.imp ,length(testSlide), 1)
  all.vars$tval95 <- qt(0.975, all.vars.df)
  return(all.vars)
}

gee_metric <- function(result.list, SlidePatient=SlidePatients, dataset=dataSchurch,DIR="~/ImageImputation/", tru.dir){
  quantile.list.df.strat <- result.list[[1]]
  trueQuantiletest <- result.list[[2]]
  trueQuantiletrain <- result.list[[3]]
  QuantileMeanImp <- result.list[[4]]
  quantile3 <- unique(quantile.list.df.strat$quantile)
  # imputation mean
  impMeans <- as.data.frame(quantile.list.df.strat %>% group_by(quantile, slides, methods) %>% reframe(means=mean(value)))
  SlidePatientGroup <- dataset %>% group_by(patients, File.Name, groups) %>% reframe()
  impMeans <- left_join(impMeans, SlidePatientGroup, by=join_by(slides==File.Name))
  trueQuantiletest <- left_join(trueQuantiletest, SlidePatient, by=join_by(slides==File.Name))
  trueQuantiletrain <- left_join(trueQuantiletrain, SlidePatient, by=join_by(slides==File.Name))
  QuantileMeanImp <- left_join(QuantileMeanImp,  SlidePatient,by=join_by(slides==File.Name))
  df.plot.geediff <- data.frame(quantile=rep(c("90%", "95%", "99%"), 5),
                                means=NA, sd=NA, stats=NA,
                                methods=rep(c("true","mean","simple","local","bySample"), each=3))
  # pool parameter estimates with Rubin's:
  # mean: mean parameter estimate
  # VarB: the variance of parameter estimate mean, of each imputation
  # VarW: the mean of parameter estimate variance, of each imputation
  quantile.list.df.strat <- left_join(quantile.list.df.strat, SlidePatient, by=join_by(slides==File.Name))
  trueData <- rbind(trueQuantiletest, trueQuantiletrain)
  meanData <- rbind(QuantileMeanImp, trueQuantiletrain)
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
      df.plot.geediff$means[idx] <- mean(m0.param$mean)
      VARB <- var(m0.param$mean)
      vars <- sqrt(var(m0.param$mean)*(n.imp+1)/n.imp+mean((m0.param$SE)^2))
      df.plot.geediff$sd[idx] <- vars
      df.plot.geediff$stats[idx] <- qt(0.975, dfCalc(VARB, vars, 10 , 200, 2))
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
}

imputation_truth <- function(marker="CD138", phenot=TRUE){
  if(phenot){
    home.dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker,"_phenotype/")
  } else {home.dir <- paste0("~/ImageImputation/metricsAllMarker/metrics_",marker,"_nophe/")}
  dir.create(home.dir)
  load("~/ImageImputation/nzNormedMarkers.RData")
  #create truth
  dataSchurch$targetMarker <- dataSchurch[,paste0("nzNorm_", marker)]
  true.quantile.df <- as.data.frame(dataSchurch %>%
                                      group_by(File.Name, groups, patients) %>%
                                      reframe(quantile = scales::percent(c(0.9, 0.95, 0.99)),
                                              value = quantile(targetMarker, c(0.9, 0.95, 0.99))))
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
}

imputation_boostrap <- function(seed.path, sample.size=30, home.dir, marker="CD138", phenot=TRUE){
  set.seed(seed.path)
  
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
  Dir <- paste0(home.dir, "rep", seed.path, "/")
  
  dir.create(Dir)
  save(patients.sampled, file=paste0(Dir, "samplePatients.RData"))
  imputation(marker=marker, trainSlides=trainSlide, testSlides=testSlide,
             dataset=subData, DIR=Dir, phenotype = phenot)
}

metrics_result_markers <- function(markers, root_dir="~/ImageImputation/metricsAllMarker"){
  full.metric <- btsp.metric <- btsp.true <- true.param <- list()
  for(marker in markers){
    ## put metrics together
    metric.list <- metric.list.nophe <- list()
    
    true.res <- get(load(paste0(root_dir, "/metrics_", marker, "_nophe/true.df.plot.geediff.RData")))
    for(rep.i in 1:500){
      Dir <- paste0(root_dir, "/metrics_",marker,"_phenotype/rep", rep.i, "/")
      load(paste0(Dir, "metric.df.RData"))
      metric.df$rep <- rep.i
      metric.list[[rep.i]] <- metric.df
      
      Dir <- paste0(root_dir, "/metrics_",marker,"_nophe/rep", rep.i, "/")
      load(paste0(Dir, "metric.df.RData"))
      metric.df$rep <- rep.i
      metric.list.nophe[[rep.i]] <- metric.df
    }
    metric.df.all <- do.call(rbind,metric.list)
    rm(list = "metric.list")
    metric.df.all$methods <- factor(metric.df.all$methods, levels=rev(c("mean", "simple","local", "bySample")))
    metric.df.all$phenotype <- "TRUE"
    metric.df.all$truth.mean <- true.res$means
    
    metric.df.nophe.all <- do.call(rbind,metric.list.nophe)
    rm(list = "metric.list.nophe")
    metric.df.nophe.all$methods <- factor(metric.df.all$methods, levels=rev(c("mean", "simple","local", "bySample")))
    metric.df.nophe.all$truth.mean <- true.res$means
    metric.df.nophe.all$phenotype <- "FALSE"
    
    # assemble true data
    true.res$marker <- marker
    true.param[[marker]] <- true.res
    
    metric.df.all <- rbind(metric.df.all, metric.df.nophe.all)
    rm(list = "metric.df.nophe.all")
    metric.df.all$marker <- marker
    full.metric[[marker]] <- metric.df.all
  }

  
  # plot 1: bias
  metric.df.all <- do.call(rbind, full.metric)
  
  metric.df.all.bias <-
    metric.df.all %>% group_by(phenotype, methods, quantile, marker) %>%
    reframe(mean.est = mean(means), se = sd(means)*(sqrt(499))/500, sd.ratio=mean(sd)/sd(means),
            btsp.bias = sqrt(mean((btspmean-truth.mean)^2)), mse.bias=sqrt(mean(bias^2)),
            coverage=mean(coverage))
  
  colnames(metric.df.all.bias)[5] <- "means"
  
  metric.df.all.bias$methods <- factor(metric.df.all.bias$methods, 
                                       levels=c("bySample", "local", "simple", "mean"),
                                       labels = c("Local\nMultiple", "Local",
                                                  "Simple", "Synthesis"))
  metric.df.all.bias$phenotype <- as.logical(metric.df.all.bias$phenotype)
  metric.df.all.bias$phenotype <- ifelse( metric.df.all.bias$phenotype, "With", "Without")
  p1.vars <- c("quantile", "means", "se", "marker", "methods", "phenotype")
  metric.df.all.bias.p1 <- metric.df.all.bias[,p1.vars]
  metric.df.all.bias.p1$type <- metric.df.all.bias$phenotype
  true.param.all <- do.call(rbind, true.param)
  custom_colors <- c("With" = "#5D7A96",
                     "Without" = "#B87D7D")
# Clean version 
p1 <- ggplot(metric.df.all.bias.p1, aes(x = methods, y = means, color = phenotype)) +
      geom_point(alpha = 0.7, position = position_dodge(0.9)) +
      geom_hline(data=true.param.all, aes(yintercept=means, linetype="True Parameter"), color="forestgreen")+
      scale_linetype_manual(values = c("True Parameter" = "dashed"), name=NULL) + labs(linetype = "True Parameter") +
      facet_wrap(quantile ~ marker, scales = "free_y", ncol=6) +
      scale_color_manual(values = custom_colors) +
      ggtitle("Expected Value of Group Difference (DII v.s. CLR)") +
      ylab("Difference in Mean 99th Quantile")+
      labs(color="Phenotype", x="Methods")+
      theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  ### plot 2: MSE of imputation estimator / SD of boostrap estimator
  custom_colors2 <- c("With" = "#5D7A96",
                   "Without" = "#B87D7D")  

  p2 <- ggplot(metric.df.all.bias)+
    geom_bar(stat="identity",aes(sd.ratio, methods, fill=phenotype), alpha=0.7, position=position_dodge(), width=0.5)+
    facet_grid(rows=vars(quantile), cols=vars(marker), scales = "free")+
    geom_vline(xintercept = 1, linetype="dashed", alpha=0.6, color="forestgreen")+
    ggtitle(paste0("Standard Error Estimation Accuracy"))+
    scale_fill_manual(values = custom_colors2)+coord_flip()+theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(fill="Phenotype", y="Methods", x="Ratio of Standard Deviation")
  
  metric.df.all.bias$mse.ratio <- metric.df.all.bias$mse.bias/metric.df.all.bias$btsp.bias
  p3 <- ggplot(metric.df.all.bias)+
    geom_bar(stat="identity",aes(mse.ratio, methods, fill=phenotype), alpha=0.7, position=position_dodge(), width=0.5)+
    facet_grid(rows=vars(quantile), cols=vars(marker), scales = "free")+
    geom_vline(xintercept = 1, linetype=3, alpha=0.6)+
    ggtitle("Relative RMSE (v.s. Full Data Estimator)")+
    scale_fill_manual(values = custom_colors2)+coord_flip()+theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    labs(fill="Phenotype", y="Methods", x="Ratio of RMSE")
  
  # plot 4: coverage
  options(na.rm=TRUE)
  library(tidyr)
  table1 <- spread(metric.df.all.bias[,c("methods", "coverage", "marker", "quantile", "phenotype")], methods, coverage, marker)
  for(i in 4:7){
    table1[[i]] <- as.numeric(table1[[i]])
  }
  table2 <- reshape2::melt(
    data.frame((table1[,4:7]), quantile=table1$quantile, phenotype=table1$phenotype, marker=table1$marker),
                 id.vars=c("quantile", "phenotype", "marker"), variable.name = "methods")
  
   table2$methods <- factor(table2$methods,
                            labels = c("Local\nMultiple", "Local","Simple", "Synthesis"))
  p4 <- ggplot(table2)+
    geom_bar(stat="identity",aes(methods, value, fill=phenotype), alpha=0.7, position=position_dodge(), width=0.5)+
    facet_grid(rows=vars(quantile), cols=vars(marker), scales = "free")+
    xlab("Percentage")+ geom_vline(xintercept = 0, linetype=2, alpha=0.8)+
    geom_hline(aes(yintercept = 0.95, linetype = "95% CI Target"), 
               color = "forestgreen") +
    scale_linetype_manual(values = c("95% CI Target" = "dashed"), name = NULL) +
    labs(linetype = "Reference Line") + ggtitle("Coverage of 95% CI")+
    scale_fill_manual(values = custom_colors2)+theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
    labs(fill="Phenotype", x="Methods", y="Coverage")
  
  return(list(p1,p2,p3,p4))
}


# # linear regression predictability check
# load("~/ImageImputation/nzNormedMarkers.RData")
# dataSchurch <- read.csv("~/ImageImputation/dataSchurchLite.csv")
# dataSchurch <- dataSchurch[,nzNormedMarkers]
# summary(lm(nzNorm_CD138~., data=dataSchurch))#Multiple R-squared:  0.5154,	Adjusted R-squared:  0.5153 
# summary(lm(nzNorm_CD20~., data=dataSchurch))#Multiple R-squared:  0.5036,	Adjusted R-squared:  0.5035 
# summary(lm(nzNorm_CD38~., data=dataSchurch))#Multiple R-squared:  0.5052,	Adjusted R-squared:  0.5051 
# summary(lm(nzNorm_CD4~., data=dataSchurch))#Multiple R-squared:  0.7437,	Adjusted R-squared:  0.7437 
# summary(lm(nzNorm_CD68~., data=dataSchurch))#Multiple R-squared:  0.6029,	Adjusted R-squared:  0.6029 
# summary(lm(nzNorm_CD8~., data=dataSchurch))#Multiple R-squared:  0.4995,	Adjusted R-squared:  0.4993 
