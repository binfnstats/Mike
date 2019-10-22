
#Boruta Feature selection
Boruta_feat<-function(dat,labs){
  PDBP_Bor<-Boruta(labs~.,data=dat,doTrace=0,maxRuns = 1000)
  #Bor_features<-
  as.vector(getSelectedAttributes(PDBP_Bor, withTentative = F))}

#Differential Expression Analysis 
Diffexp_top<-function(dat,nfeatin=1000,nfeatout=100){
  design <- model.matrix(~0 + dat$labs, data = dat$labs)
  fit <- lmFit(t(dat[,2:nfeatin]),design)
  fit <- eBayes(fit,trend=TRUE)
  plotSA(fit, main="Probe-level")
  as.vector(rownames(topTable(fit,adjust = "fdr",sort.by = "B",number = nfeatout)))}

# load data.
load_PDBP <- function(dat,labs,feats="all",featnum=100,exp="BF") {
  mod_data<-readRDS(dat)
  mod_data<-mod_data[mod_data$Q2600==exp,]
  pd55_sub_labels<-readRDS(labs)
  #BF_dat_labs<-readRDS("BF_dat_labs.rds")
  labs<-pd55_sub_labels
  df=cbind(mod_data,labs)
  bor=Boruta_feat(df,labs)
  diff=Diffexp_top(df,1000,featnum)
  bb_ncv_feat=bb_ncv_feat("model_ncv_allfeat.rds",0.1) #set threshold
    if (feats=="all") {features=1:featnum}
    if (feats=="bor") {features=c(bor)}
    if (feats=="diff") {features=c(diff)}
    if (feats=="bb") {features=c(bb_ncv_feat)}
 
  data_mod<-mod_data[,features]
  saveRDS(data_mod,file="data_mod_current_run.RDS")
  cbind(labs,data_mod)
}


#2. Run Caret Models

#change to df_rand for randomised samples

traintest <-function(df){
  set.seed(107)
  inTrain <- createDataPartition(y = df$diagnosis, p = .6, list = FALSE)
  df[ inTrain,]
  testing <- df[-inTrain,]  
}


bplot <-function(res_mlist){
  #transform dataset to view boxplot
  res_mlist=as.data.frame(res_mlist$values)
  res_mlist=gather(res_mlist,model,value,-Resample)
  res_mlist2=res_mlist[grep("ROC", res_mlist$model) ,]
  # Boxplot of model performance
  res_mlist2 %>%
    ggplot( aes(x=model, y=value, fill=model)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    #theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("AUC values for models run") +
    xlab("")}


mod_predict <-function(df,df_test,nih){
  set.seed(1951)  # set the seed
  my_control <- trainControl(
    method = "repeatedcv",number =3,
    ## repeated ten times
    repeats = 10,
    savePredictions="final",
    classProbs=TRUE,
    index=createResample(df$labs, 20),
    summaryFunction=twoClassSummary)
  
  model_list <- caretList(
    labs~., data=df,
    trControl=my_control,
    methodList=c("rf", "glm","svmLinear","rpart","treebag","nnet"))#,"nnet" "glmboost", "nnet", "treebag", "svmLinear", "rpart"))
  
  greedy_ensemble <- caretEnsemble(
    model_list, 
    metric="ROC",
    trControl=trainControl(
      method = "repeatedcv",number =3,
      ## repeated ten times
      repeats = 10,
      savePredictions="final",
      classProbs=TRUE,
      index=createResample(df$labs, 20),
      summaryFunction=twoClassSummary))
  
  model_preds <- lapply(model_list, predict, newdata=df_test, type="prob")
  model_preds <- data.frame(model_preds)
  ens_preds <- predict(greedy_ensemble, newdata=df_test, type="prob")
  model_preds$ensemble <- ens_preds
  
  nih_preds <- lapply(model_list, predict, newdata=nih*9.35/8.25, type="prob")
  nih_preds <- data.frame(nih_preds)
  ens_preds_nih <- predict(greedy_ensemble, newdata=nih*9.35/8.25, type="prob")
  nih_preds$ensemble <- ens_preds_nih
  
  
  saveRDS(model_preds,file="model_preds.RDS")
  saveRDS(nih_preds,file="nih_preds.RDS")
  saveRDS(model_list,file="model_list.RDS")
  res_mlist=resamples(model_list)
  saveRDS(res_mlist,file="res_mlist.RDS")
  #show boxplot
  bplot(res_mlist)
  target<-model_preds
}

#3. Run blkbox Models


#4. Assess Performance

mod_diag <-function(model_preds,test,Label_graph){
  y =ordered(test$labs,levels = c("PD","HC"),labels = c(0, 1))
  mlist_rf <- plot.roc(y, model_preds$rf.PD,
                       main=Label_graph,
                       percent=TRUE,
                       col="#1c61b6",print.auc.y = 60,print.auc = TRUE)
  
  
  mlist_glm <- plot.roc(y, model_preds$glm.PD, 
                        percent=TRUE, 
                        col="#008600",print.auc = TRUE,print.auc.y = 55,add=TRUE)
  mlist_nnet <- plot.roc(y, model_preds$nnet.PD,
                         percent=TRUE,
                         col="black",print.auc = TRUE,print.auc.y = 50,add=TRUE)
  mlist_svm <- plot.roc(y, model_preds$svmLinear.PD, 
                        percent=TRUE, 
                        col="red",print.auc = TRUE,print.auc.y = 45,add=TRUE)
  legend("bottomright", legend=c("Random Forest", "GLM",#"NNET",
                                 "SVM Linear"), col=c("#1c61b6", "#008600","black","red"), lwd=2)
  
  mlist_treebag <- plot.roc(y, model_preds$treebag.PD,
                            main="Treebag",
                            percent=TRUE,
                            col="#1c61b6",print.auc.y = 60,print.auc = TRUE)
  legend("bottomright", legend=c("Treebag"), col=c("#1c61b6"), lwd=2)
  
 
  
  mlist_ens <- plot.roc(y, model_preds$ensemble,
                        main="Ensemble model",
                        percent=TRUE,
                        col="#1c61b6",print.auc.y = 60,print.auc = TRUE)
  legend("bottomright", legend=c("Ensemble"), col=c("#1c61b6"), lwd=2)
  
  aucvals=as.data.frame(list(auc(y, model_preds$glm.PD),auc(y, model_preds$nnet.PD)
                             ,auc(y, model_preds$rf.PD),
                             auc(y, model_preds$rpart.PD),
                             auc(y, model_preds$svmLinear.PD),auc(y, model_preds$treebag.PD),
                             auc(y, model_preds$ensemble)))
 
  names(aucvals) <- c( "glm","nnet"
                       ,"rf","rpart","svm","treebag","ensemble")
  saveRDS(aucvals,file="aucvals.RDS")
  
  
  }
  


#Extreme Gradient Boosting and GBM

gbmxgb_pred <-function(training,testing){
  trainX <-training[,names(training)!="labs"]        # Pull out the dependent variable
  testX <- testing[,names(training)!="labs"]
  ctrl <- trainControl(method = "repeatedcv",   # 10fold cross validation
                       number = 5,							# do 5 repititions of cv
                       summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                       classProbs=TRUE,
                       allowParallel = TRUE)
  grid <- expand.grid(interaction.depth=c(1,2), # Depth of variable interactions
                      n.trees=c(10,20),	        # Num trees to fit
                      shrinkage=c(0.01,0.1),		# Try 2 values for learning rate 
                      n.minobsinnode = 5)
  
  set.seed(1951)  # set the seed
  gbm.tune <- train(x=trainX,y=training$labs,
                    method = "gbm",
                    metric = "ROC",
                    trControl = ctrl,
                    tuneGrid=grid,
                    verbose=FALSE)
  gbm.pred <- predict(gbm.tune,testX)
  #saveRDS(gbm.pred,"gbm_pred.rds")
  #saveRDS(gbm.tune,"gbm_tune.rds")
  
  gbm.probs <- predict(gbm.tune,testX,type="prob")
  
  
  set.seed(1951)
  
  xgbGrid <- expand.grid(nrounds = c(100,200),  # this is n_estimators in the python code above
                         max_depth = c(10, 15, 20, 25),
                         colsample_bytree = seq(0.5, 0.9, length.out = 5),
                         ## The values below are default values in the sklearn-api. 
                         eta = 0.1,
                         gamma=0,
                         min_child_weight = 1,
                         subsample = 1)
  
  xgb.tune <-train(x=trainX,y=training$labs,
                   method="xgbTree",
                   metric="ROC",
                   trControl=ctrl,
                   tuneGrid=xgbGrid)
  
  xgb.probs <- predict(xgb.tune,testX,type="prob")
  #head(xgb.probs)
  
  rValues <- resamples(list(xgb=xgb.tune,gbm=gbm.tune))
  
  GBM_model <- plot.roc(testing$labs, gbm.probs$PD,
                        main="XGB and GBM Models",
                        percent=TRUE,
                        col="#1c61b6",print.auc.y = 60,print.auc = TRUE)
  XGB_model <- plot.roc(testing$labs, xgb.probs$PD,
                        percent=TRUE,
                        col="#008600",print.auc.y = 80,print.auc = TRUE,add=TRUE)
  legend("bottomright", legend=c("GBM", "XGB"), col=c("#1c61b6", "#008600"), lwd=2)
  
  saveRDS(gbm.probs,file="gbm.RDS")
  saveRDS(xgb.probs,file="xgb.RDS")
  
  bwplot(rValues,metric="ROC",main="GBM vs xgboost")	# boxplot
  dotplot(rValues,metric="ROC",main="GBM vs xgboost")}	# dotplot}


#BLKBOX MODELS

blkbox_run<-function(data=data_mod,labels=labs){
  #labs<-as.data.frame(labels)
  pd55_sub_part=Partition(data,labels=labels,size=.6,seed=1234)
  bb <- blkbox(data = pd55_sub_part,exclude="xgboost")
  # Performance
  bb_perf = Performance(bb)
  bb_ncv<- blkboxNCV(data = data,
                     labels = labels,
                     Method = "randomforest",exclude = "xgboost",
                     AUC = 0.8,metric="AUROC")
  # Standard ROC curve of each model
  blkboxROC(bb_perf)
  #plot ncv averages for blkbox
  ncv.plot(bb_ncv, metric = "AUROC",  title = "")}

#Extract Features from blkbox algorithm
bb_ncv_feat<-function(dat,imp_thresh=0.15){
  model_ncv_allfeat=readRDS(dat)
  fi_ncv=data.frame(model_ncv_allfeat$MeanFeatureTable)
  fi_ncv=spread(fi_ncv,algorithm,Importance)
  fi_ncv$overallimp <- rowMeans(fi_ncv[c('GLM', 'nnet','party','randomforest','SVM')], na.rm=TRUE)
  feat_set_rf_ncv=fi_ncv[fi_ncv$overallimp>imp_thresh,]$feature
  as.vector(feat_set_rf_ncv)
  
}
