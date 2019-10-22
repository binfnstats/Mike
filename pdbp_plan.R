### The purpose of this is to run caret models and blkbox (using Drake) on different feature sets and store the outputs.
#Still need to get xgboost and gb to produce output
#separately need to create output dataset and predict on that


setwd(proc_data)
plan <- drake_plan(
  #load and transform data
  data_mod=load_PDBP("mod_data.rds","PD_55_sub_labels.rds","bb",20,"BF"),
  #Bor_features=Boruta_feat(data_mod),
  inTrain = createDataPartition(y = data_mod$labs, p = .7, list = FALSE),
  training = data_mod[ inTrain,],
  testing = data_mod[ -inTrain,],
  
  #nih data bring in
  mod_data=readRDS("mod_data.rds"),
  nih_data=mod_data[mod_data$Q2600=="nih",],
  #run caret models, predict and score nih data
  mod_preds=mod_predict(training,testing,nih_data),
  gbmxgb_preds=gbmxgb_pred(training,testing),
  #show caret diagnostics
  plot_roc=mod_diag(mod_preds,testing,"ROC features"),
  mod=plot_roc
  #run blkbox and show diagnostics
  #,blkbox_run(data_mod,BF_dat_labs)
  )
# diagnose(mod_preds)
#make plan and configure etc.
make(plan)
  config <- drake_config(plan)
vis_drake_graph(config)

nih_preds=readRDS("nih_preds.RDS")
write.csv(nih_preds,"nih_preds.csv")
View(nih_preds)

#pulling the key AUC values out and generating boxplots
aucvals=readRDS("aucvals.RDS")
aucvals=gather(aucvals,model,value)
View(aucvals)

res_mlist=readRDS("res_mlist.RDS")
bplot(res_mlist)

#blkbox code - should be able to get blkbox into a function, failed so far.

data_mod=load_PDBP("mod_data.rds","PD_55_sub_labels.rds","all",2599,"BF")
bb_ncv_feat=bb_ncv_feat("model_ncv_allfeat.rds",0.15)
mod_data=data_mod[,c(bb_ncv_feat)]
mod_data=mod_data[,colnames(mod_data)!="labs"]
#labels=mod_data[,colnames(mod_data)=="labs"]
pd55_sub_labels<-readRDS("PD_55_sub_labels.rds")
pd55_sub_part=Partition(mod_data,labels=pd55_sub_labels,size=.6,seed=1234)
bb <- blkbox(data = pd55_sub_part,exclude="xgboost")
# Performance
bb_perf = Performance(bb)
blkboxROC(bb_perf)
bb_ncv<- blkboxNCV(data = mod_data,
                   labels = pd55_sub_labels,
                   Method = "randomforest",exclude = "xgboost",
                   AUC = 0.8,metric="AUROC")
ncv.plot(bb_ncv, metric = "AUROC",  title = "")


##### NIH output data

setwd(proc_data)
nih_data<-readRDS("mod_data.rds")
nih_data<-mod_data[mod_data$Q2600=="nih",]
nih_preds <- lapply(model_list, predict, newdata=nih_val, type="prob")



