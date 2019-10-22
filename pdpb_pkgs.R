#Packages

library(caret)
library(corrplot)			# plot correlations
library(doParallel)		# parallel processing
library(dplyr)        # Used by caret
library(gbm)				  # GBM Models
library(pROC)				  # plot the ROC curve
library(xgboost)      # Extreme Gradient Boosting
library("caret")
library("mlbench")
library("pROC")
library(ROCR)
library("rpart")
library("caretEnsemble")
library("mlbench")
library("randomForest")
library("nnet")
library("caTools")
library("gbm")
library(Boruta)
library(blkbox)
library(drake)
library(limma)
library(tidyr)
library(tidyverse)
library(hrbrthemes)
install.packages("hrbrthemes")
library(viridis)
library(ggplot2)
proj<-"/Users/michaelallwright/Dropbox (Sydney Uni)/michael_PhD/Projects/Parkinson's Longitudinal Study/"
Analysis<-"/Users/michaelallwright/Dropbox (Sydney Uni)/michael_PhD/Projects/Parkinson's Longitudinal Study/Analysis/"
proc_data="/Users/michaelallwright/Dropbox (Sydney Uni)/michael_PhD/Projects/Parkinson's Longitudinal Study/Data/Processed/"
mod_runs="/Users/michaelallwright/Dropbox (Sydney Uni)/michael_PhD/Projects/Parkinson's Longitudinal Study/Data/Processed/Model Runs/"
setwd(proc_data)