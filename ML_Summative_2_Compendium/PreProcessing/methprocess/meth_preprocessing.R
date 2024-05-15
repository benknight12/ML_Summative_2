library(parallel)
library(praznik)
library(mlr3pipelines)
library(mlr3filters)
library(mlr3learners)
library(mlr3tuning)
library(caret)
library(mlr3extralearners)
library(mlr3viz)
library(precrec) ## mlr3 plots
library(rpart) ## classification trees
library(ranger) ## random forest
library(glmnet) ## lasso
library(igraph) ## mlr3 pipelines
library(data.table)
library(neuralnet)
library(pROC)
library(ggplot2)

methylation <- readRDS("../../preprocessed/methylation.Rda")
dataset <- methylation
majority_class <- names(sort(table(dataset$pfi), decreasing = TRUE))[1]
dataset$pfi[is.na(dataset$pfi)] <- majority_class
# Ensure factor levels are valid R variable names
dataset$pfi <- as.factor(dataset$pfi)
levels(dataset$pfi) <- make.names(levels(dataset$pfi))
split <- createDataPartition(dataset$pfi, p = 0.75, list = FALSE)
train_data <- dataset[split, ]
test_data <- dataset[-split, ]
saveRDS(train_data,"../../preprocessed/meth_training_full.Rda")
saveRDS(test_data,"../../preprocessed/meth_testing_full.Rda")
