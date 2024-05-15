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

## Load the preprocessed datasets (linked to pfi and redundant columns removed)
methylation<-readRDS("../../preprocessed/meth_training_full.Rda")

plot_roc_curve <- function(predictions, actual, type, best) {
  roc_curve <- roc(response = actual, predictor = as.numeric(predictions), levels = rev(levels(actual)))
  auc_val <- auc(roc_curve)
  plot <- ggroc(roc_curve) +
    ggtitle(paste("ROC Curve for", best, "with", type)) +
    annotate("text", x = 0.6, y = 0.2, label = paste("AUC:", round(auc_val, 2)), size = 5) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "red")
  return(plot)
}
print("1b")
# Function for fitting models
dataset <- cbind(methylation[,c(300000:394980)],pfi =methylation[,"pfi"])

majority_class <- names(sort(table(dataset$pfi), decreasing = TRUE))[1]
dataset$pfi[is.na(dataset$pfi)] <- majority_class
  # Ensure factor levels are valid R variable names
dataset$pfi <- as.factor(dataset$pfi)
levels(dataset$pfi) <- make.names(levels(dataset$pfi))
split <- createDataPartition(dataset$pfi, p = 0.75, list = FALSE)
train_data <- dataset[split, ]
test_data <- dataset[-split, ]
train.gene.task <- TaskClassif$new(
  backend = train_data,
  target = "pfi",
  id = "Breast cancer PFI")
test.gene.task <- TaskClassif$new(
  backend = test_data,
  target = "pfi",
  id = "Breast cancer PFI")
print("2")
# Prepare data
setup_pipe_train <- (
  po("imputesample") %>>%
    po("scale") %>>%
    po("filter", flt("mrmr"), filter.nfeat = 13) %>>%
    po("classbalancing",
       id = "oversample", adjust = "major",
       reference = "major", shuffle = FALSE, ratio = 2) %>>%
    po("nop")
)
setup_pipe_test <- (
  po("imputesample") %>>%
    po("scale") %>>%
    po("nop")
)
print("3")

  ### Training
  # Obtain a processed training set
setup_pipe_train$train(train.gene.task)
x<-setup_pipe_train$train(train.gene.task)$nop.output
processed_df <- as.data.frame(x$data())[,c(x$feature_types$id[x$feature_types$type=="numeric"],"pfi")]
saveRDS(processed_df,paste0("../../preprocessed/MethData_following_feature_selection/TrimData_3.Rda"))


