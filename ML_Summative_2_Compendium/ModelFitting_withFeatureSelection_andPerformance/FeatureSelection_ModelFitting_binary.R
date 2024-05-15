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
# These are the 5 we actually use
mutation<-readRDS("../preprocessed/mutations.Rda")
mrna<-readRDS("../preprocessed/minmrna.Rda")
clinical<-readRDS("../preprocessed/clinical.Rda")
geneprotlinked <- readRDS("../preprocessed/gene_prot_joined_set.Rda")

# Not used (directly) in analysis
# mirna<-readRDS("preprocessed/minmirna.Rda")
# protein<-readRDS("preprocessed/protein.Rda")
# V Large so considered separately
# methylation<-readRDS("preprocessed/methylation.Rda")

print("1")
# Gene feature selection

# Assuming `mrna` is your dataset and `pfi` is a binary outcome


## Function for plotting roc curves
plot_roc_curve <- function(predictions, actual, type, best) {
  roc_curve <- roc(response = actual, predictor = as.numeric(predictions), levels = rev(levels(actual)))
  auc_val <- auc(roc_curve)
  plot <- ggroc(roc_curve) + 
    ggtitle(paste("ROC Curve for", best, "with", type)) +
    annotate("text", x = 0.6, y = 0.2, label = paste("AUC:", round(auc_val, 2)), size = 5) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "red")
  return(plot)
}
# Function for fitting models
fit_models <- function(dataset, plot_roc=TRUE,set_number){
  set.seed(1) 
  list_of_sets <- c("mrna", "geneprotlinked", "clinical", "mutation")
  set_name <- list_of_sets[set_number]
  dataset$pfi <- as.factor(dataset$pfi)
  majority_class <- names(sort(table(dataset$pfi), decreasing = TRUE))[1]
  dataset$pfi[is.na(dataset$pfi)] <- majority_class
  # Ensure factor levels are valid R variable names
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
      po("filter", flt("mrmr"), filter.nfeat = 50) %>>%
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
  if (set_name == "clinical"){
    feat_names <- x$feature_types$id
  }else{
    feat_names <- x$feature_types$id[x$feature_types$type=="numeric"]
  }
  processed_df <- as.data.frame(setup_pipe_train$predict(train.gene.task)$nop.output$data())[,c(feat_names,"pfi")]
  print(names(processed_df))
  #Define formula from feature selection
  formula <- as.formula(paste("pfi ~", paste(names(processed_df)[names(processed_df) != "pfi"], collapse = " + ")))
  # Train an elnet model
  elnet_model = train(
    formula, data = processed_df,
    method = "glmnet",
    trControl = trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary),
    tuneLength = 10,
    metric = "ROC"
  )
  # Train an xgBoost
  print("set_name")
  xgb_model <- train(
    formula, data = processed_df,      
    method = "xgbTree",   
    trControl = trainControl(method = "cv", number = 2, classProbs = TRUE),
    tuneLength = 4,
    # Assign higher weight to minority class
    weights = ifelse(processed_df$pfi == 1, 5, 1)
  )
  print(set_name)
  layers <- c(2,2)
  # Train a Nerual Net
  nn_model <- neuralnet(
    formula, data = processed_df,
    hidden = layers,
    linear.output = FALSE)
  print("4")
  
  ### Testing
  #Process test data to fit our model to
  setup_pipe_test$train(test.gene.task)
  test_processed_df <- as.data.frame(setup_pipe_test$predict(test.gene.task)$nop.output$data())[,names(processed_df)]
  #Predict from an elnet model
  elnet_predictions_prob <- predict(elnet_model, test_processed_df[names(test_processed_df) != "pfi"], type = "prob")[, 2]
  # Adjust to 0.25 as form of 'model tuning'
  elnet_predictions <- ifelse(elnet_predictions_prob > 0.2, 1, 0)
  levels(test_processed_df$pfi)<-c(0,1)
  levels(processed_df$pfi)<-c(0,1)
  elnet_conf_matrix <- confusionMatrix(factor(elnet_predictions), 
                                    factor(test_processed_df$pfi))
  elnet_accuracy <- sum(diag(elnet_conf_matrix$table)) / sum(elnet_conf_matrix$table)
  elnet_recall <- sum(elnet_conf_matrix$table[2,2]) / sum(elnet_conf_matrix$table[2,])
  elnet_precision <- sum(elnet_conf_matrix$table[2,2]) / sum(elnet_conf_matrix$table[,2])
  elnet_f1_score <- (2*elnet_precision*elnet_recall)/(elnet_precision+elnet_recall)
  print(elnet_conf_matrix$table)
  print(paste("Elastic Net Accuracy:", elnet_accuracy, ", F1-Score: ", elnet_f1_score))
  # Predict from xgboost model
  xgb_predictions_prob <- predict(xgb_model, test_processed_df[names(test_processed_df) != "pfi"], type = "prob")[, 2]
  xgb_predictions <- ifelse(xgb_predictions_prob > 0.5, 1, 0)
  xgb_conf_matrix <- confusionMatrix(factor(xgb_predictions), 
                                       factor(test_processed_df$pfi))
  xgb_accuracy <- sum(diag(xgb_conf_matrix$table)) / sum(xgb_conf_matrix$table)
  xgb_recall <- sum(xgb_conf_matrix$table[2,2]) / sum(xgb_conf_matrix$table[2,])
  xgb_precision <- sum(xgb_conf_matrix$table[2,2]) / sum(xgb_conf_matrix$table[,2])
  xgb_f1_score <- (2*xgb_precision*xgb_recall)/(xgb_precision+xgb_recall)
  print(xgb_conf_matrix$table)
  print(paste("XGBoost Accuracy:", xgb_accuracy, ", F1-Score: ", xgb_f1_score))
  
  # Predict from neuralnet model
  nn_predictions <- compute(nn_model, test_processed_df[names(test_processed_df) != "pfi"])
  nn_predicted_values <- nn_predictions$net.result
  #Use 0.5 as threshold for prediction
  nn_predicted_classes <- ifelse(nn_predicted_values > 0.5, 1, 0)
  #Fit confusion matrix and calculate accuracy
  nn_conf_matrix <- confusionMatrix(factor(nn_predicted_classes[,2]), 
                                 factor(test_processed_df$pfi))
  nn_accuracy <- sum(diag(nn_conf_matrix$table)) / sum(nn_conf_matrix$table)
  nn_recall <- sum(nn_conf_matrix$table[2,2]) / sum(nn_conf_matrix$table[2,])
  nn_precision <- sum(nn_conf_matrix$table[2,2]) / sum(nn_conf_matrix$table[,2])
  nn_f1_score <- (2*nn_precision*nn_recall)/(nn_precision+nn_recall)
  print(nn_conf_matrix$table)
  print(paste("MLP neural net Accuracy:", nn_accuracy, ", F1-Score: ", nn_f1_score))
  print("5")
  print(paste("Factors selected", list(names(processed_df)[-1]))) 
  if (plot_roc==TRUE){
    f1s <- c(elnet=elnet_f1_score,xgboost = xgb_f1_score,neuralnet =nn_f1_score) 
    best_model <- names(which.max(f1s[!is.nan(f1s)]))
    if(is.null(best_model)){
      best_model <- "elnet"
    }
    print(paste("Best model is:", best_model))
    #define true class labels
    trueclasses <- factor(test_processed_df$pfi)
    trueclassestrain <- factor(processed_df$pfi) 
    # Plot ROC based on the best model
    if (best_model == "elnet") {
      elnet_train_predictions <- predict(elnet_model, processed_df[names(processed_df) != "pfi"], type = "prob")[, 2]
      plot_train <- plot_roc_curve(as.numeric(elnet_train_predictions), trueclassestrain, "train", best = best_model)
      plot_test <- plot_roc_curve(as.numeric(elnet_predictions_prob), trueclasses, "test", best = best_model)
    } else if (best_model == "xgboost") {
      xgb_train_predictions <- predict(xgb_model, processed_df[names(processed_df) != "pfi"], type = "prob")[, 2]
      plot_train <- plot_roc_curve(as.numeric(xgb_train_predictions), trueclassestrain, "train", best = best_model)
      plot_test <- plot_roc_curve(as.numeric(xgb_predictions_prob), trueclasses, "test", best = best_model)
    } else if (best_model == "neuralnet") {
      nn_train_predictions <- compute(nn_model, processed_df[names(processed_df) != "pfi"])
      nn_train_predicted_values <- nn_train_predictions$net.result
      plot_train <- plot_roc_curve(as.numeric(nn_train_predicted_values[, 2]), trueclassestrain, "train", best = best_model)
      plot_test <- plot_roc_curve(nn_predicted_values[, 1], trueclasses, "test", best = best_model)
    }
    ggsave(filename = paste0("../plot/",set_name, "train_ROC.png"), plot = plot_train)
    ggsave(filename = paste0("../plot/",set_name, "test_ROC.png"), plot = plot_test)
  } 
  return((c(elacc = elnet_accuracy, el_f1 = elnet_f1_score, xgbacc = xgb_accuracy,xgb_f1=xgb_f1_score,nnacc = nn_accuracy,nn_f1=nn_f1_score)))
  print("6")
}



# Set up list of datasets so we can parellalise
datasets <- list(as.data.frame(mrna),
                 as.data.frame(geneprotlinked),
                 as.data.frame(clinical),
                 as.data.frame(mutation)
                 )
list_of_sets <- c("mrna", "geneprotlinked", "clinical", "mutation")


## Attempt to parallelise this task
results <-cbind(as.data.frame(mclapply(1:4, function(x) {print(list_of_sets[x])
  fit_models(datasets[[x]],plot_roc = TRUE,set_number = x)}, mc.cores=4)))
fit_models(clinical,plot_roc = T,set_number=3)
names(clinical)
# Marker
print("endish")

results<-as.data.frame(rbind(results[,1],results[,2],results[,3]))
write.csv(results, "../plot/Result_Output.csv")
