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

dat0<-readRDS("../preprocessed/MethData_following_feature_selection/TrimData_0.Rda")
dat1<-readRDS("../preprocessed/MethData_following_feature_selection/TrimData_1.Rda")
dat2<-readRDS("../preprocessed/MethData_following_feature_selection/TrimData_2.Rda")
dat3<-readRDS("../preprocessed/MethData_following_feature_selection/TrimData_3.Rda")
dim(dat1)
dat1<-dat1[,-ncol(dat1)]
dat2<-dat2[,-ncol(dat2)]
dat3<-dat3[,-ncol(dat3)]
dim(dat1)
processed_df<- cbind(dat0,dat1,dat2,dat3)

plot_roc_curve <- function(predictions, actual, type, best) {
  roc_curve <- roc(response = actual, predictor = as.numeric(predictions), levels = rev(levels(actual)))
  auc_val <- auc(roc_curve)
  plot <- ggroc(roc_curve) + 
    ggtitle(paste("ROC Curve for", best, "with", type)) +
    annotate("text", x = 0.6, y = 0.2, label = paste("AUC:", round(auc_val, 2)), size = 5) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "red")
  return(plot)
}

f1Score <- function(data, lev = NULL, model = NULL) {
  precision <- posPredValue(data$pred, data$obs, positive = lev[2])
  recall <- sensitivity(data$pred, data$obs, positive = lev[2])
  f1 <- ifelse((precision + recall) > 0, (2 * precision * recall) / (precision + recall), NA)
  return(c(F1 = f1))
}

dim(processed_df)
processed_df <- processed_df[, !duplicated(colnames(processed_df), fromLast = TRUE)] 

processed_df$pfi <- as.factor(processed_df$pfi)
majority_class <- names(sort(table(processed_df$pfi), decreasing = TRUE))[1]
processed_df$pfi[is.na(processed_df$pfi)] <- majority_class
# Ensure factor levels are valid R variable names
levels(processed_df$pfi) <-c("X0","X1")
levels(processed_df$pfi) <- make.names(levels(processed_df$pfi))

head(processed_df)
table(is.na(processed_df))
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
table(is.na(processed_df$cg00573880))
# Train an xgBoost
xgb_model <- train(
  formula, data = processed_df,
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 2, classProbs = TRUE),
  tuneLength = 4,
  # Assign higher weight to minority class
  weights = ifelse(processed_df$pfi == 1, 5, 1)
)  
layers <- c(2,2)
# Train a Nerual Net
nn_model <- neuralnet(
  formula, data = processed_df,
  hidden = layers,
  linear.output = FALSE)
print("4")

### Testing
## Create raw test data of just features identified from training data feature selection
raw_test_data<-readRDS("testing_full.Rda")
test_data<-raw_test_data[,names(processed_df)]
test.gene.task <- TaskClassif$new(
  backend = test_data,
  target = "pfi",
  id = "Breast cancer PFI")
setup_pipe_test <- (
  po("imputesample") %>>%
    po("scale") %>>%
    po("nop")
)
#Process test data to fit our model to
setup_pipe_test$train(test.gene.task)
test_processed_df <- as.data.frame(setup_pipe_test$predict(test.gene.task)$nop.output$data())
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
xgb_predictions <- ifelse(xgb_predictions_prob > 0.3, 1, 0)
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
nn_predicted_classes <- ifelse(nn_predicted_values > 0.3, 1, 0)
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

best_model <- names(which.max(c(elnet = elnet_accuracy, xgboost = xgb_accuracy, neuralnet = nn_accuracy)))
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
ggsave(filename = "../plot/methtrain_ROC.png", plot = plot_train)
ggsave(filename = "..plot/methtest_ROC.png", plot = plot_test)



results<-c(elacc = elnet_accuracy, el_f1 = elnet_f1_score, xgbacc = xgb_accuracy,xgb_f1=xgb_f1_score,nnacc = nn_accuracy,nn_f1=nn_f1_score)

print("endish")
write.csv(results, "../plot/Meth_Result_Output.csv")
