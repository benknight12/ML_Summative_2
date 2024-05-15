

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

methylation<-readRDS("../preprocessed/methylation.Rda")
mutation<-readRDS("../preprocessed/mutations.Rda")
mrna<-readRDS("../preprocessed/minmrna.Rda")
linked <- readRDS("../preprocessed/gene_prot_joined_set.Rda")
clinical<-readRDS("../preprocessed/clinical.Rda")


## Conducting PCA on the full dataset would take far too long so we will perform it on a shrunk version using the features selected from MIFS to see if we can obtain better feature selection
methylation_names <-c("cg00573880", "cg00910503", "cg01227027", "cg01254170", 
"cg01881939", "cg02033323", "cg02841571", "cg03519441", "cg04120542", "cg04144589",
 "cg04295632", "cg04864609", "cg05302420", "cg06561166", "cg08080498", "cg08128789",
 "cg08606356", "cg09113768", "cg09118017", "cg09551172", "cg10874502", "cg11989257",
 "cg11992171", "cg12278631", "cg12496657", "cg13002957", "cg13904493", "cg13982688",
 "cg14284187", "cg14574905", "cg15703756", "cg16101101", "cg16278877", "cg16376691",
 "cg16995290", "cg17985124", "cg19083143", "cg19685285", "cg19800407", "cg20539961",
 "cg20816889", "cg21460053", "cg21915935", "cg22864414", "cg24140362", "cg24540003",
 "cg24553417", "cg24920145", "cg25498838", "cg25704407", "cg27316970", "cg27547442")



mutation_names <- c("ABCB10", "ARVCF", "BOD1L1", "C5orf52", "CARD6", "CATSPERB", "CCT2", 
"CDHR4", "CLEC17A", "COL6A4P1", "CTNNB1", "CX3CR1", "CYP2A13", "DCST1", "DNMT3A", "DSG3",
 "ERBB2IP", "EVC2", "FBXO44", "GALNT6", "GMPS", "GNB2L1", "GPATCH4", "GPR1", "HIPK2", 
"IGSF9B", "IL12RB1", "ITGA7", "KYNU", "LCMT2", "LEMD2", "LRFN2", "MECOM", "MTHFD1L", 
"MTR", "NCBP2L", "NLGN1", "PHF14", "PTPRO", "SDPR", "SI", "TACC2", "TCEA3", "TIMM50",
 "TRADD", "WDR7", "ZDHHC1", "ZNF280B", "ZNF432", "ZNF595")

mrna_names <- c("ABCG8",       "ACCSL",       "ADH1B",       "ADIG",     "ALX4",
"APOC2",           "ART3",        "AVPR1B",	 "BARHL1",    "CCDC117",
"CDH4",        "CGB2",        "CNR2",        "CNTD2",     "CRISP2",
"CRISP3",          "FLJ44054",    "HRASLS2",     "HTR1F",     "IGF1R",
"KCNK2",           "KRTAP19.1",   "KRTAP4.11",   "LASS3",     "LOC146336",
"LRRC69",          "NAV3",        "NCRNA00028",  "NRXN2",     "NXPH1",
"PCDHB1",          "PGLYRP2",     "PLGLA",	 "PSD2",      "RIMS3",
"SAG",         "SCUBE2",         "SLC26A2",     "SLC38A11",   "SLPI",
"SLURP1",          "SNORA14B",    "SNORD115.13", "SPRR2E",  "TBC1D3H",
"TBPL2",           "TRPC7",	 "VNN3" ,    "WDR72",       "YBX2")

linked_names <- c("ACVRL1.x", "ARID1A.x", "ARID1A.y", "ASNS.x",   "ATM.x",    "ATM.y",
"BRCA2.y","BRD4.x",   "CDK1.x",   "CDK1.y",   "COG3.x",  "DIRAS3.x",
"DUSP4.x",  "EGFR.x",   "EGFR.y",   "EPPK1.x",  "EPPK1.y","ERCC1.x",
"ERCC1.y",  "ERCC5.x",  "ERCC5.y",  "FASN.x",   "G6PD.x",  "G6PD.y",
"GAB2.x",   "GAB2.y",   "IGFBP2.x", "IGFBP2.y", "INPP4B.x","IRS1.x",
"IRS1.y",   "MYH11.x",  "NF2.x",    "NF2.y",    "PCNA.y",  "PDCD4.x",
"PDK1.x",   "PEA15.x",  "PEA15.y",  "PRDX1.x",  "PREX1.x", "PTEN.x",
"RBM15.x",  "SCD.x",    "SCD.y",    "TAZ.x",    "TFRC.x",  "XBP1.y",
"XRCC1.x",  "XRCC1.y")

mrna <- mrna[,mrna_names]
methylation <- methylation[,methylation_names]
mutation <- mutation[,mutation_names]
linked <- linked[,linked_names]


names(mrna)
print(dim(mrna))

names(methylation) <- lapply(names(methylation), function(x) x <- paste0(x,"meth"))
names(mutation) <- lapply(names(mutation), function(x) x <- paste0(x,"mut"))
names(mrna) <- lapply(names(mrna), function(x) x <- paste0(x,"mrna"))
names(linked) <- lapply(names(linked), function(x) x <- paste0(x,"linked"))
names(clinical) <- lapply(names(clinical), function(x) x <- paste0(x,"clin"))

methylation$names<- rownames(methylation)
mutation$names<- rownames(mutation)
linked$names<- rownames(linked)
mrna$names<- rownames(mrna)
clinical$names<- rownames(clinical)


noclinical_full_merge<-(merge(mrna, merge(mutation,merge(clinical[,c("pficlin","names")],merge(methylation,linked, by = "names"), by = "names"), by = "names"), by="names"))
full_merge<-(merge(mrna, merge(mutation,merge(clinical,merge(methylation,linked, by = "names"), by = "names"), by = "names"), by="names"))

print(dim(full_merge))
dataset<-full_merge
dataset$pfi <- as.factor(dataset$pficlin)
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

# Prepare data
setup_pipe_train <- (
  po("imputesample") %>>%
    po("scale") %>>%
    po("pca") %>>%
    po("nop")
)
setup_pipe_test <- (
  po("imputesample") %>>%
    po("scale") %>>%
    po("nop")
)


### Training
# Obtain a processed training set
setup_pipe_train$train(train.gene.task)
x<-setup_pipe_train$train(train.gene.task)$nop.output
features<-cbind(lapply(1:30,function(x) paste0("PC",x)))
paste(rbind(features,"pfi"))
processed_df <- as.data.frame(setup_pipe_train$predict(train.gene.task)$nop.output$data())[,paste(rbind(features,"pfi"))]
print(names(processed_df))
head(processed_df)
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
xgb_model <- train(
  formula, data = processed_df,
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 2, classProbs = TRUE),
  tuneLength = 4,
  # Assign higher weight to minority class
  weights = ifelse(processed_df$pfi == 1, 5, 1)
)
layers <- c(2,1)
if(names(dataset)[1]==names(methylation)[1]){
  layers <- c(1,1)
}
# Train a Nerual Net

nn_model <- neuralnet(
  formula, data = processed_df,
  hidden = layers,
  linear.output = FALSE)


### Testing


#Process test data to fit our model to
transformed_task_test <- setup_pipe_train$predict(test.gene.task)
test_processed_df <-as.data.frame(transformed_task_test$nop.output$data())[,names(processed_df)]#Predict from an elnet model
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
print(paste("Factors selected", list(names(processed_df)[-1])))

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
ggsave(filename = "../plot/merged_train_ROC.png", plot = plot_train)
ggsave(filename = "../plot/merged_test_ROC.png", plot = plot_test)


answers<-as.data.frame(c(elacc = elnet_accuracy, el_f1 = elnet_f1_score, xgbacc = xgb_accuracy,xgb_f1=xgb_f1_score,nnacc = nn_accuracy,nn_f1=nn_f1_score))


## Fit without clinical data
print(dim(noclinical_full_merge))
dataset<-noclinical_full_merge
dataset$pfi <- as.factor(dataset$pficlin)
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

# Prepare data
setup_pipe_train <- (
  po("imputesample") %>>%
    po("scale") %>>%
    po("pca") %>>%
    po("nop")
)
setup_pipe_test <- (
  po("imputesample") %>>%
    po("scale") %>>%
    po("nop")
)


### Training
# Obtain a processed training set
setup_pipe_train$train(train.gene.task)
x<-setup_pipe_train$train(train.gene.task)$nop.output
features<-cbind(lapply(1:30,function(x) paste0("PC",x)))
paste(rbind(features,"pfi"))
processed_df <- as.data.frame(setup_pipe_train$predict(train.gene.task)$nop.output$data())[,paste(rbind(features,"pfi"))]
print(names(processed_df))
head(processed_df)
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
xgb_model <- train(
  formula, data = processed_df,
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 2, classProbs = TRUE),
  tuneLength = 4,
  # Assign higher weight to minority class
  weights = ifelse(processed_df$pfi == 1, 5, 1)
)
layers <- c(2,1)
if(names(dataset)[1]==names(methylation)[1]){
  layers <- c(1,1)
}
# Train a Nerual Net
nn_model <- neuralnet(
  formula, data = processed_df,
  hidden = layers,
  linear.output = FALSE)


### Testing


#Process test data to fit our model to
transformed_task_test <- setup_pipe_train$predict(test.gene.task)
test_processed_df <-as.data.frame(transformed_task_test$nop.output$data())[,names(processed_df)]#Predict from an elnet model
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
print(paste("Factors selected", list(names(processed_df)[-1])))

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
ggsave(filename = "../plot/merged_train_ROC_noclin.png", plot = plot_train)
ggsave(filename = "../plot/merged_test_ROC_noclin.png", plot = plot_test)


answers<-as.data.frame(c(elacc = elnet_accuracy, el_f1 = elnet_f1_score, xgbacc = xgb_accuracy,xgb_f1=xgb_f1_score,nnacc = nn_accuracy,nn_f1=nn_f1_score))


