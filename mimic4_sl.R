#SL analysis
load("~/Projects/Research/MIMIC4/code/mimic4_sl.RData")
set.seed(42)

ct_sl <- SuperLearner(Y = final_rf_ben1$contr_ct, X = final_rf_ben1[,c(7:9,11:18,20:58,68:69)], family = binomial(),
                         SL.library = c("SL.ranger","SL.glm","SL.ksvm","SL.xgboost","SL.glmnet"))

mri_sl <- SuperLearner(Y = final_rf_ben1$contr_mri, X = final_rf_ben1[,c(7:9,11:18,20:58,68:69)], family = binomial(),
                      SL.library = c("SL.ranger","SL.glm","SL.ksvm","SL.xgboost","SL.glmnet"))

mort_ct_sl <- SuperLearner(Y = final_rf_ben1$death, X = final_rf_ben1[,c(7:9,11:18,20:58,68:69,79)], family = binomial(),
                      SL.library = c("SL.ranger","SL.glm","SL.ksvm","SL.xgboost","SL.glmnet"))

mort_mri_sl <- SuperLearner(Y = final_rf_ben1$death, X = final_rf_ben1[,c(7:9,11:18,20:58,68:69,80)], family = binomial(),
                       SL.library = c("SL.ranger","SL.glm","SL.ksvm","SL.xgboost","SL.glmnet"))

adm_ct_sl <- SuperLearner(Y = final_rf_ben1$readmit_num, X = final_rf_ben1[,c(7:9,11:18,20:58,68:69,79)], family = gaussian(),
                           SL.library = c("SL.ranger","SL.glm","SL.ksvm","SL.xgboost","SL.glmnet"))

adm_mri_sl <- SuperLearner(Y = final_rf_ben1$readmit_num, X = final_rf_ben1[,c(7:9,11:18,20:58,68:69,80)], family = gaussian(),
                            SL.library = c("SL.ranger","SL.glm","SL.ksvm","SL.xgboost","SL.glmnet"))

save.image("P:/pacce/s4k/MIMIC4/mimic4_sl.RData")

sl <- ct_sl
X <- final_rf_ben1[,c(7:9,11:18,20:58,68:69)]
Y <- final_rf_ben1$contr_ct



#test prediction
pred = predict(sl, X, onlySL = TRUE)
compiled_accuracy <- mean(round(pred$pred)==Y) #90.2%

#filter prognostic factors
library(vip)
imp_fun <- function(object, newdata) { # for permutation-based VI scores
  predict(object, newdata = newdata)$pred
}
set.seed(42)
var_imp1 <- vip::vi(sl, method = "permute", train = X, target = Y, metric = "rmse",
                    pred_wrapper = imp_fun, nsim = 5)

barplot(var_imp1$Importance[1:5], main="CT Use Variable Importance", horiz=TRUE,
        names.arg=c("Past Meds","Heart Rate","Systolic BP","Pain","New Meds"),cex.names=1.1,col="brown3", xlab = "Variable Importance")

pred_rocr = ROCR::prediction(pred$pred, Y)
auc = ROCR::performance(pred_rocr, measure = "auc", x.measure = "cutoff")@y.values[[1]] #
perf <- ROCR::performance(pred_rocr, "tpr", "fpr")
plot(perf)


#ROC
library(ROCR)

g_predictions <- cbind(sl$SL.predict, sl$library.predict[,1:5])
g_pred.mat <- prediction(g_predictions, labels = matrix(Y,
                                                   nrow = length(Y),
                                                   ncol = 6))
g_perf.mat <-performance(g_pred.mat,"tpr","fpr")

plot(g_perf.mat, main = "CT with Contrast Prediction ROC" , col = as.list(c("black","red","darkgreen","blue","orange","purple")), xlab = "1 - Specificity", ylab = "Sensitivity")
legend("bottomright",legend = c("SuperLearner","Random Forest","Logistic Regression",
                                "KSVM","XGBoost","Lasso Regression"),
       col = c("black","red","darkgreen","blue","orange","purple"), cex = 1, lty = 1)