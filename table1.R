#Table 1
library(tableone)

final_rf_ben2 <- final_rf_ben1
final_rf_ben2$contrast_imaging <- ifelse(final_rf_ben1$contr_ct + final_rf_ben1$contr_mr > 0, 1,0)
var_names <- colnames(final_rf_ben2[,c(7,9,11:18,20:58,68:69,71:72,81)])

factor_vars <- c("race_adj")
tab1 <- CreateTableOne(vars = var_names, data = final_rf_ben2,factorVars = factor_vars, strata = "contrast_imaging",addOverall = TRUE)
tab1Mat <- print(tab1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(tab1Mat, file = "P:/pacce/s4k/MIMIC4/table1.csv")
