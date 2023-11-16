library(forestplot)
library(dplyr)
# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <- structure(list(mean  = c(NA, NA, 0.97, 0.76, NA, 3.51, 2.87), 
                                      lower = c(NA, NA, 0.43, 0.27, NA, 3.32, 2.67),
                                      upper = c(NA, NA, 1.51, 1.25, NA, 3.71, 3.06)),
                                 .Names = c("mean", "lower", "upper"), 
                                 row.names = c(NA, -11L), 
                                 class = "data.frame")

tabletext <- cbind(c("", "Risk of Mortality", "          CT with Contrast Use", "          MRI with Contrast Use",
                     "Number of Readmissions", "          CT with Contrast Use", "          MRI with Contrast Use"),
                   c("Difference Estimate", "", "0.97%", "0.76%", "", "3.51", "2.87"),
                   c("95% CI", "", "(0.43%, 1.51%)", "(0.27%, 1.25%)","", "(3.32, 3.71)", "(2.67, 3.06)"))

cochrane_from_rmeta %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(rep(TRUE, 2), rep(FALSE, 2), TRUE,rep(FALSE, 2)),
             clip = c(0, 4),
             xticks = seq(from = 0, to = 4, by = 0.25),
             xlog = FALSE,
             xlab = "Difference Estimate",
             txt_gp = fpTxtGp(ticks=gpar(cex=0.7),xlab=gpar(cex=0.8)),
             boxsize = 0.07,
             align = c("l","c","c"),
             col = fpColors(box = "black",
                            line = "grey",
                            summary = "black"))