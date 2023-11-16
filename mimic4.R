library(foreign)
library(corrplot)
library(comorbidity)
#library(mice)
library(tmle)


###LOAD files in
diag <- read.csv("P:/pacce/s4k/MIMIC4/ed/diagnosis.csv")
edstays <- read.csv("P:/pacce/s4k/MIMIC4/ed/edstays.csv")
medrecon <- read.csv("P:/pacce/s4k/MIMIC4/ed/medrecon.csv")
pyxis <- read.csv("P:/pacce/s4k/MIMIC4/ed/pyxis.csv")
triage <- read.csv("P:/pacce/s4k/MIMIC4/ed/triage.csv")
#vitals <- read.csv("P:/pacce/s4k/MIMIC4/ed/vitalsign.csv")

admissions <- read.csv("P:/pacce/s4k/MIMIC4/ed/admissions.csv")
drg <- read.csv("P:/pacce/s4k/MIMIC4/ed/drgcodes.csv")
rad <- read.csv("P:/pacce/s4k/MIMIC4/ed/radiology_detail.csv")

###PREPROCESSING

#imaging codes (A, C, D, E, G, S = x-ray, Q = CT, T = MRI, J, U = US, F, I, V = Fluoro)
#number of total imaging
rad_code <- subset(rad, field_name == "exam_code")
rad_code$one <- 1
rad_num <- aggregate(rad_code[,c("one")],by = list(rad_code$subject_id), FUN = sum)
colnames(rad_num) <- c("subject_id", "rad_num")

#number of total X-ray
rad_code$type <- substring(rad_code$field_value,1,1)
x_ray <- subset(rad_code,type %in% c("A", "C","D","E","G","S"))
x_ray_num <- aggregate(x_ray[,c("one")],by = list(x_ray$subject_id), FUN = sum)
colnames(x_ray_num) <- c("subject_id", "x_ray_num")

#number of total CT
ct <- subset(rad_code,type=="Q")
ct_num <- aggregate(ct[,c("one")],by = list(ct$subject_id), FUN = sum)
colnames(ct_num) <- c("subject_id", "ct_num")

#number of total MRI
mri <- subset(rad_code,type=="T")
mri_num <- aggregate(mri[,c("one")],by = list(mri$subject_id), FUN = sum)
colnames(mri_num) <- c("subject_id", "mri_num")

#number of total US
us <- subset(rad_code,type %in% c("J","U"))
us_num <- aggregate(us[,c("one")],by = list(us$subject_id), FUN = sum)
colnames(us_num) <- c("subject_id", "us_num")

#number of total fluoro
fluoro <- subset(rad_code,type %in% c("F","I","V"))
fluoro_num <- aggregate(fluoro[,c("one")],by = list(fluoro$subject_id), FUN = sum)
colnames(fluoro_num) <- c("subject_id", "fluoro_num")

#filter out ct w contrast and mr w contrast and create flags by subject id
ct_contr_desc <- unique(rad$field_value[grepl("contrast",rad$field_value,ignore.case = TRUE)&
                                          grepl("ct",rad$field_value,ignore.case = TRUE)])

ct_contr_desc1 <- ct_contr_desc[c(2:7,9,15,19,20,23,24,26:30,33,35:40)]

ct_contr_id <- unique(rad$subject_id[rad$field_value %in% ct_contr_desc1])

mr_contr_desc <- unique(rad$field_value[grepl("contrast",rad$field_value,ignore.case = TRUE)&
                                          grepl("mr",rad$field_value,ignore.case = TRUE)])

mr_contr_desc1 <- mr_contr_desc[c(1,3,4,6,9,10,13:15,19:24,26,28:30,32:35,
                                  37,40,43:49,51,55:59,62:69,71:73,76,78,80,82,
                                  83,85,86,89,91:93)]

mr_contr_id <- unique(rad$subject_id[rad$field_value %in% mr_contr_desc1])


#medrecon
#compile number of existing meds
medrecon$one <- 1
medrecon_num <- aggregate(medrecon[,c("one")],by = list(medrecon$subject_id, medrecon$stay_id), FUN = sum)
colnames(medrecon_num) <- c("subject_id","stay_id", "medrecon_num")

#pyxis
#compile number of new meds
pyxix_num <- aggregate(pyxis[,c("med_rn")],by = list(pyxis$subject_id, pyxis$stay_id), FUN = max)
colnames(pyxix_num) <- c("subject_id","stay_id", "pyxis_num")

##diagnoses
#first sanity check stay ids to ensure uniqueness
#for (i in 1:length(unique(diag$stay_id))){
#  id <- unique(diag$stay_id)[i]
#  sub <- length(unique(diag$subject_id[diag$stay_id == id]))
#  if (sub > 1){
#    print(sub)
#    print(id)
#  }
#}

#number of diagnoses
diag$one <- 1
diag_num <- aggregate(diag[,c("one")],by = list(diag$subject_id, diag$stay_id), FUN = sum)
colnames(diag_num) <- c("subject_id","stay_id", "diag_num")

#classify diagnoses by comorbidity package
diag_9 <- subset(diag, icd_version==9)
diag_10 <- subset(diag, icd_version==10)

diag9_comorb <- comorbidity(diag_9,id="stay_id",code="icd_code",
                            map = "elixhauser_icd9_quan", assign0=TRUE)
diag10_comorb <- comorbidity(diag_10,id="stay_id",code="icd_code",
                            map = "elixhauser_icd10_quan", assign0=TRUE)
diag_comorb <- rbind(diag9_comorb,diag10_comorb)

#compute elixhauser scores
elix_comorb <- diag_comorb
elix_comorb$w_elix <- score(diag_comorb,weights = "vw",assign0 = TRUE)
elix_comorb$uw_elix <- score(diag_comorb,weights = NULL,assign0 = TRUE)

#drg aggregation
sub_drg_sev <- drg[!is.na(drg$drg_severity),]
drg_severity <- aggregate(sub_drg_sev[,c("drg_severity")],by = list(sub_drg_sev$subject_id, sub_drg_sev$hadm_id), FUN = mean)
colnames(drg_severity) <- c("subject_id","hadm_id", "drg_severity")
drg_mortality <- aggregate(sub_drg_sev[,c("drg_mortality")],by = list(sub_drg_sev$subject_id, sub_drg_sev$hadm_id), FUN = mean)
colnames(drg_mortality) <- c("subject_id","hadm_id", "drg_mortality")

#merge ed datasets and compile useful columns (subject_id, stay_id)
combined_ed <- merge(edstays, triage, by = c("subject_id","stay_id"))
combined_ed <- merge(combined_ed, drg_severity, by = c("subject_id","hadm_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, drg_mortality, by = c("subject_id","hadm_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, elix_comorb, by = c("stay_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, diag_num, by = c("subject_id","stay_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, pyxix_num, by = c("subject_id","stay_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, medrecon_num, by = c("subject_id","stay_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, rad_num, by = c("subject_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, x_ray_num, by = c("subject_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, ct_num, by = c("subject_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, mri_num, by = c("subject_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, us_num, by = c("subject_id"), all.x= TRUE)
combined_ed <- merge(combined_ed, fluoro_num, by = c("subject_id"), all.x= TRUE)

#check covariate columns
#adjust race
combined_ed$race_adj <- combined_ed$race
combined_ed$race_adj <- ifelse(combined_ed$race %in% c("WHITE","WHITE - EASTERN EUROPEAN",
                                                       "WHITE - RUSSIAN", "WHITE - BRAZILIAN",
                                                       "WHITE - OTHER EUROPEAN") ,"WHITE",combined_ed$race_adj)
combined_ed$race_adj <- ifelse(combined_ed$race %in% c("BLACK/AFRICAN","BLACK/AFRICAN AMERICAN",
                                                       "BLACK/CAPE VERDEAN","BLACK/CARIBBEAN ISLAND") ,"BLACK",combined_ed$race_adj)
combined_ed$race_adj <- ifelse(combined_ed$race %in% c("HISPANIC OR LATINO","HISPANIC/LATINO - CENTRAL AMERICAN",
                                                       "HISPANIC/LATINO - COLUMBIAN","HISPANIC/LATINO - CUBAN",
                                                       "HISPANIC/LATINO - DOMINICAN","HISPANIC/LATINO - GUATEMALAN",
                                                       "HISPANIC/LATINO - HONDURAN","HISPANIC/LATINO - MEXICAN",
                                                       "HISPANIC/LATINO - PUERTO RICAN","HISPANIC/LATINO - SALVADORAN",
                                                       "PORTUGUESE","SOUTH AMERICAN") ,"HISPANIC",combined_ed$race_adj)
combined_ed$race_adj <- ifelse(combined_ed$race %in% c("ASIAN","ASIAN - ASIAN INDIAN","ASIAN - CHINESE",
                                                       "ASIAN - KOREAN","ASIAN - SOUTH EAST ASIAN",
                                                       "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER") ,"AAPI",combined_ed$race_adj)
combined_ed$race_adj <- ifelse(combined_ed$race %in% c("MULTIPLE RACE/ETHNICITY","PATIENT DECLINED TO ANSWER",
                                                       "UNABLE TO OBTAIN","UNKNOWN") ,"OTHER",combined_ed$race_adj)

#adjust transport
combined_ed$arrival_transport <- ifelse(combined_ed$arrival_transport == "UNKNOWN","OTHER",combined_ed$arrival_transport)

#adjust pain
combined_ed$pain <- as.numeric(combined_ed$pain)

#convert NA to 0 for num_cols
num_cols <- colnames(combined_ed)[54:62]
num_df <- combined_ed[,num_cols]
num_df[is.na(num_df)] <- 0

final_df <- cbind(combined_ed[,-c(54:62)],num_df)
  
####OUTCOMES
final_df <- merge(final_df, admissions[,c("subject_id","hadm_id","admittime","dischtime","hospital_expire_flag",
                                        "insurance","marital_status")], 
                  by = c("subject_id","hadm_id"), all.x= TRUE)

#drg severity/mortality
final_df$drg_severity[is.na(final_df$drg_severity)] <- 0
final_df$drg_mortality[is.na(final_df$drg_mortality)] <- 0

#insurance, marital status
final_df$insurance[is.na(final_df$insurance)] <- "Other"
final_df$marital_status[is.na(final_df$marital_status)] <- ""

#remove missing covariate rows-- keep only complete cases
final_df <- final_df[complete.cases(final_df[,c(6:8,10:17,19:57,67:68)]), ]

#create mortality flags
final_df$er_death <- ifelse(final_df$disposition == "EXPIRED",1,0)
final_df$death <- ifelse(final_df$hospital_expire_flag == 1,1,0)
final_df$death[is.na(final_df$death)] <- 0

#create readmission count
admissions$one <- 1
readmit_num <- aggregate(admissions[,c("one")],by = list(admissions$subject_id), FUN = sum)
colnames(readmit_num) <- c("subject_id", "readmit_num")
final_df <- merge(final_df, readmit_num, by = c("subject_id"), all.x= TRUE)
final_df$readmit_num[is.na(final_df$readmit_num)] <- 0

#create LOS (for admissions excluding death)
final_admission_df <- subset(final_df,disposition == "ADMITTED")
final_admission_df <- subset(final_admission_df,hospital_expire_flag != 1)
final_admission_df$los <- as.Date(final_admission_df$dischtime) - as.Date(final_admission_df$admittime)

#convert treatments to binary flags
final_df$rad_flag <- ifelse(final_df$rad_num > 0, 1, 0) 
final_df$x_ray_flag <- ifelse(final_df$x_ray_num > 0, 1, 0)
final_df$ct_flag <- ifelse(final_df$ct_num > 0, 1, 0)
final_df$mri_flag <- ifelse(final_df$mri_num > 0, 1, 0)
final_df$us_flag <- ifelse(final_df$us_num > 0, 1, 0)
final_df$fluoro_flag <- ifelse(final_df$fluoro_num > 0, 1, 0)

#subset renal patients
rf_diag_ids <- unique(diag$subject_id[(grepl("renal",diag$icd_title,ignore.case = TRUE)|
                                         grepl("kidney",diag$icd_title,ignore.case = TRUE))&
                                        !(grepl("adrenal",diag$icd_title,ignore.case = TRUE))])
#use drg descriptions
rf_drg_ids <- unique(drg$subject_id[(grepl("renal",drg$description,ignore.case = TRUE)|
                                grepl("kidney",drg$description,ignore.case = TRUE))&
                                !(grepl("adrenal",drg$description,ignore.case = TRUE))])

rf_ids <- unique(c(rf_diag_ids,rf_drg_ids))

final_rf <- subset(final_df, subject_id %in% rf_ids)

#create contrast ct and mr flags
final_rf$contr_ct <- ifelse(final_rf$subject_id %in% ct_contr_id,1,0)
final_rf$contr_mr <- ifelse(final_rf$subject_id %in% mr_contr_id,1,0)

#create LOS (for admissions excluding death)
final_admission_rf <- subset(final_rf,disposition == "ADMITTED")
final_admission_rf <- subset(final_admission_rf,hospital_expire_flag != 1)
final_admission_rf$los <- as.Date(final_admission_rf$dischtime) - as.Date(final_admission_rf$admittime)

#aggregate admissions to patient level by taking last admission
final_rf$indate <- as.POSIXlt(final_rf$intime)
final_rf_id <- aggregate(final_rf[,c("indate")],by = list(final_rf$subject_id), FUN = max)
colnames(final_rf_id) <- c("subject_id", "indate")
final_rf_ben <- merge(final_rf_id, final_rf, by = c("subject_id","indate"),all.x=TRUE, all.y = FALSE)

#filter repeats manually
final_rf_ben_dup_ids <- final_rf_ben$subject_id[duplicated(final_rf_ben$subject_id)]
final_rf_ben_dups <- final_rf_ben[final_rf_ben$subject_id %in% final_rf_ben_dup_ids,c("subject_id","stay_id","intime")]
final_rf_ben_dups$intime1 <- as.POSIXlt(final_rf_ben_dups$intime)
final_rf_ben_dups_rm_id <- aggregate(final_rf_ben_dups[,c("intime1")],by = list(final_rf_ben_dups$subject_id), FUN = min)
colnames(final_rf_ben_dups_rm_id) <- c("subject_id", "intime1")
final_rf_ben_dups_rm <- merge(final_rf_ben_dups_rm_id, final_rf_ben_dups, by = c("subject_id","intime1"),all.x=TRUE, all.y = FALSE)

final_rf_ben1 <- final_rf_ben[!(final_rf_ben$stay_id %in% c(final_rf_ben_dups_rm$stay_id,"32263202")),]

#####Analysis
#tmle analysis
sl_library <- c("SL.ranger","SL.glm","SL.ksvm","SL.xgboost","SL.biglasso")

readmit_ct_tmle <- tmle(Y=final_rf_ben1$readmit_num,
                        A=final_rf_ben1$contr_ct,
                        W=final_rf_ben1[,c(7:9,11:18,20:58,68:69)],
                        Q.SL.library = sl_library,
                        g.SL.library = sl_library,
                        family = "gaussian",verbose=TRUE)

save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle_rf1.RData")

readmit_mri_tmle <- tmle(Y=final_rf_ben1$readmit_num,
                         A=final_rf_ben1$contr_mr,
                         W=final_rf_ben1[,c(7:9,11:18,20:58,68:69)],
                         Q.SL.library = sl_library,
                         g.SL.library = sl_library,
                         family = "gaussian",verbose=TRUE)

save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle_rf1.RData")

mort_ct_tmle <- tmle(Y=final_rf_ben1$death,
                     A=final_rf_ben1$contr_ct,
                     W=final_rf_ben1[,c(7:9,11:18,20:58,68:69)],
                     Q.SL.library = sl_library,
                     g.SL.library = sl_library,
                     family = "binomial",verbose=TRUE)

save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle_rf1.RData")

mort_mri_tmle <- tmle(Y=final_rf_ben1$death,
                      A=final_rf_ben1$contr_mr,
                      W=final_rf_ben1[,c(7:9,11:18,20:58,68:69)],
                      Q.SL.library = sl_library,
                      g.SL.library = sl_library,
                      family = "binomial",
                      verbose=TRUE)

save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle_rf1.RData")


##OLD TMLE
# mort_ct_tmle <- tmle(Y=final_df$death,
#                       A=final_df$ct_flag,
#                       W=final_df[,c(6:8,10:17,19:57,67:68)],
#                       Q.SL.library = sl_library,
#                       g.SL.library = sl_library,
#                       family = "binomial",
#                       V=5,verbose=TRUE)
# 
# save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle.RData")
# 
# mort_mri_tmle <- tmle(Y=final_df$death,
#                   A=final_df$mri_flag,
#                   W=final_df[,c(6:8,10:17,19:57,67:68)],
#                   Q.SL.library = sl_library,
#                   g.SL.library = sl_library,
#                   family = "binomial",
#                   V=5,verbose=TRUE)
# 
# save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle.RData")
# 
# mort_us_tmle <- tmle(Y=final_df$death,
#                       A=final_df$us_flag,
#                       W=final_df[,c(6:8,10:17,19:57,67:68)],
#                       Q.SL.library = sl_library,
#                       g.SL.library = sl_library,
#                       family = "binomial",
#                       V=5,verbose=TRUE)
# 
# save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle.RData")
# 
# readmit_mri_tmle <- tmle(Y=final_df$readmit_num,
#                      A=final_df$mri_flag,
#                      W=final_df[,c(6:8,10:17,19:57,67:68)],
#                      Q.SL.library = sl_library,
#                      g.SL.library = sl_library,
#                      family = "gaussian",
#                      V=5,verbose=TRUE)
# 
# save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle.RData")
# 
# readmit_ct_tmle <- tmle(Y=final_df$readmit_num,
#                          A=final_df$ct_flag,
#                          W=final_df[,c(6:8,10:17,19:57,67:68)],
#                          Q.SL.library = sl_library,
#                          g.SL.library = sl_library,
#                          family = "gaussian",
#                          V=5,verbose=TRUE)
# 
# save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle.RData")
# 
# readmit_us_tmle <- tmle(Y=final_df$readmit_num,
#                         A=final_df$us_flag,
#                         W=final_df[,c(6:8,10:17,19:57,67:68)],
#                         Q.SL.library = sl_library,
#                         g.SL.library = sl_library,
#                         family = "gaussian",
#                         V=5,verbose=TRUE)
# 
# save.image("P:/pacce/s4k/MIMIC4/mimic4_tmle.RData")

