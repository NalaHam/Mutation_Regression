library(dplyr)

#organize data
rep_mut <- filter(full_mut_count, Project %in% c("TCGA-ACC", "TCGA-BLCA",
                                                    "TCGA-COAD", "TCGA-GBM", "TCGA-HNSC", "TCGA-KIRC", 
                                                    "TCGA-LUAD", "TCGA-CHOL", "TCGA-KICH", "TCGA-KIRP", 
                                                    "TCGA-LIHC", "TCGA-LUSC", "TCGA-MESO", "TCGA-READ", 
                                                    "TCGA-STAD", "TCGA-THCA", "TCGA-UVM", "TCGA-DLBC", 
                                                    "TCGA-LGG", "TCGA-SARC","TCGA-THYM", "TCGA-PCPG", "TCGA-LAML",
                                                    "TCGA-PAAD", "TCGA-SKCM"))


rep_mut$total_mut <- rowSums(rep_mut[,3:14], na.rm = FALSE)

rep_mut <- rep_mut[, c(15,2,3,4,5,6,7,8,9,10,11,12,13,14,24,16,17,18,19,20,21,22,23)]

write.csv(rep_mut, file = "reproductive_mutation_df.csv")


no_na_rep <- rep_mut[complete.cases(rep_mut[ , 3:14]),]

#see is anything is a factor
is.factor(no_na_rep$race)
is.factor(no_na_rep$cancer_level)
is.factor(no_na_rep$Project)
is.factor(no_na_rep$gender)
is.factor(no_na_rep$ethnicity)

#regression on uncontrolled data
regmult <- lm(total_mut~factor(gender)+age,
              data = no_na_rep)

summary(regmult)

mean(no_na_rep$total_mut) #234
range(no_na_rep$total_mut) #1 and 22065. The upper part is extremely high for the mean

sum(no_na_rep$total_mut >= 2000) #108

#drop outliners
Q <- quantile(no_na_rep$total_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(no_na_rep$total_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

eliminated <- subset(no_na_rep, no_na_rep$total_mut > (Q[1] - 1.5*iqr) & no_na_rep$total_mut < (Q[2]+1.5*iqr))


#testing gender
regmult_1 <- lm(total_mut~factor(gender),
              data = eliminated)

summary(regmult_1)

mean(eliminated$total_mut) #94
range(eliminated$total_mut) #1 and 410

nobs(regmult_1) #6154

#adding age
regmult_2 <- lm(total_mut~factor(gender)+age,
                data = eliminated)

summary(regmult_2)

nobs(regmult_2) #6111
#most significance out of all tests. R^2 value is 0.106

#adding race
sum(is.na(no_na_rep$race))

regmult_3 <- lm(total_mut~factor(gender)+age+factor(race),
                data = eliminated)

summary(regmult_3)#6111

nobs(regmult_3)
#has no effect on the varibles influence, except females is no longer significant

#adding cancer stage
regmult_4 <- lm(total_mut~factor(gender)+age+factor(race)+factor(cancer_level),
                data = eliminated)

summary(regmult_4)

nobs(regmult_4) #4529
#not useful for pan cancer analysis, should do by cancer

#adding ethnicity
regmult_5 <- lm(total_mut~factor(gender)+age+factor(race)+factor(ethnicity),
                data = eliminated)

summary(regmult_5)

nobs(regmult_5) #6111
#same as race

#adding years smoked
is.numeric(no_na_rep$years_smoked)#check to see if numeric

regmult_6 <- lm(total_mut~factor(gender)+age+factor(race)+factor(ethnicity)+years_smoked,
                data = eliminated)

summary(regmult_6)

nobs(regmult_6)#485
#male results don't make sense and this has a low n for testing
# R^2 value at multi = 0.17 and adj= 0.159

regmult_a <- lm(total_mut~factor(gender)+age+years_smoked,
                data = eliminated)

summary(regmult_a)
nobs(regmult_a)#485
#results make some sense but there is still weirdness and n is small
#R^2 at 0.1329 and 0.1275

#adding cigarettes per day
is.numeric(no_na_rep$cigarettes_per_day)
regmult_7 <- lm(total_mut~factor(gender)+age+cigarettes_per_day,
                data = eliminated)

summary(regmult_7)

nobs(regmult_7)#1135
#age and cigs make sense but males has a negative effect. and R^2 is 0.04

#adding alcohol history
regmult_8 <- lm(total_mut~factor(gender)+age+factor(alcohol_history),
                data = eliminated)

summary(regmult_8)

nobs(regmult_8)#6111
#female, male, and age stay significant. R^2 is 0.106

#add everything but cancer stage

regmult_9 <- lm(total_mut~factor(gender)+age+factor(alcohol_history)+
                  cigarettes_per_day+years_smoked+factor(race)+factor(ethnicity),
                data = eliminated)

summary(regmult_9)

nobs(regmult_9)#452
#weirdest results but strongest R^2 value at around 0.25

#try without years_smoked
regmult_10 <- lm(total_mut~factor(gender)+age+factor(alcohol_history)+
                  cigarettes_per_day+factor(race)+factor(ethnicity),
                data = eliminated)

summary(regmult_10)

nobs(regmult_10)#1135
#weird results and R^2 is 0.14

#try without cigarettes_per_day
regmult_11 <- lm(total_mut~factor(gender)+age+factor(alcohol_history)+
                  factor(race)+factor(ethnicity),
                data = eliminated)

summary(regmult_11)

nobs(regmult_11)#6111
#okay results, no significance for females. R^2 values are 0.107






















