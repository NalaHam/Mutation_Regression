
unique(log_eliminate$Project)

#GBM
gbm <- filter(no_na_rep, Project == "TCGA-GBM")

hist(gbm$log_mut)

Q <- quantile(gbm$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(gbm$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

gbm <- subset(gbm, gbm$log_mut > (Q[1] - 1.5*iqr) & gbm$log_mut < (Q[2]+1.5*iqr))

#regression
unique(gbm$alcohol_history) #no data 
unique(gbm$race)
unique(gbm$ethnicity) 
unique(gbm$cancer_level) #no data 
unique(gbm$years_smoked) #no data 
unique(gbm$cigarettes_per_day) #no data 


reg_gbm <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity),
               data = gbm)

summary(reg_gbm) #female and age are significant and R^2 = 0.12

nobs(reg_gbm)#356

res<- resid(reg_gbm) 

plot(fitted(reg_gbm), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part



#LUAD
luad <- filter(full_mut_count, Project == "TCGA-LUAD")

luad <- luad[complete.cases(luad[ , 3:15]),]

skewness(luad$total_mut, na.rm = TRUE)


luad$log_mut <- log10(luad$total_mut) #transformed the mutations into log10() of them, can be viewed as percentages

hist(luad$log_mut)

Q <- quantile(luad$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(luad$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

luad <- subset(luad, luad$log_mut > (Q[1] - 1.5*iqr) & luad$log_mut < (Q[2]+1.5*iqr))

#regression
unique(luad$alcohol_history) #no data 
unique(luad$race)
unique(luad$ethnicity) 
unique(luad$cancer_level) 
unique(luad$years_smoked) 
unique(luad$cigarettes_per_day)


reg_luad <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
              data = luad)

summary(reg_luad) #female is significant and R^2 = 0.05

nobs(reg_luad)#475

res<- resid(reg_luad) 

plot(fitted(reg_luad), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird data point 



#LUSC
lusc <- filter(full_mut_count, Project == "TCGA-LUSC")

lusc <- lusc[complete.cases(lusc[ , 3:15]),]

skewness(lusc$total_mut, na.rm = TRUE)


lusc$log_mut <- log10(lusc$total_mut) #transformed the mutations into log10() of them, can be viewed as percentages

hist(lusc$log_mut)

Q <- quantile(lusc$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(lusc$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

lusc <- subset(lusc, lusc$log_mut > (Q[1] - 1.5*iqr) & lusc$log_mut < (Q[2]+1.5*iqr))
#regression
unique(lusc$alcohol_history) #no data 
unique(lusc$race)
unique(lusc$ethnicity) 
unique(lusc$cancer_level) 
unique(lusc$years_smoked) 
unique(lusc$cigarettes_per_day)


reg_lusc <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = lusc)

summary(reg_lusc) #female is significant and R^2 = 0.06

nobs(reg_lusc)#450

res<- resid(reg_lusc) 

plot(fitted(reg_lusc), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line 




#BLCA
blca <- filter(log_eliminate, Project == "TCGA-BLCA")
hist(blca$log_mut)


#regression
unique(blca$alcohol_history) #no data 
unique(blca$race)
unique(blca$ethnicity) 
unique(blca$cancer_level) 
unique(blca$years_smoked) #no data
unique(blca$cigarettes_per_day)


reg_blca <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = blca)

summary(reg_blca) #female and white are significant and R^2 = 0.09

nobs(reg_blca)#403

res<- resid(reg_blca) 

plot(fitted(reg_blca), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part 




#PAAD
paad <- filter(log_eliminate, Project == "TCGA-PAAD")
hist(paad$log_mut) #not normal

paad <- filter(full_mut_count, Project == "TCGA-PAAD")

paad <- paad[complete.cases(paad[ , 3:15]),]

skewness(paad$total_mut, na.rm = TRUE)


paad$log_mut <- log10(paad$total_mut) #transformed the mutations into log10() of them, can be viewed as percentages

hist(paad$log_mut)#still not normal

Q <- quantile(paad$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(paad$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

paad <- subset(paad, paad$log_mut > (Q[1] - 1.5*iqr) & paad$log_mut < (Q[2]+1.5*iqr))

hist(paad$log_mut) #great

#regression
unique(paad$alcohol_history) 
unique(paad$race)
unique(paad$ethnicity) 
unique(paad$cancer_level) 
unique(paad$years_smoked) #not alot of data
unique(paad$cigarettes_per_day)#not alot of data


reg_paad <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = paad)

summary(reg_paad) #female and male are significant and R^2 = 0.1

nobs(reg_paad)#152

res<- resid(reg_paad) 

plot(fitted(reg_paad), res)
abline(0,0) #res plot looks very random = normal

qqnorm(res)
qqline(res) #follows the line mostly



#KIRP
kirp <- filter(log_eliminate, Project == "TCGA-KIRP")
hist(kirp$log_mut)

kirp <- filter(full_mut_count, Project == "TCGA-KIRP")

kirp <- kirp[complete.cases(kirp[ , 3:15]),]

skewness(kirp$total_mut, na.rm = TRUE)


kirp$log_mut <- log10(kirp$total_mut) #transformed the mutations into log10() of them, can be viewed as percentages

hist(kirp$log_mut)#still not normal

Q <- quantile(kirp$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(kirp$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

kirp <- subset(kirp, kirp$log_mut > (Q[1] - 1.5*iqr) & kirp$log_mut < (Q[2]+1.5*iqr))

hist(kirp$log_mut) #better

#regression
unique(kirp$alcohol_history) #no data 
unique(kirp$race)
unique(kirp$ethnicity) 
unique(kirp$cancer_level) 
unique(kirp$years_smoked) #not a lot of data
unique(kirp$cigarettes_per_day)#not alot of data


reg_kirp <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = kirp)

summary(reg_kirp) #female, male, and age are significant and R^2 = 0.3

nobs(reg_kirp)#241

res<- resid(reg_kirp) 

plot(fitted(reg_kirp), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line 




#LIHC
lihc <- filter(log_eliminate, Project == "TCGA-LIHC")
hist(lihc$log_mut)
lihc <- filter(full_mut_count, Project == "TCGA-LIHC")

lihc <- lihc[complete.cases(lihc[ , 3:15]),]

skewness(lihc$total_mut, na.rm = TRUE)


lihc$log_mut <- log10(lihc$total_mut) #transformed the mutations into log10() of them, can be viewed as percentages

hist(lihc$log_mut)#still not normal

Q <- quantile(lihc$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(lihc$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

lihc <- subset(lihc, lihc$log_mut > (Q[1] - 1.5*iqr) & lihc$log_mut < (Q[2]+1.5*iqr))

hist(lihc$log_mut) #better
#regression
unique(lihc$alcohol_history) #no data 
unique(lihc$race)
unique(lihc$ethnicity) 
unique(lihc$cancer_level) 
unique(lihc$years_smoked) #no data
unique(lihc$cigarettes_per_day) #no data


reg_lihc <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = lihc)

summary(reg_lihc) #female and age are significant and R^2 = 0.22

nobs(reg_lihc)#324

res<- resid(reg_lihc) 

plot(fitted(reg_lihc), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line




#SARC
SARC <- filter(no_na_rep, Project == "TCGA-SARC")
hist(SARC$log_mut)

Q <- quantile(SARC$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(SARC$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

SARC <- subset(SARC, SARC$log_mut > (Q[1] - 1.5*iqr) & SARC$log_mut < (Q[2]+1.5*iqr))

hist(SARC$log_mut) #better


#regression
unique(SARC$alcohol_history) #no data 
unique(SARC$race)
unique(SARC$ethnicity) 
unique(SARC$cancer_level) #none
unique(SARC$years_smoked) #none
unique(SARC$cigarettes_per_day) #none


reg_SARC <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity),
               data = SARC)

summary(reg_SARC) #female and age are significant and R^2 = 0.2

nobs(reg_SARC)#222

res<- resid(reg_SARC) 

plot(fitted(reg_SARC), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line 




#THYM
THYM <- filter(no_na_rep, Project == "TCGA-THYM")
hist(THYM$log_mut)

Q <- quantile(THYM$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(THYM$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

THYM <- subset(THYM, THYM$log_mut > (Q[1] - 1.5*iqr) & THYM$log_mut < (Q[2]+1.5*iqr))

hist(THYM$log_mut) #better


#regression
unique(THYM$alcohol_history) #no data 
unique(THYM$race)
unique(THYM$ethnicity) 
unique(THYM$cancer_level) 
unique(THYM$years_smoked) #none
unique(THYM$cigarettes_per_day) #none


reg_THYM <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = THYM)

summary(reg_THYM) #age is significant and R^2 = 0.3

nobs(reg_THYM)#117

res<- resid(reg_THYM) 

plot(fitted(reg_THYM), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird data point 




#MESO doesn't have over 100 samples
MESO <- filter(no_na_rep, Project == "TCGA-MESO")
hist(MESO$log_mut)

Q <- quantile(MESO$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(MESO$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

MESO <- subset(MESO, MESO$log_mut > (Q[1] - 1.5*iqr) & MESO$log_mut < (Q[2]+1.5*iqr))

hist(MESO$log_mut) # better?

#regression
unique(MESO$alcohol_history) #no data 
unique(MESO$race)
unique(MESO$ethnicity) #only white
unique(MESO$cancer_level) 
unique(MESO$years_smoked) #none
unique(MESO$cigarettes_per_day) #none


reg_MESO <- lm(log_mut~factor(gender)+age+factor(race)+factor(cancer_level),
               data = MESO)

summary(reg_MESO) #female and age are significant and R^2 = 0.3

nobs(reg_MESO)#75

res<- resid(reg_MESO) 

plot(fitted(reg_MESO), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part




#COAD
COAD <- filter(no_na_rep, Project == "TCGA-COAD")
hist(COAD$log_mut)

Q <- quantile(COAD$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(COAD$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

COAD <- subset(COAD, COAD$log_mut > (Q[1] - 1.5*iqr) & COAD$log_mut < (Q[2]+1.5*iqr))

hist(COAD$log_mut) #much better!

#regression
unique(COAD$alcohol_history) #no data 
unique(COAD$race)
unique(COAD$ethnicity) 
unique(COAD$cancer_level) 
unique(COAD$years_smoked) #none
unique(COAD$cigarettes_per_day) #none


reg_COAD <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = COAD)

summary(reg_COAD) #female and age are significant and R^2 = 0.2

nobs(reg_COAD)#319

res<- resid(reg_COAD) 

plot(fitted(reg_COAD), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but weird lower end




#STAD
STAD <- filter(no_na_rep, Project == "TCGA-STAD")
hist(STAD$log_mut)

Q <- quantile(STAD$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(STAD$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

STAD <- subset(STAD, STAD$log_mut > (Q[1] - 1.5*iqr) & STAD$log_mut < (Q[2]+1.5*iqr))

hist(STAD$log_mut) # better

#regression
unique(STAD$alcohol_history) #no data 
unique(STAD$race)
unique(STAD$ethnicity) 
unique(STAD$cancer_level) 
unique(STAD$years_smoked) #none
unique(STAD$cigarettes_per_day)#none


reg_STAD <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = STAD)

summary(reg_STAD) #female is significant and R^2 = 0.1

nobs(reg_STAD)#366

res<- resid(reg_STAD) 

plot(fitted(reg_STAD), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but right end is bad!




#SKCM 
SKCM <- filter(no_na_rep, Project == "TCGA-SKCM")
hist(SKCM$log_mut)

Q <- quantile(SKCM$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(SKCM$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

SKCM <- subset(SKCM, SKCM$log_mut > (Q[1] - 1.5*iqr) & SKCM$log_mut < (Q[2]+1.5*iqr))

hist(SKCM$log_mut) # better

#regression
unique(SKCM$alcohol_history) #no data 
unique(SKCM$race) #only has asian or white
unique(SKCM$ethnicity) 
unique(SKCM$cancer_level) 
unique(SKCM$years_smoked) #none
unique(SKCM$cigarettes_per_day)#none


reg_SKCM <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = SKCM)

summary(reg_SKCM) #female and white are significant and R^2 = 0.15

nobs(reg_SKCM)#412

res<- resid(reg_SKCM) 

plot(fitted(reg_SKCM), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird data point 




#CHOL
CHOL <- filter(no_na_rep, Project == "TCGA-CHOL")
hist(CHOL$log_mut)

Q <- quantile(CHOL$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CHOL$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

CHOL <- subset(CHOL, CHOL$log_mut > (Q[1] - 1.5*iqr) & CHOL$log_mut < (Q[2]+1.5*iqr))

hist(CHOL$log_mut) # better
#regression
unique(CHOL$alcohol_history) #no data 
unique(CHOL$race)
unique(CHOL$ethnicity) 
unique(CHOL$cancer_level) 
unique(CHOL$years_smoked) #none
unique(CHOL$cigarettes_per_day)#none


reg_CHOL <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = CHOL)

summary(reg_CHOL) #female is significant and R^2 = 0.3 and 0.09

nobs(reg_CHOL)#42

res<- resid(reg_CHOL) 

plot(fitted(reg_CHOL), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part




#KIRC
KIRC <- filter(no_na_rep, Project == "TCGA-KIRC")
hist(KIRC$log_mut)

Q <- quantile(KIRC$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(KIRC$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

KIRC <- subset(KIRC, KIRC$log_mut > (Q[1] - 1.5*iqr) & KIRC$log_mut < (Q[2]+1.5*iqr))

hist(KIRC$log_mut) # much better
#regression
unique(KIRC$alcohol_history) #no data 
unique(KIRC$race)
unique(KIRC$ethnicity) 
unique(KIRC$cancer_level) 
unique(KIRC$years_smoked) #not that much data
unique(KIRC$cigarettes_per_day) #not that much data


reg_KIRC <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = KIRC)

summary(reg_KIRC) #female and age are significant and R^2 = 0.2

nobs(reg_KIRC)#313

res<- resid(reg_KIRC) 

plot(fitted(reg_KIRC), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird at the ends




#THCA
THCA <- filter(no_na_rep, Project == "TCGA-THCA")
hist(THCA$log_mut)

Q <- quantile(THCA$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(THCA$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

THCA <- subset(THCA, THCA$log_mut > (Q[1] - 1.5*iqr) & THCA$log_mut < (Q[2]+1.5*iqr))

hist(THCA$log_mut) # better?
#regression
unique(THCA$alcohol_history) #no data 
unique(THCA$race)
unique(THCA$ethnicity) 
unique(THCA$cancer_level) 
unique(THCA$years_smoked) #none
unique(THCA$cigarettes_per_day)#none


reg_THCA <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = THCA)

summary(reg_THCA) #female, age, and stage IVA are significant and R^2 = 0.09

nobs(reg_THCA)#486

res<- resid(reg_THCA) 

plot(fitted(reg_THCA), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but weird on ends




#HNSC
HNSC <- filter(no_na_rep, Project == "TCGA-HNSC")
hist(HNSC$log_mut)

Q <- quantile(HNSC$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(HNSC$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

HNSC <- subset(HNSC, HNSC$log_mut > (Q[1] - 1.5*iqr) & HNSC$log_mut < (Q[2]+1.5*iqr))

hist(HNSC$log_mut) #much better

#regression
unique(HNSC$alcohol_history) #has data
unique(HNSC$race)
unique(HNSC$ethnicity) 
unique(HNSC$cancer_level) 
unique(HNSC$years_smoked) 
unique(HNSC$cigarettes_per_day)


reg_HNSC <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = HNSC)

summary(reg_HNSC) #female and age are significant and R^2 = 0.05

nobs(reg_HNSC)#475

res<- resid(reg_HNSC) 

plot(fitted(reg_HNSC), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird data point 




#LAML
LAML <- filter(no_na_rep, Project == "TCGA-LAML")
hist(LAML$log_mut)

Q <- quantile(LAML$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(LAML$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

LAML <- subset(LAML, LAML$log_mut > (Q[1] - 1.5*iqr) & LAML$log_mut < (Q[2]+1.5*iqr))

hist(LAML$log_mut) # still very weird

#regression
unique(LAML$alcohol_history) #no data 
unique(LAML$race)
unique(LAML$ethnicity) 
unique(LAML$cancer_level) #none
unique(LAML$years_smoked) #none
unique(LAML$cigarettes_per_day)#none


reg_LAML <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity),
               data = LAML)

summary(reg_LAML) #female and age are significant and R^2 = 0.05

nobs(reg_LAML)#119

res<- resid(reg_LAML) 

plot(fitted(reg_LAML), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #bad!




#READ
READ <- filter(no_na_rep, Project == "TCGA-READ")
hist(READ$log_mut)

Q <- quantile(READ$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(READ$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

READ <- subset(READ, READ$log_mut > (Q[1] - 1.5*iqr) & READ$log_mut < (Q[2]+1.5*iqr))

hist(READ$log_mut) # better

#regression
unique(READ$alcohol_history) #no data 
unique(READ$race)
unique(READ$ethnicity) 
unique(READ$cancer_level) 
unique(READ$years_smoked) #none
unique(READ$cigarettes_per_day) #none


reg_READ <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = READ)

summary(reg_READ) #female and age are significant and R^2 = 0.21 and 0.08

nobs(reg_READ)#120

res<- resid(reg_READ) 

plot(fitted(reg_READ), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part 




#LGG
LGG <- filter(no_na_rep, Project == "TCGA-LGG")
hist(LGG$log_mut)

Q <- quantile(LGG$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(LGG$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

LGG <- subset(LGG, LGG$log_mut > (Q[1] - 1.5*iqr) & LGG$log_mut < (Q[2]+1.5*iqr))

hist(LGG$log_mut) #much better

#regression
unique(LGG$alcohol_history) #no data 
unique(LGG$race)
unique(LGG$ethnicity) 
unique(LGG$cancer_level) #none
unique(LGG$years_smoked) #none
unique(LGG$cigarettes_per_day) #none


reg_LGG <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity),
               data = LGG)

summary(reg_LGG) #female and age are significant and R^2 = 0.33

nobs(reg_LGG)#485

res<- resid(reg_LGG) 

plot(fitted(reg_LGG), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part




#DLBC
DLBC <- filter(no_na_rep, Project == "TCGA-DLBC")
hist(DLBC$log_mut)

Q <- quantile(DLBC$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(DLBC$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

DLBC <- subset(DLBC, DLBC$log_mut > (Q[1] - 1.5*iqr) & DLBC$log_mut < (Q[2]+1.5*iqr))

hist(DLBC$log_mut) # better

#regression
unique(DLBC$alcohol_history) #no data 
unique(DLBC$race)
unique(DLBC$ethnicity) 
unique(DLBC$cancer_level) 
unique(DLBC$years_smoked) #none
unique(DLBC$cigarettes_per_day)#none


reg_DLBC <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = DLBC)

summary(reg_DLBC) #female is significant and R^2 = 0.2 and -0.03

nobs(reg_DLBC)#40

res<- resid(reg_DLBC) 

plot(fitted(reg_DLBC), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but kinda weird




#KICH
KICH <- filter(no_na_rep, Project == "TCGA-KICH")
hist(KICH$log_mut)

Q <- quantile(KICH$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(KICH$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

KICH <- subset(KICH, KICH$log_mut > (Q[1] - 1.5*iqr) & KICH$log_mut < (Q[2]+1.5*iqr))

hist(KICH$log_mut) # better

#regression
unique(KICH$alcohol_history) #no data 
unique(KICH$race)
unique(KICH$ethnicity) 
unique(KICH$cancer_level) 
unique(KICH$years_smoked) 
unique(KICH$cigarettes_per_day)


reg_KICH <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = KICH)

summary(reg_KICH) #female and age are significant and R^2 = 0.09

nobs(reg_KICH)#475

res<- resid(reg_KICH) 

plot(fitted(reg_KICH), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird data point 




#ACC
ACC <- filter(no_na_rep, Project == "TCGA-ACC")
hist(ACC$log_mut)

Q <- quantile(ACC$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(ACC$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

ACC <- subset(ACC, ACC$log_mut > (Q[1] - 1.5*iqr) & ACC$log_mut < (Q[2]+1.5*iqr))

hist(ACC$log_mut) # better
#regression
unique(ACC$alcohol_history) #no data 
unique(ACC$race)
unique(ACC$ethnicity) 
unique(ACC$cancer_level) 
unique(ACC$years_smoked) #none
unique(ACC$cigarettes_per_day) #none


reg_ACC <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = ACC)

summary(reg_ACC) #female, age, and stage 4 are significant and R^2 = 0.3

nobs(reg_ACC)#85

res<- resid(reg_ACC) 

plot(fitted(reg_ACC), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird data point 




#PCPG
PCPG <- filter(no_na_rep, Project == "TCGA-PCPG")
hist(PCPG$log_mut)

Q <- quantile(PCPG$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(PCPG$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

PCPG <- subset(PCPG, PCPG$log_mut > (Q[1] - 1.5*iqr) & PCPG$log_mut < (Q[2]+1.5*iqr))

hist(PCPG$log_mut) # better
#regression
unique(PCPG$alcohol_history) #no data 
unique(PCPG$race)
unique(PCPG$ethnicity) 
unique(PCPG$cancer_level) #none
unique(PCPG$years_smoked) #none
unique(PCPG$cigarettes_per_day)#none


reg_PCPG <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity),
               data = PCPG)

summary(reg_PCPG) #female and age are significant and R^2 = 0.2

nobs(reg_PCPG)#175

res<- resid(reg_PCPG) 

plot(fitted(reg_PCPG), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but a little weird on the ends




#UVM
UVM <- filter(no_na_rep, Project == "TCGA-UVM")
hist(UVM$log_mut)

Q <- quantile(UVM$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(UVM$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

UVM <- subset(UVM, UVM$log_mut > (Q[1] - 1.5*iqr) & UVM$log_mut < (Q[2]+1.5*iqr))

hist(UVM$log_mut) # better
#regression
unique(UVM$alcohol_history) #no data 
unique(UVM$race)
unique(UVM$ethnicity) 
unique(UVM$cancer_level) 
unique(UVM$years_smoked) #none
unique(UVM$cigarettes_per_day)#none


reg_UVM <- lm(log_mut~factor(gender)+age+factor(race)+factor(ethnicity)+factor(cancer_level),
               data = UVM)

summary(reg_UVM) #female and age are significant and R^2 = 0.17 and 0.04

nobs(reg_UVM)#76

res<- resid(reg_UVM) 

plot(fitted(reg_UVM), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but weird at the ends











