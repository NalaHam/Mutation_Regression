library(moments)
library(dplyr)


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

nobs(reg_luad)#473

res<- resid(reg_luad) 

plot(fitted(reg_luad), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird data point 











