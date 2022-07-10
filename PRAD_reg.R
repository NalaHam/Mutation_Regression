#PRAD regression

#organize full_mut_count df 
full_mut_count$total_mut <- rowSums(full_mut_count[,3:14], na.rm = FALSE)

full_mut_count <- full_mut_count[, c(15,2,3,4,5,6,7,8,9,10,11,12,13,14,24,16,17,18,19,20,21,22,23)]

write.csv(full_mut_count, file = "all_mutation_df.csv")

prostate_cancer <- filter(full_mut_count, Project == "TCGA-PRAD")

prostate_cancer <- prostate_cancer[complete.cases(prostate_cancer[ , 3:14]),] #490 people

#regression
unique(prostate_cancer$alcohol_history) #no data on years_smoked, cigs a day, or alcohol_history

reg_prad <- lm(total_mut~age+factor(race)+factor(ethnicity),
                data = prostate_cancer)

summary(reg_prad) #age are significant and R^2 = 0.002

nobs(reg_prad)#490

res<- resid(reg_prad) 

plot(fitted(reg_prad), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but one weird data point 

# ^not great stuff, try with removing outlies

Q <- quantile(prostate_cancer$total_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(prostate_cancer$total_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

prad_eliminate <- subset(prostate_cancer, prostate_cancer$total_mut > (Q[1] - 1.5*iqr)
                         & prostate_cancer$total_mut < (Q[2]+1.5*iqr))

#redo regression
reg_prad <- lm(total_mut~age+factor(race)+factor(ethnicity),
               data = prad_eliminate)

summary(reg_prad) #age is significant and R^2 = 0.07

nobs(reg_prad)#458

res<- resid(reg_prad) 

plot(fitted(reg_prad), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part 
















