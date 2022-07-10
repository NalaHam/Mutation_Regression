


no_na_rep$ln_mut <- log(no_na_rep$total_mut)

hist(no_na_rep$ln_mut)

#can skip (hist looks very normal)
Q <- quantile(luad$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(luad$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

luad <- subset(luad, luad$log_mut > (Q[1] - 1.5*iqr) & luad$log_mut < (Q[2]+1.5*iqr))
#end of skip

ln_reg_non_rep <- lm(ln_mut~factor(gender)+age+factor(alcohol_history)+factor(race)+
                       factor(ethnicity),
                data = no_na_rep)

summary(ln_reg_non_rep)

nobs(ln_reg_non_rep) #6908

res<- resid(ln_reg_non_rep) 

plot(fitted(ln_reg_non_rep), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #follows the line for the most part but a little weird on ends 

#understanding the use of log() on the coefficients 

(exp(coef(ln_reg_non_rep)["factor(gender)male"]) - 1 ) *100 # 38.83% increase in mutations if you are a male 


(exp(coef(ln_reg_non_rep)["age"]) - 1 ) *100 # 3.74% increase in mutations for each year older someone is  


range(no_na_rep$age, na.rm = TRUE) #14 to 90













