#try with removing outlies
#drop outliners
Q <- quantile(no_na_rep$log_mut, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(no_na_rep$log_mut)

up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

log_eliminate <- subset(no_na_rep, no_na_rep$log_mut > (Q[1] - 1.5*iqr) & no_na_rep$log_mut < (Q[2]+1.5*iqr))

regmult_c <- lm(log_mut~factor(gender)+age+factor(alcohol_history),
                data = log_eliminate)

summary(regmult_c) #females,males,age are significant and R^2 = 0.16

nobs(regmult_c)#6752

res_a <- resid(regmult_c) 

plot(fitted(regmult_c), res_a)
abline(0,0) #res plot looks random = normal

qqnorm(res_a)
qqline(res_a) #follows the line pretty darn closely. 

#try with more variables
regmult_d <- lm(log_mut~factor(gender)+age+factor(alcohol_history)+factor(race),
                data = log_eliminate)

summary(regmult_d) #females,males,age are significant and R^2 = 0.16

nobs(regmult_d)#6752

res <- resid(regmult_d) 

plot(fitted(regmult_d), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) 

#try with even more variables
regmult_e <- lm(log_mut~factor(gender)+age+factor(alcohol_history)+factor(race)+factor(ethnicity),
                data = log_eliminate)

summary(regmult_e) #females,males,age are significant and R^2 = 0.162

nobs(regmult_e)#6752

res <- resid(regmult_e) 

plot(fitted(regmult_e), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) 

#try with even more variables
regmult_f <- lm(log_mut~factor(gender)+age+factor(alcohol_history)+factor(race)+
                  factor(ethnicity)+factor(alcohol_history),
                data = log_eliminate)

summary(regmult_f) #females,males,age are significant and R^2 = 0.162

nobs(regmult_f)#6752

res <- resid(regmult_f) 

plot(fitted(regmult_f), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) 
#same result as last two

#try with even more variables
regmult_g <- lm(log_mut~factor(gender)+age+factor(alcohol_history)+factor(race)+
                  factor(ethnicity)+factor(alcohol_history)+years_smoked+cigarettes_per_day,
                data = log_eliminate)

summary(regmult_g) #females,alcohol_history=not reported, years_smoked are significant and R^2 = 0.22

nobs(regmult_g)#559

res <- resid(regmult_g) 

plot(fitted(regmult_g), res)
abline(0,0) #res plot looks very random = normal

qqnorm(res)
qqline(res)  #at the ends this does not follow the line


summary(log_eliminate$total_mut)
var(log_eliminate$total_mut)
sd(log_eliminate$total_mut)
range(log_eliminate$total_mut)

