
#choose to not remove outliers 
regmult_b <- lm(log_mut~factor(gender)+age+factor(alcohol_history),
                data = no_na_rep)

summary(regmult_b) #female,male,and age are all significant. R^2 = 0.153

nobs(regmult_b)#6908

res <- resid(regmult_b) 

plot(fitted(regmult_b), res)
abline(0,0) #res plot looks random = normal

qqnorm(res)
qqline(res) #kinda follows the line but not at the ends

#try with more
regmult_ <- lm(log_mut~factor(gender)+age+factor(alcohol_history),
                data = no_na_rep)

summary(regmult_)

range(no_na_rep$total_mut)


