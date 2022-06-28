
regmult_8 <- lm(total_mut~factor(gender)+age+factor(alcohol_history),
                data = eliminated)

summary(regmult_8)

nobs(regmult_8)#6111
#female, male, and age stay significant. R^2 is 0.106

res <- resid(regmult_8)

#see plot of residuals of regression, for it to be good data there needs to be random placement of points around y=o
plot(fitted(regmult_8), res)
abline(0,0)

#Q-Q normal plots show if the data is normal or not. If it follows the line it is normal
qqnorm(res)
qqline(res)
#it doesn't follow the line on the ends ==> it is not normal

#do a distribution plot to also check for normacy, if it is bell shaped then it is normal
plot(density(res))

#can also do a histogram to see if normally distributed
hist(eliminated$total_mut)

#plot of residuals and normal Q-Q
par(mfrow = c(2,2))
plot(regmult_8)



