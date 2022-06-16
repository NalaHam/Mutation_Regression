
no_na_nonsense <- mut_count[!is.na(mut_count$Nonsense_Mutation), ]

par(mfrow=c(1,2))

plot(no_na_nonsense$Nonsense_Mutation~no_na_nonsense$age, xlab = "age in years", 
     ylab= "Nonsense Mutations")

plot(mut_count$Nonsense_Mutation~mut_count$age, xlab = "age in years", 
     ylab= "Nonsense Mutations")

sum(is.na(mut_count$Nonsense_Mutation))

sum(!is.na(mut_count$Nonsense_Mutation))

1256+9874

line_fit <- lm(Nonsense_Mutation~age, data = mut_count)
line_fit #slope is 0.08422, meaning that there is a 0.08422 increase in Nonsense mutations with every 1 year increase of age

line_fit_1 <- lm(Nonsense_Mutation~years_smoked, data = mut_count)
line_fit_1 # slope is 0.0898, meaning that there is a 0.0898 increase in Nonsense mutations with every 1 year increase of smoking
















    