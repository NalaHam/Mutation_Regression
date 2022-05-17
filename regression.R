library(dplyr)
library(ggplot2)
library(plyr)
library(broom)
library(ggpubr)

unique(BLCA_Mut$sample_id)
unique(BLCA_Mut$Variant_Classification) #see how many types of mutations there are

nonsense <- BLCA_Mut %>% filter(Variant_Classification == "Nonsense_Mutation") #make data with just nonsense mutations

Nonsense <- subset(nonsense, select = c("X", "sample_id", "Variant_Classification",
                                        "cancer_level", "gender", "age_at_index")) #subset data columns

df <- count(Nonsense, 'sample_id') #get number of nonsense mutations 

df <- merge(df, Nonsense, by = "sample_id", all = T) #add the count of nonsense mutations to each of the individuals in df

df_A <- df[, -3] #get rid of X column 

df_A <- df_A %>% 
  group_by(sample_id) %>%
  summarise_all(~list(unique(.[!is.na(.)]))) %>% 
  print.data.frame #compress data frame to just be a list of the unique people 

unique(df_A$cancer_level) #see what kinds of different stages there are

data <- filter(df_A, df_A$cancer_level == "Stage I" | df_A$cancer_level == "Stage II") #filter data so that it only has stage I and stage II people

data$sex.f <- ifelse(data$gender == 'male', 1, 0) #make male = 1 and female = 0 

data$X <- seq.int(nrow(data)) #make a column for the number or each row (not necessary, I just did it)

data$age_at_index <- as.numeric(data$age_at_index)#change the age from a character to a value
data$freq <- as.numeric(data$freq) #change number of mutations from a character to a value

cor(data$sex.f, data$age_at_index) #see the correlated between the two varibles. A small cor is needed for a regresssion

hist(data$freq)# done to see if data is normal or not. It is not. a normal is needed for a simple regression

plot(freq ~ sex.f, data = data)
plot(freq ~ age_at_index, data = data)

mutation.lm <- lm(freq ~ sex.f + age_at_index, data = data)
summary(mutation.lm)

par(mfrow = c(2,2))
plot(mutation.lm)
