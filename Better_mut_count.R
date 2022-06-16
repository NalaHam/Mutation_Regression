library(dplyr)
library(plyr)

n_distinct(pan_cancers_NA$sample_id)
unique(pan_cancers_NA$alcohol_intensity)

unique(pan_cancers_NA$Variant_Classification) #find the different types of mutations

#get info for each of the mutations. Gives TRUE or FALSE result
df <- summarise(pan_cancers_NA, Project, sample_id, Variant_Classification == "Silent", 
                Variant_Classification == "Nonsense_Mutation", 
                Variant_Classification == "Missense_Mutation", 
                Variant_Classification == "Splice_Region", 
                Variant_Classification == "Frame_Shift_Del",
                Variant_Classification == "Splice_Site",
                Variant_Classification == "Frame_Shift_Ins",
                Variant_Classification == "Nonstop_Mutation",
                Variant_Classification == "RNA",
                Variant_Classification == "In_Frame_Del",
                Variant_Classification == "Translation_Start_Site",
                Variant_Classification == "In_Frame_Ins",
                Variant_Classification == "NA", gender, cancer_level, age_at_index,
                race, ethnicity, years_smoked, cigarettes_per_day, alcohol_history) 

#Rename columns
names(df) <- c("Project", "sample_id", "Silent", "Nonsense_Mutation", "Missense_Mutation", "Splice_Region",
               "Frame_Shift_Del", "Splice_Site", "Frame_Shift_Ins", "Nonstop_Mutation",
               "RNA", "In_Frame_Del", "Translation_Start_Site", "In_Frame_Ins", 
               "unknown", "gender", "cancer_level", "age_at_index", 
               "race", "ethnicity", "years_smoked", "cigarettes_per_day", "alcohol_history") 

#change TRUE or FALSE to 1 or 0
df$Silent <- as.integer(as.logical(df$Silent))
df$Nonsense_Mutation <- as.integer(as.logical(df$Nonsense_Mutation))
df$Missense_Mutation <- as.integer(as.logical(df$Missense_Mutation))
df$Splice_Region <- as.integer(as.logical(df$Splice_Region))
df$Frame_Shift_Del <- as.integer(as.logical(df$Frame_Shift_Del))
df$Splice_Site <- as.integer(as.logical(df$Splice_Site))
df$Frame_Shift_Ins <- as.integer(as.logical(df$Frame_Shift_Ins))
df$Nonstop_Mutation <- as.integer(as.logical(df$Nonstop_Mutation))
df$RNA <- as.integer(as.logical(df$RNA))
df$In_Frame_Del <- as.integer(as.logical(df$In_Frame_Del))
df$Translation_Start_Site <- as.integer(as.logical(df$Translation_Start_Site))
df$In_Frame_Ins <- as.integer(as.logical(df$In_Frame_Ins))


mut_count <- ddply(df, "sample_id", numcolwise(sum)) #Make df of the sum of the columns with respect to the individual sample_ids. 

df_3 <- df[, c(1,2,16,17,18,19, 20, 21, 22, 23)] #Make df with just project, sample_id, and demographics

mut_count <- merge(mut_count, df_3, by = "sample_id") #combine the sum of mutations and the demographics by sample_id. This gives repeated rows of the same person.

mut_count <- distinct(mut_count, .keep_all = TRUE) #keep only the unique rows in mut_count

mut_count <- mut_count[, -c(14,15,16)] #remove extra rows

n_distinct(mut_count$sample_id)

sum(!is.na(mut_count$years_smoked.y))

names(mut_count)[, c(17, 20, 21)] <- c("age", "years_smoked", "cigarettes_per_day")

names(mut_count)[17] <- "age"
names(mut_count)[20] <- "years_smoked"
names(mut_count)[21] <- "cigarettes_per_day"

write.csv(mut_count, file = "full_mut_count.csv")






