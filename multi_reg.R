

mut_count$total_mut <- rowSums(mut_count[,2:13], na.rm = TRUE)

mut_count <- mut_count[, c(1,2,3,4,5,6,7,8,9,10,11,12,13,23,14,15,16,17,18,19,20,21,22)]

is.factor(mut_count$race)
is.factor(mut_count$cancer_level)
is.factor(mut_count$Project)
is.factor(mut_count$gender)
is.factor(mut_count$ethnicity)

regmult <- lm(total_mut~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(regmult)

nonse_reg <- lm(Nonsense_Mutation~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(nonse_reg)

miss_reg <- lm(Missense_Mutation~factor(gender)+age+factor(race)+factor(ethnicity)
                +years_smoked+factor(alcohol_history),
                data = mut_count)
summary(miss_reg)

names(mut_count)

sil_reg <- lm(Silent~factor(gender)+age+factor(race)+factor(ethnicity)
                +years_smoked+factor(alcohol_history),
                data = mut_count)
summary(sil_reg)

spl_reg <- lm(Splice_Region~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(spl_reg)

fram_d_reg <- lm(Frame_Shift_Del~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(fram_d_reg)

spl_s_reg <- lm(Splice_Site~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(spl_s_reg)

fram_i_reg <- lm(Frame_Shift_Ins~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(fram_i_reg)

nonst_reg <- lm(Nonstop_Mutation~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(nonst_reg)

RNA_reg <- lm(RNA~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(RNA_reg)

in_fd_reg <- lm(In_Frame_Del~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(in_fd_reg)

trans_reg <- lm(Translation_Start_Site~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(trans_reg)

in_fi_reg <- lm(In_Frame_Ins~factor(gender)+age+factor(race)+factor(ethnicity)
              +years_smoked+factor(alcohol_history),
              data = mut_count)
summary(in_fi_reg)

new_reg <- lm(total_mut~factor(gender)+age+factor(race)+factor(ethnicity)
                +years_smoked+factor(alcohol_history)+factor(Project),
                data = mut_count)
summary(new_reg)


