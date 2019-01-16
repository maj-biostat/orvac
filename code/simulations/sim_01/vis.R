library(dplyr)
library(tidyr)
library(ggplot2)
options(max.print=100)

df <- as.data.frame(readRDS("results-2019-01-06-11:09.RDS"))
colnames(df)
head(df,50)



# how often do we stop?

# select the first instance of stop_trial_immu = 1 by simidx

df_sero_stop <- df %>%
  dplyr::filter(look == 1:5) %>%
  dplyr::group_by(simidx, trial_stop = stop_trial_immu==1) %>% 
  dplyr::arrange(simidx) %>%
  dplyr::mutate(count= row_number(), 
         first_in_sim = stop_trial_immu & count==1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-count, -trial_stop) %>%
  dplyr::filter(first_in_sim == T) %>%
  dplyr::select(simidx, look, n_obs, n_max, n_max_sero, immu_ppos_n, immu_ppos_max, stop_trial_immu) 

df_sero_stop %>% print(n = 25)

nrow(df_sero_stop)/max(df$simidx)



# tte

df_tte_stop <- df %>%
  dplyr::group_by(simidx, trial_stop = stop_trial_clin==1) %>% 
  dplyr::arrange(simidx) %>%
  dplyr::mutate(count= row_number(), 
                first_in_sim = stop_trial_clin & count==1) %>%
  dplyr::ungroup() %>%
  dplyr::select(-count, -trial_stop) %>%
  dplyr::filter(first_in_sim == T) %>%
  dplyr::select(simidx, look, n_obs, n_max, n_max_sero, clin_ppos_n, clin_ppos_max, stop_trial_clin) 

df_tte_stop %>% print(n = 25)

nrow(df_tte_stop)/max(df$simidx)
