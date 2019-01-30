library(dplyr)
library(tidyr)
library(ggplot2)
options(max.print=100)

d1 <- readRDS("out/res-2019-01-29-18-18-41.RDS")

df_res <- as.data.frame(d1$results)

nsim <- d1$cfg$nsims
# true values for median surv (mnth)
baseclin <- med_surv(d1$cfg$b0tte)
trtclin <- med_surv(d1$cfg$b0tte + d1$cfg$b1tte)
# true values for seroconversion prob
basesero <- d1$cfg$baselineprobsero
trtsero <- d1$cfg$trtprobsero
# accrual
accrual <- d1$cfg$people_per_interim_period
info_delay <- d1$cfg$sero_info_delay

outcome_labs <- c("Futile", "Superior", "Inconclusive")
df_res$trial_outcome <- NA
df_res$trial_outcome <- ifelse(df_res$stop_i_fut == 1 | df_res$stop_c_fut == 1, 1, df_res$trial_outcome)
df_res$trial_outcome <- ifelse(df_res$stop_c_sup == 1, 2, df_res$trial_outcome)
df_res$trial_outcome <- ifelse(df_res$inconclusive == 1, 3, df_res$trial_outcome)

df_res$trial_outcome <- factor(df_res$trial_outcome, levels = 1:3, labels = outcome_labs)

# df_res$trial_outcome <- ifelse(stop_v_samp == 1, 3, df_res$trial_outcome)

# trials are either futile or superior 
df_outcome <- df_res %>%
  dplyr::filter(trial_outcome %in% c("Futile", "Superior", "Inconclusive") )

df_stopv <- df_res %>%
  dplyr::filter(stop_v_samp == 1) 

n <- c(table(df_res$trial_outcome),   nrow(df_stopv))
prob <- c(prop.table(table(df_res$trial_outcome)),   nrow(df_stopv)/max(df_res$idxsim))
names(prob) <- c(names(prob)[1:3], "StopV")

df_ss <- df_outcome %>% 
  dplyr::group_by(trial_outcome) %>%
  dplyr::summarise(ss_mean = mean(n_obs), 
                   ss_sd = sd(n_obs)) %>%
  tidyr::complete(trial_outcome, fill = list(ss_mean = 0))



dtmp <- data.frame(idsim = d1$cfg$idsim,
                   nsim = nsim,
                   
                   basesero = basesero,
                   trtsero = trtsero,
                   baseclin = baseclin,
                   trtclin = trtclin,
                   accrual = accrual,
                   info_delay = info_delay,
                   
                   
                   prob_sup = prob["Superior"],
                   prob_fut = prob["Futile"],
                   prob_incon = prob["Inconclusive"],
                   prob_stopv = prob["StopV"],
                   
                   n_sup = n[2],
                   n_fut = n[1],
                   n_incon = n[3],
                   n_stopv = n[4],
                   
                   ss_fut_mean = unlist(df_ss[df_ss$trial_outcome == "Futile", "ss_mean"]), 
                   ss_fut_sd = unlist(df_ss[df_ss$trial_outcome == "Futile", "ss_sd"]), 
                   
                   ss_sup_mean = unlist(df_ss[df_ss$trial_outcome == "Superior", "ss_mean"]), 
                   ss_sup_sd = unlist(df_ss[df_ss$trial_outcome == "Superior", "ss_sd"]), 
                   
                   ss_incon_mean = unlist(df_ss[df_ss$trial_outcome == "Inconclusive", "ss_mean"]), 
                   ss_incon_sd = unlist(df_ss[df_ss$trial_outcome == "Inconclusive", "ss_sd"]), 
                   
                   ss_stopv_mean = mean(df_stopv$n_obs),
                   ss_stopv_sd = sd(df_stopv$n_obs))

dtmp <- dtmp %>%
  dplyr::mutate(f_prob_sup = sprintf("%.3f (%d)", prob_sup, n_sup),
                f_prob_fut = sprintf("%.3f (%d)", prob_fut, n_fut),
                f_prob_incon = sprintf("%.3f (%d)", prob_incon, n_incon),
                f_prob_stopv = sprintf("%.3f (%d)", prob_stopv, n_stopv),
                
                f_ss_sup = sprintf("%.0f (%.1f)", ss_sup_mean, ss_sup_sd),
                f_ss_fut = sprintf("%.0f (%.1f)", ss_fut_mean, ss_fut_sd),
                f_ss_incon = sprintf("%.0f (%.1f)", ss_incon_mean, ss_incon_sd),
                f_ss_stopv = sprintf("%.0f (%.1f)", ss_stopv_mean, ss_stopv_sd) )

df_0 <- dtmp %>%
  dplyr::arrange(-accrual, baseclin, basesero) %>%
  dplyr::mutate(row = 1:n()) %>%
  dplyr::select(row, nsim, 
                basesero, trtsero, baseclin, trtclin, accrual,info_delay,
                f_prob_sup, f_prob_fut, f_prob_incon, f_prob_stopv, 
                f_ss_sup, f_ss_fut, f_ss_incon, f_ss_stopv,
                idsim) %>%
  dplyr::select(-row)


digits <- c(2, 2, 0, 0, 0, 2,
            0, 0, 0, 0,
            0, 0, 0, 0)
options(knitr.kable.NA = '-')
#stopifnot(nrow(df_0) == 24)
# Null01 SeroTrt01 ClinTrt01 ClinTrt02 ClinSeroTrt01 ClinSeroTrt02 ClinSeroTrt03 ClinSeroTrt04
kable(df_0 %>% dplyr::filter(idsim == "Test_00") %>% dplyr::select(-idsim, -nsim) %>%
        dplyr::arrange(basesero, baseclin, -accrual, info_delay),
      caption = "Table 1. Simulation Null Scenarios",
      col.names = c("base", "trt", "base", "trt", "accrual", "info delay",
                    "superiority", "futility", "inconclusive", "stop v samp",
                    "superiority", "futility", "inconclusive", "stop v samp"),
      digits = digits) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = F, position = "left",
                font_size = 11,
                latex_options = "hold_position") %>%
  add_header_above(c("Seroconversion" = 2, "Clinical" = 2,
                     " " = 1, " " = 1,
                     "Probability (number of trials)" = 4,
                     "Mean sample size (sd)" = 4))


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
