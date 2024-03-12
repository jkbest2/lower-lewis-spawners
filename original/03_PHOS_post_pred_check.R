library(tidyverse)
library(magrittr)
library(rjags)
library(lubridate)

#Process Data
source("original/00_PHOS_DataPreparation.R")

#Load Model Code and Functions
source("original/01_PHOS_functions_and_JAGSmodel.R")

#Load saved model object
m.PHOS = readRDS("data/pHOS_mcmc_differentRE_04292018.rds")

# Assign mcmc draws to data frame
draws <- as_data_frame(as.matrix(m.PHOS$mcmc))

# Get sample size from processed data
N = nrow(obs_JAGS)

# randomly select 1000 posterios draws
samps <- sample(1:nrow(draws), 1000, replace = F)

# run posterior simulation for single posterior draw
pp_pHOS = post_pred_pHOS(N, draws[samps[1],]) %>% 
  mutate(iter = 1)

# plus 999 more
for( i in 2:1000)
  pp_pHOS = bind_rows(pp_pHOS, post_pred_pHOS(N, draws[samps[i],])  %>% 
                        mutate(iter = i))

# Summarize posterior simmulations and compare to observed data

Post_Pred <- left_join(left_join(pp_pHOS %>% 
                                   group_by(event, iter) %>% 
                                   summarise(N_tangle = sum(n_return + n_spawn),
                                             N_trap   = sum(n_trap)) %>% 
                                   ungroup() %>% 
                                   group_by(event) %>% 
                                   summarise(avg_tangle = mean(N_tangle),
                                             tangle_10  = quantile(N_tangle, .1),
                                             tangle_90  = quantile(N_tangle, .9),
                                             avg_trap   = mean(N_trap),
                                             trap_10    = quantile(N_trap, .1),
                                             trap_90    = quantile(N_trap, .9)) %>% 
                                   mutate(event = recode(event, 
                                                         "Event_00" = "Event 00",
                                                         "Event_01" = "Event 01",
                                                         "Event_02" = "Event 02",
                                                         "Event_03" = "Event 03",
                                                         "Event_04" = "Event 04",
                                                         "Event_05" = "Event 05",
                                                         "Event_06" = "Event 06",
                                                         "Event_07" = "Event 07",
                                                         "Event_08" = "Event 08",
                                                         "Event_09" = "Event 09",
                                                         "Event_10" = "Event 10",
                                                         "Event_11" = "Event 11"
                                   )),
                                 
                                 trap_observed %>% 
                                   gather(-ORIGIN, -GENDER, -LIFE, -PIT_ID,
                                          key = event, value = obs_state) %>% 
                                   group_by(event) %>%
                                   summarise(n_trap_obs = sum(obs_state == 4)) %>% 
                                   ungroup() %>% 
                                   mutate(event = recode(event, 
                                                         "Event_0" = "Event 00",
                                                         "Event_1" = "Event 01",
                                                         "Event_2" = "Event 02",
                                                         "Event_3" = "Event 03",
                                                         "Event_4" = "Event 04",
                                                         "Event_5" = "Event 05",
                                                         "Event_6" = "Event 06",
                                                         "Event_7" = "Event 07",
                                                         "Event_8" = "Event 08",
                                                         "Event_9" = "Event 09",
                                                         "Event_10" = "Event 10",
                                                         "Event_11" = "Event 11"
                                   ))),
                       
                       observed_state %>% 
                         gather(-ORIGIN, -GENDER, -LIFE, -PIT_ID,
                                key = event, value = obs_state) %>% 
                         group_by(event) %>%
                         summarise(n_tangle_obs = sum(obs_state %in% c(2, 3))) %>% 
                         ungroup() %>% 
                         mutate(event = recode(event, 
                                               "Event_0" = "Event 00",
                                               "Event_1" = "Event 01",
                                               "Event_2" = "Event 02",
                                               "Event_3" = "Event 03",
                                               "Event_4" = "Event 04",
                                               "Event_5" = "Event 05",
                                               "Event_6" = "Event 06",
                                               "Event_7" = "Event 07",
                                               "Event_8" = "Event 08",
                                               "Event_9" = "Event 09",
                                               "Event_10" = "Event 10",
                                               "Event_11" = "Event 11")))


g_1 <- ggplot(Post_Pred,
              aes(x = event, y = n_tangle_obs)) +
  geom_bar(stat = "identity", fill = "slategrey") +
  geom_pointrange(aes(y = avg_tangle, 
                      ymin = tangle_10, 
                      ymax = tangle_90)) + 
  ylab("Number Captured per Tangle Net Event \n Observed and Replicated") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())


g_2 <- ggplot(Post_Pred,
              aes(x = event, y = n_trap_obs)) +
  geom_bar(stat = "identity", fill = "slategrey") +
  geom_pointrange(aes(y = avg_trap, 
                      ymin = trap_10, 
                      ymax = trap_90)) + 
  ylab("Merwin Trap Number Captured \n Between each Tangle Net Event \n Observed and Replicated") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank())


cowplot::plot_grid(g_1, g_2, align = "v", ncol = 1)
