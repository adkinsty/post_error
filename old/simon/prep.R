library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(qualtRics)

setwd('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/simon/')


dat <- read_csv('dat/all_test_data.csv') %>% 
  mutate(exp = 'simon', rew = 0, RT = trial_resp.rt * 1000, congr = congruencyN1, congr0 = congruencyN0, subj = id) %>%
  filter(PT > 50 & PT < 700) %>% filter(RT >= 1400 & RT <= 1600) %>%
  mutate(
    phase = ifelse(exp != 'color_reward_control', 'none', phase),
    y = accuracy,
    err0 = ifelse(accuracyN0 == 0, .5, -.5), 
    congr = ifelse(congr == 'congruent', .5, -.5), 
    congr0 = ifelse(congr0 == 'incongruent', .5, -.5), 
    rew = ifelse(rew == 5, .5, ifelse(rew == 1, -.5, 0)),
    time = PT / 1000) %>%
  filter(!is.na(congr0))
  
write_csv(dat, 'dat.csv')



