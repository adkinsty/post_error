library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(qualtRics)

setwd('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/color_reward_control')
library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

dat <- read_csv("dat/all_test_data.csv") %>% 
  mutate(exp = "color_reward_control", subj = id) %>%
  filter(PT > 50 & PT < 700) %>% filter(RT >= 1400 & RT <= 1600) %>%
  mutate(
    phase = ifelse(phase == 'location', .5, -.5),
    y = accuracy,
    err0 = ifelse(accuracyN0 == 0, .5, -.5), 
    rew = ifelse(rew == 5, .5, ifelse(rew == 1, -.5, 0)),
    time = PT / 1000) %>%
  filter(!is.na(accuracyN0))
  
write_csv(dat, "dat.csv")



dat %>% ggplot(aes(x = PT, group = phase, color = phase, fill = phase)) + geom_density()
dat %>% ggplot(aes(x = RT, group = phase, color = phase, fill = phase)) + geom_density()
