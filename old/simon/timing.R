library(tidyverse)

dat <- read_csv("all_dat.csv") %>%
  mutate(
    too_slow = ifelse(RT > 1600, 1, 0), 
    too_fast = ifelse(RT < 1400, 1, 0),
    outcome = ifelse(too_slow, 'too_slow', ifelse(too_fast, 'too_fast', 'well_timed'))) %>%
  group_by(id) %>% 
  mutate(outcome0 = lag(outcome)) %>%
  filter(!is.na(outcome0))

dat %>% group_by(congruencyN0, congruencyN1) %>% 
  summarise(
    accuracy = mean(accuracy),
    too_slow = mean(outcome == 'too_slow'),
    too_fast = mean(outcome == 'too_fast'),
    well_timed = mean(outcome == 'well_timed'))

dat %>% group_by(congr, err0) %>% 
  summarise(too_slow = mean(too_slow), 
    too_fast = mean(too_fast))

dat %>% filter(PT >=0 & PT <= 1000) %>%
  ggplot(aes(x = PT, y = too_slow, colour = factor(congr), fill = factor(congr))) +
  geom_smooth() + 
  facet_wrap(err0 ~ .) +
  scale_alpha_manual(values = c(.25, .5))

dat %>% filter(PT >=0 & PT <= 1000) %>%
  ggplot(aes(x = PT, y = too_fast, colour = factor(congr), fill = factor(congr))) +
  geom_smooth() + 
  facet_wrap(err0 ~ .) +
  scale_alpha_manual(values = c(.25, .5))
