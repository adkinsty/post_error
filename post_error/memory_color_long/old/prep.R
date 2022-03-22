library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(qualtRics)

exp_name <- 'memory_color_long'

setwd(sprintf('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/%s/', exp_name))
source("func.R")

get_dat <- function() {
  # adapted from Han's analysis script

  df <- list.files(path='dat/',pattern = '*.csv',full.names=T) %>% 
    map_df(~read_csv(., col_types = cols(trial_resp_timing.keys=col_character()))) %>% as.data.frame() %>% 
    filter(feedback_cond == 'long')

  tim <- df[,c('participant','ontime')] %>% 
    drop_na() %>% 
    group_by(participant) %>% 
    summarise(values =  rle(ontime)$values, run = rle(ontime)$lengths) %>%
    filter(values!='Ontime' & run >= 10)

  id <- tim$participant %>% unique()
  bad_subs = unique(c(id))

  exp_dat <- df %>%
    mutate(RT = forcedRT_resp.rt*1000, 
          PT = (forcedRT_resp.rt - target_onset)*1000, 
          accuracy = forcedRT_resp.corr,
          key = forcedRT_resp.keys, 
          stim = correct_key) %>% 
    group_by(participant) %>%
    mutate(acc0 = lag(accuracy),
      time0 = lag(PT), 
      key0 = lag(key),
      stim0 = lag(stim)) %>%
    filter(!participant %in% bad_subs & block_id >= 1 & ontime=='Ontime') %>%
    filter(PT > 0 & PT < 1000) %>%
    filter(!is.na(forcedRT_resp.rt)) %>%
    mutate(
      y = accuracy,
      iti = feedback_cond,
      exp = 'memory_color_long',
      err0 = ifelse(acc0 == 0, .5, -.5), 
      time = PT / 1000,
      key_rep = key0 == key,
      stim_rep = stim0 == stim,
      key_stim_rep = key0 == stim, 
      time_rep = time0 <= time + 100 & time >= time - 100) %>% tibble() %>%
    filter(!is.na(err0))
    
  return(exp_dat)
}

dat <- get_dat()

write_csv(dat, 'dat.csv')


# make directories
if (!dir.exists('fit')) {
  dir.create('fit')
}
if (!dir.exists('fig')) {
  dir.create('fig')
}
if (!dir.exists('loo')) {
  dir.create('loo')
}
if (!dir.exists('tab')) {
  dir.create('tab')
}
