library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
library(qualtRics)

setwd('/Users/adkinsty/Dropbox (University of Michigan)/from_box/side_projects/distraction/post_error/stroop/')
library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

get_dat <- function() {
  # adapted from Han's analysis script

  df <- list.files(path='dat/',pattern = '*.csv',full.names=T) %>% 
    map_df(~read_csv(.)) %>% as.data.frame()

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
          text_key = str_sub(text, 1, 1),
          color_key = str_sub(letterColor, 1, 1)) %>% 
    group_by(participant) %>%
    mutate(acc0 = lag(accuracy)) %>%
    filter(!participant %in% bad_subs & block_id >= 1 & ontime=='Ontime') %>%
    filter(PT > 50 & PT < 700) %>%
    filter(!is.na(forcedRT_resp.rt)) %>%
    #filter(!finger_err) %>%
    mutate(
      y = accuracy,
      exp = 'stroop', 
      err0 = ifelse(acc0 == 0, .5, -.5), 
      congr = ifelse(trial_type == 'congruent', .5, -.5), 
      time = PT / 1000) %>% tibble() %>%
    filter(!is.na(err0))
    
  return(exp_dat)
}

dat <- get_dat()

write_csv(dat, 'dat.csv')

