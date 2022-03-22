library(tidyverse) # Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686
#library(qualtRics)

exp_name <- 'memory'

setwd(sprintf('/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error/%s/', exp_name))
source("../../../forced_response/func.R")

get_dat <- function() {
  # adapted from Han's analysis script

  df <- list.files(path='dat/',pattern = '*.csv',full.names=T) %>% 
    map_df(~read_csv(., col_types = cols(trial_resp_timing.keys=col_character()))) %>% as.data.frame()

  tim <- df[,c('participant','ontime','block_id')] %>% 
    drop_na() %>% 
    group_by(participant, block_id) %>% 
    summarise(values = rle(ontime)$values, run = rle(ontime)$lengths) %>%
    filter(values!='Ontime' & run >= 10)

  kerr <- df[,c('participant','forcedRT_resp.corr', 'forcedRT_resp.rt', 'target_onset','block_id')] %>% 
    drop_na() %>% 
    mutate(pt = forcedRT_resp.rt - target_onset,
      bad_error = !forcedRT_resp.corr & pt > .5) %>%
    group_by(participant, block_id) %>% 
    summarise(values = rle(bad_error)$values, run = rle(bad_error)$lengths) %>%
    filter(values & run >= 5)

  id_tim <- tim$participant %>% unique()
  id_err <- kerr$participant %>% unique()
  bad_subs = unique(c(id_tim, id_err))

  exp_dat <- df %>%
    filter(block_id >= 1) %>%
    #filter(!(participant %in% bad_subs) & block_id >= 1) %>%
    filter(!is.na(forcedRT_resp.rt)) %>%
    mutate(RT = forcedRT_resp.rt*1000, 
          PT = (forcedRT_resp.rt - target_onset)*1000 / 2, #divide by 2 to equate with other experiments...
          PT_scale = 2,
          accuracy = forcedRT_resp.corr,
          key = forcedRT_resp.keys, 
          stim = correct_key) %>% 
    group_by(participant) %>%
    mutate(
      acc00 = lag(accuracy, n = 2),
      acc0 = lag(accuracy),
      time0 = lag(PT), 
      key0 = lag(key),
      stim0 = lag(stim)) %>%
    mutate(
      y = accuracy,
      exp = 'memory', 
      err00 = ifelse(acc00 == 0, .5, -.5),
      err0 = ifelse(acc0 == 0, .5, -.5), 
      time = PT / 1000,
      key_rep = key0 == key,
      stim_rep = stim0 == stim,
      key_stim_rep = key0 == stim, 
      time_rep = time0 <= time + 100 & time >= time - 100) %>% tibble() %>%
    filter(!is.na(err0)) %>%
    filter(RT >= 1900 & RT <= 2100) %>%
    filter(PT > 0 & PT < 1000)

    
  return(exp_dat)
}

dat <- get_dat()

write_csv(dat, 'dat.csv')


dat %>% group_by(participant) %>%
    summarise(exp = exp[1]) %>%
    write_csv(sprintf("%s_participants.csv", exp_name))

# df %>% ggplot(aes(x = forcedRT_resp.rt)) +
#   geom_histogram() +
#   scale_x_continuous(breaks = seq(1, 3, .1), limits = c(1, 3)) +
#   geom_vline(xintercept = 1.9, colour = 'green') +
#   geom_vline(xintercept = 2.1, colour = 'green')

# dat %>% ggplot(aes(x = PT/1000)) +
#   geom_histogram()# +
#   #scale_x_continuous(breaks = seq(-1, 2, .2), limits = c(-1, 2))


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
