library(tidyverse);library(qualtRics)
source("../resources.R")

# read simon data
df <-
  list.files(path='data/',pattern = "*.csv",full.names=T) %>% 
  map_df(~read_csv(., col_types = cols(trial_resp_timing.keys=col_character()))) %>% as.data.frame()

# remove practice blocks
df = df %>% filter(forcedRT_resp.rt > 0) %>% select(
  participant,
  block_id,
  forcedRT_resp.rt,
  forcedRT_resp.corr,
  target_onset,
  ontime,
  feedback_cond,
) 

# check subjects' avg performance
df.avg = df %>% group_by(participant) %>% summarise(medianRT_forced = median(forcedRT_resp.rt, na.rm=T),
                                                    meanACC_forced = mean(forcedRT_resp.corr, na.rm=T),
                                                    n = n())
a = df.avg$participant[abs(df.avg$medianRT_forced-2)>=.1]
b = df.avg$participant[df.avg$meanACC_forced<=.6]

bad_subs = unique(c(a, b))

# remove bad subjects
df = df %>% filter(!participant %in% bad_subs)

# identify blocks of trials where Ss were repeatedly slow or fast
repeat_offtime = df[,c('participant','block_id','ontime')] %>% 
  drop_na() %>% 
  group_by(participant, block_id) %>% 
  summarise(val1 =  rle(ontime)$values, run1 = rle(ontime)$lengths) %>%
  filter(val1!='Ontime' & run1 >= 10)

# identify blocks of trials where Ss repeatedly made errors
repeat_errors = df[,c('participant','block_id', 'forcedRT_resp.corr')] %>% 
  drop_na() %>% 
  group_by(participant, block_id) %>% 
  summarise(val2 =  rle(forcedRT_resp.corr)$values, run2 = rle(forcedRT_resp.corr)$lengths) %>%
  filter(val2!=1 & run2 >= 10)

# remove those blocks
df = plyr::join_all(list(df, repeat_offtime, repeat_errors), by=c('participant','block_id'),type='left') %>%
  filter(is.na(run1) & is.na(run2)) %>% select(-c(val1, run1, val2, run2))

# get a head count
df %>% group_by(feedback_cond) %>% summarise(n = length(unique(participant)))

# mark error trials
df = df %>% group_by(participant, block_id) %>% 
  mutate(previous_correct = dplyr::lag(forcedRT_resp.corr, n = 1, default = NA) %>% recode('0'='Incorrect','1'='Correct')) %>%
  drop_na()


# Forced-response plot
df %>%
  mutate(RT = forcedRT_resp.rt*1000, 
         PT = (forcedRT_resp.rt - target_onset)*1000, 
         accuracy = forcedRT_resp.corr) %>% 
  filter(!participant %in% bad_subs & block_id >= 1 & ontime=='Ontime')  %>% 
  group_by(previous_correct, feedback_cond) %>%
  summarise(t = 0:1000,
            mu = sw_smooth(x = PT, y = accuracy, FUN = mean, width = 100),
            er = sw_smooth(x = PT, y = accuracy, FUN = se, width = 100)) %>%
  ggplot(aes(x = t, y = mu, colour = previous_correct, fill = previous_correct)) +
  geom_ribbon(aes(ymin = mu - er, ymax = mu + er), colour = NA, alpha = .2) +
  geom_line(size = .8) + 
  facet_wrap(~feedback_cond)+
  scale_y_continuous("Response Accuracy", breaks = seq(0, 1, .25), 
                     labels = c("0", ".25", ".5", ".75", "1")) +
  scale_x_continuous("Preparation Time (ms)", breaks = seq(0, 1000, 100)) +
  scale_colour_manual(values = c(wong[6], wong[7])) +
  scale_fill_manual(values = c(wong[6], wong[7])) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1), expand = F) +

  theme_bw(base_size = 12, base_family = "sans") +
  theme(legend.position = c(.8,.4), legend.key.size = unit(.8, "cm"),
        axis.ticks.length = unit(0.1, "cm"), axis.ticks = element_line(size = .5), 
        axis.line = element_line(size = .5), panel.grid.minor = element_blank())

ggsave('memory_color.png', width = 5, height = 5)
