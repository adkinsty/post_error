library(tidyverse);library(qualtRics)
source("../resources.R")

# read simon data
df <-
  list.files(path='data/',pattern = "*.csv",full.names=T) %>% 
  map_df(~read_csv(., col_types = cols(trial_resp_timing.keys=col_character()))) %>% as.data.frame()

# check subjects' avg performance
df.avg = df %>% group_by(participant) %>% summarise(medianRT_forced = median(forcedRT_resp.rt, na.rm=T),
                                                    meanACC_forced = mean(forcedRT_resp.corr, na.rm=T),
                                                    n = n())

# check if they were repeatedly slow or fast
t = df[,c('participant','ontime')] %>% 
  drop_na() %>% 
  group_by(participant) %>% 
  summarise(values =  rle(ontime)$values, run = rle(ontime)$lengths) %>%
  filter(values!='Ontime' & run >= 10)

# check % of on-time trials
k = df %>% filter(forcedRT_resp.rt > 0) %>% group_by(participant) %>% summarise(on_time = mean(ontime=='Ontime'))

a = df.avg$participant[abs(df.avg$medianRT_forced-2)>=.1]
b = df.avg$participant[df.avg$meanACC_forced<=.6]
c = t$participant %>% unique()
d = k$participant[k$on_time < .5]

bad_subs = unique(c(a, b, c))

# mark error trials
df = df %>% filter(forcedRT_resp.rt > 0) %>% select(
  participant,
  block_id,
  forcedRT_resp.rt,
  forcedRT_resp.corr,
  target_onset,
  ontime,
) %>% group_by(participant, block_id) %>% 
  mutate(previous_correct = dplyr::lag(forcedRT_resp.corr, n = 1, default = NA) %>% recode('0'='Incorrect','1'='Correct')) %>%
  drop_na()


# Forced-response plot
df %>%
  mutate(RT = forcedRT_resp.rt*1000, 
         PT = (forcedRT_resp.rt - target_onset)*1000, 
         accuracy = forcedRT_resp.corr) %>% 
  filter(!participant %in% bad_subs & block_id >= 1 & ontime=='Ontime')  %>% 
  group_by(previous_correct) %>%
  summarise(t = 0:1000,
            mu = sw_smooth(x = PT, y = accuracy, FUN = mean, width = 100),
            er = sw_smooth(x = PT, y = accuracy, FUN = se, width = 100)) %>%
  ggplot(aes(x = t, y = mu, colour = previous_correct, fill = previous_correct)) +
  geom_ribbon(aes(ymin = mu - er, ymax = mu + er), colour = NA, alpha = .2) +
  geom_line(size = .8) + 
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

ggsave('simon_arrow.png', width = 5, height = 5)
