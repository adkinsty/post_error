library(tidyverse)

setwd("/Users/adkinsty/Dropbox (University of Michigan)/LeeLab/Experiments/Exp_files/forced_response/post_error")

memory_id <- read_csv("memory/memory_participants.csv")
memory2_id <- read_csv("memory2/memory2_participants.csv")
memory_color_short_id <- read_csv("memory_color_short/memory_color_short_participants.csv")
memory_color_long_id <- read_csv("memory_color_long/memory_color_long_participants.csv")

participant_ids <- bind_rows(memory_id, memory2_id, memory_color_short_id, memory_color_long_id)

write_csv(participant_ids, "participant_ids.csv")
