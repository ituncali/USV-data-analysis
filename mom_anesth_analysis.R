momanesth_counts_all <- count_frame_2 %>% filter(recording == "MA")

ma.per <- momanesth_counts_all %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(usv.count = sum(total.counts)) %>%
  group_by(strain) %>% mutate(percent=usv.count/sum(usv.count)*100) %>%
  filter(percent>2)

label.to.keep <- as.vector(unique(ma.per$categories.allowed))

momanesth_counts <- momanesth_counts_all %>% filter(categories.allowed %in% label.to.keep)


ma.per <- momanesth_counts %>% 
  group_by(strain, rat.id, file.name) %>%
  mutate(usv.per=total.counts/sum(total.counts)) 


##frequencies
momanesth_freqs <- data_freqs %>% filter(recording == "MA") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label %in% label.to.keep)


##durations
momanesth_durs <- data_durs %>% filter(recording == "MA") %>% 
  group_by(strain, label, file.name,recording, rat.id) %>% 
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label %in% label.to.keep)

