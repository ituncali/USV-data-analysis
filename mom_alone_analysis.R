maternal_counts_all <- count_frame_2 %>% filter(recording == "MomAlone" | recording == "PupsSep")

maternal.per <- maternal_counts_all %>% 
  group_by(strain, categories.allowed, recording) %>% 
  summarise(usv.count = sum(total.counts)) %>%
  group_by(strain) %>% mutate(percent=usv.count/sum(usv.count)*100) %>%
  filter(percent>2)

label.to.keep.maternal <- as.vector(unique(maternal.per$categories.allowed))

maternal_counts <- maternal_counts_all %>% filter(categories.allowed %in% label.to.keep.maternal)

maternal.per <- maternal_counts %>% 
  group_by(strain, rat.id, file.name, recording) %>% 
  mutate(usv.per = total.counts/sum(total.counts))

##frequencies
maternal_freqs <- data_freqs %>% filter(recording == "MomAlone" | recording == "PupsSep") %>% 
  group_by(strain, label, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label %in% label.to.keep.maternal) %>%
  filter(!label=="unclear")


maternal_durs <- data_durs %>% filter(recording == "MomAlone" | recording == "PupsSep") %>% 
  group_by(strain, label, rat.id) %>% 
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label %in% label.to.keep.maternal)


