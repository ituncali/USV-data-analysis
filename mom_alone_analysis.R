maternal_counts_all <- count_frame_2 %>% filter(recording == "MomAlone" | recording == "PupsSep")

maternal.per <- maternal_counts_all %>% 
  group_by(strain, categories.allowed, recording) %>% 
  summarise(usv.count = sum(total.counts)) %>%
  group_by(strain) %>% mutate(percent=usv.count/sum(usv.count)*100) %>%
  filter(percent>2)

label.to.keep <- as.vector(unique(maternal.per$categories.allowed))

maternal_counts <- maternal_counts_all %>% filter(categories.allowed %in% label.to.keep)

maternal.per <- maternal_counts %>% 
  group_by(strain, rat.id, file.name, recording) %>% 
  mutate(usv.per = total.counts/sum(total.counts))
