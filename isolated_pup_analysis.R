pup_counts_all <- count_frame_2 %>% filter(recording == "Mpupiso"|
                                             recording=="Fpupiso")
pup.per <- pup_counts_all %>% 
  group_by(strain, categories.allowed, recording) %>% 
  summarise(usv.count = sum(total.counts)) %>%
  group_by(strain) %>% mutate(percent=usv.count/sum(usv.count)*100) %>%
  filter(percent>2)

label.to.keep.pup <- as.vector(unique(pup.per$categories.allowed))

pup_counts <- pup_counts_all %>% filter(categories.allowed %in% label.to.keep.pup)
pup_counts$rat.id.fix <- as.factor(ifelse(pup_counts$recording=="Fpupiso",paste(pup_counts$rat.id,"F"),paste(pup_counts$rat.id,"M")))


lme_pups <- lme(total.counts ~ strain * recording * categories.allowed, 
                random = ~1| rat.id.fix/categories.allowed,
                data = pup_counts)
anova.lme(lme_pups)

#percentages
pup.per <- pup_counts %>% 
  group_by(strain, recording, rat.id, rat.id.fix) %>% 
  mutate(usv.per = total.counts/sum(total.counts))
pup.per.lme <- lme(usv.per~strain*categories.allowed, 
                   random=~1|rat.id.fix/categories.allowed, 
                   data=pup.per)
anova.lme(mbt.per.lme)


##Frequencies
pup_freqs <- data_freqs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label %in% label.to.keep.pup)
pup_freqs$rat.id.fix <- ifelse(pup_freqs$recording=="Fpupiso",paste(pup_freqs$rat.id,"F"),paste(pup_freqs$rat.id,"M"))
lme_pup_freq <- lme(mean.freq ~ strain * recording * label, 
                    random = ~1|rat.id.fix,
                    data = pup_freqs)
anova(lme_pup_freq)

##Durations
pup_durs <- data_durs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  group_by(strain, label,recording,file.name,rat.id) %>%
  summarise(mean.dur = mean(duration), SEM = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label %in% label.to.keep.pup)
pup_durs$rat.id.fix <- ifelse(pup_durs$recording=="Fpupiso",paste(pup_durs$rat.id,"F"),paste(pup_durs$rat.id,"M"))




