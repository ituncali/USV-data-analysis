mbt_counts_all <- count_frame_2 %>% filter(recording == "MBT")

mbt.per <- mbt_counts_all %>% 
       group_by(strain, categories.allowed) %>% 
       summarise(usv.count = sum(total.counts)) %>%
       group_by(strain) %>% mutate(percent=usv.count/sum(usv.count)*100) %>%
  filter(percent>2)

label.to.keep <- as.vector(unique(mbt.per$categories.allowed))

mbt_counts <- mbt_counts_all %>% filter(categories.allowed %in% label.to.keep)

library(nlme)
lme_mbt <- lme(total.counts ~ strain * categories.allowed, 
               random = ~1| rat.id/categories.allowed,
               data = mbt_counts)
anova.lme(lme_mbt)
#percentages
mbt.per <- mbt_counts %>% 
  group_by(strain, rat.id, file.name) %>% 
  mutate(usv.per = total.counts/sum(total.counts))
mbt.per.lme <- lme(usv.per~strain*categories.allowed, 
                   random=~1|rat.id/categories.allowed, 
                   data=mbt.per)
anova.lme(mbt.per.lme)
mbt.per.sum <- summary(lsmeans(mbt.per.lme, pairwise~strain|categories.allowed, adjust="Tukey"))[["contrasts"]]
View(mbt.per.sum[mbt.per.sum$p.value<0.05,])



##frequencies
mbt_freqs <- data_freqs %>% filter(recording == "MBT") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label %in% label.to.keep) %>%
  filter(!label=="short-c")
  #they are normally distributed within groups!
#need to take out short-c

lme_mbt_freq <- lme(mean.freq ~ strain * label, 
                    random = ~1|rat.id/label, 
                    data = mbt_freqs)
anova.lme(lme_mbt_freq)
mbt.lme.freq.sum <- summary(lsmeans(lme_mbt_freq, pairwise ~ strain |label, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.freq.sum[mbt.lme.freq.sum$p.value < 0.05,])


##durations
mbt_durs <- data_durs %>% filter(recording == "MBT") %>% group_by(strain, label,recording,file.name,rat.id) %>%
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label %in% label.to.keep) %>% filter(!label=="short-c")
#need to remove short-c again
lme_mbt_dur <- lme(mean.dur ~ strain * label,
                   random = ~1|rat.id/label,
                   data = mbt_durs)
anova.lme(lme_mbt_dur)
mbt.lme.dur.sum <- summary(lsmeans(lme_mbt_dur, pairwise~strain|label, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.dur.sum[mbt.lme.dur.sum$p.value < 0.05,])


