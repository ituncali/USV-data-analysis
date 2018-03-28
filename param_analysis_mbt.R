data_counts <- read.csv(file.choose(),stringsAsFactors = F)
source("C:/Users/ituncali/Documents/Master's Thesis/USV-Data-Management/src/count_total.R")
library(stringr)
library(dplyr)

mbt_counts <- count_frame_3 %>% filter(recording == "MBT")

##counts analysis
library(nlme)
lme_mbt <- lme(total.counts ~ strain * categories.allowed, 
                random = ~1| rat.id,
                data = mbt_counts[mbt_counts$categories.allowed!="flat",])
anova.lme(lme_mbt)

library(lsmeans)
mbt.lme.sum <- summary(lsmeans(lme_mbt, pairwise ~ categories.allowed|strain, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.sum[mbt.lme.sum$p.value < 0.05,])

#frequency data

mbt_freqs <- data_freqs %>% filter(recording == "MBT") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label == "flat" | label == "flat-z" | label == "short")

#shapiro wilks test of normality
mbt_freqs_sw <- data_freqs %>% filter(recording == "MBT") %>% 
  select(strain, duration, m.freq, label)
stargazer(mbt_freqs_sw, type = "text")

shapiro.test(log(pup_freqs_sh[mbt_freqs_sw$strain == "SD" & mbt_freqs_sw$label == "flat",]$m.freq))


##frequency analysis
lme_mbt_freq <- lme(mean.freq ~ strain * label, 
                    random = ~1|rat.id, 
                    data = mbt_freqs)
anova.lme(lme_mbt_freq)
mbt.lme.freq.sum <- summary(lsmeans(lme_mbt_freq, pairwise ~ strain |label, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.freq.sum[mbt.lme.freq.sum$p.value < 0.05,])


#duration data test normality
mbt_durs_sw <- data_durs %>% filter(recording == "MBT")

#duration data
mbt_durs <- data_durs %>% filter(recording == "MBT") %>% group_by(strain, label,recording,file.name,rat.id) %>%
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label == "flat" | label == "flat-z" | label == "short")


##duration analysis
lme_mbt_dur <- lme(mean.dur ~ strain * label,
                   random = ~1|rat.id,
                   data = mbt_durs)
anova.lme(lme_mbt_dur)
mbt.lme.dur.sum <- summary(lsmeans(lme_mbt_dur, pairwise~strain|label, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.dur.sum[mbt.lme.dur.sum$p.value < 0.05,])


##maternal behavior data and analysis
mbt_behavior_unmelted <- read.csv("data/Exp1_MB_Scores.csv",stringsAsFactors = F)
library(reshape2)
mbt_behavior <- melt(data=mbt_behavior_unmelted,id = c("strain","rat.id", "file.name"),
                     variable.name="behavior",value.name = "score")

#behavioral counts
mbt.beh.counts.data <- mbt_behavior %>% filter(behavior == "retrieval" |
                                                 behavior == "mouthing" |
                                                 behavior == "corporal" |
                                                 behavior == "anogenital" |
                                                 behavior == "nest.building")
lme_mbt_beh_counts <- lme(score ~ strain * behavior,
                          random = ~1|rat.id,
                          data = mbt.beh.counts.data)
anova.lme(lme_mbt_beh_counts)
mbt.lme.beh.counts.sum <- summary(lsmeans(lme_mbt_beh_counts, pairwise~strain * behavior, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.beh.counts.sum[mbt.lme.beh.counts.sum$p.value < 0.05,])


#behavioral latencies
mbt.beh.lats.data <- mbt_behavior %>% filter(behavior == "lat.retrieve" |
                                               behavior == "lat.group" |
                                               behavior == "lat.hover" |
                                               behavior == "lat.nurse")
lme_mbt_beh_lats <- lme(score ~ strain * behavior,
                        random = ~1|rat.id,
                        data = mbt.beh.lats.data)
anova.lme(lme_mbt_beh_lats)
mbt.lme.beh.lats.sum <- summary(lsmeans(lme_mbt_beh_lats, pairwise~strain |behavior, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.beh.lats.sum[mbt.lme.beh.lats.sum$p.value < 0.05,])

#behavioral durations
mbt.beh.durs.data <- mbt_behavior %>% filter(behavior == "dur.hover" |
                                               behavior == "dur.nurse")

lme_mbt_beh_durs <- lme(score ~ strain * behavior,
                        random = ~1|rat.id,
                        data = mbt.beh.durs.data)
anova.lme(lme_mbt_beh_durs)


#usvs during behaviors
beh_times_data <- read.csv("data/Exp1_MBT_beh_times.csv",stringsAsFactors = F)
names(beh_times_data) <- c("rat.id","file.name","strain","behavior","start","end","beh.dur")
beh_times_data$behavior <- gsub("HKP","LKP",beh_times_data$behavior)
breaks <- c(rbind(beh_times_data$start,beh_times_data$end))

start.mbt <- start.rows %>% filter(recording == "MBT")
merged.data <- left_join(start.mbt, beh_times_data, by = c("file.name","rat.id","strain"))

beh.usv.counts.data <- merged.data %>% 
  ungroup() %>%
  filter(start.time > start & start.time < end) %>%
  group_by(label, rat.id, strain, behavior) %>%
  summarise(usv.count = length(label)) %>%
  filter(behavior == "Retrieval" |
           behavior == "Hover" |
           behavior == "LKP")

#total usv counts by behavior
beh.total.usv.counts <- beh.usv.counts.data %>% group_by(strain, file.name, rat.id, behavior) %>%
  summarise(tot.usv = sum(usv.counts)) %>% filter(behavior == "Retrieval" |
                                                    behavior == "Hover" |
                                                    behavior == "LKP")

lme.tot.usv.beh <- lme(usv.count ~ strain * behavior * label,
                       random = ~1|rat.id,
                       data = beh.usv.counts.data)
anova.lme(lme.tot.usv.beh)

lme.beh.tot.usv.sum <- summary(lsmeans(lme.tot.usv.beh, pairwise~strain*behavior|label, adjust = "Tukey"))[["contrasts"]]
View(lme.beh.tot.usv.sum[lme.beh.tot.usv.sum$p.value < 0.05,])

#usv counts by type by behavior
beh.usv.type.counts <- beh.usv.counts.data %>% group_by(strain, rat.id, behavior, label) %>%
  summarise(usv.count = sum(usv.count)) %>% filter(label == "flat" | label == "flat-mz" |
                                                      label == "flat-z" |
                                                      label == "short")
lme.usv.type.beh <- lme(usv.count ~ strain * behavior * label,
                        random = ~1|rat.id,
                        data = beh.usv.type.counts)
anova.lme(lme.usv.type.beh)
lme.beh.usv.type.sum <- summary(lsmeans(lme.usv.type.beh, pairwise~strain * behavior * label, adjust="Tukey"))[["contrasts"]]
View(lme.beh.usv.type.sum[lme.beh.usv.type.sum$p.value < 0.05,])



#before and after grouping with label
grouping.types.unmelted <- mbt.grouping.data %>% filter(label == "complex" | label == "flat" |
                                                          label == "flat-z" | label == "flat-mz" |
                                                          label == "short" | label == "short-c" |
                                                          label == "trill" | label == "trill-c") %>%
  group_by(strain, rat.id, file.name, score, label) %>%
  mutate(count.b4.group = ifelse(start.time < score, no.counts,0), 
         count.after.group = ifelse(start.time > score, no.counts, 0)) %>%
  summarise(b4.group = sum(count.b4.group), after.group = sum(count.after.group))

grouping.types <- melt(data = grouping.types.unmelted, id = c("strain","rat.id", "file.name", "score", "label"),
                                        variable.name="grouping.time",value.name = "usv.count")

lme.group.type <- lme(usv.count ~ strain * grouping.time * label,
                      random = ~1|rat.id,
                      data = grouping.types)
anova.lme(lme.group.type)

lme.group.type.sum <- summary(lsmeans(lme.group.type, pairwise~grouping.time * label, adjust="Tukey"))[["contrasts"]]
View(lme.group.type.sum[lme.group.type.sum$p.value < 0.05,])
