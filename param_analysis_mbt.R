data_counts <- read.csv(file.choose(),stringsAsFactors = F)
source("C:/Users/ituncali/Documents/Master's Thesis/USV-Data-Management/src/count_total.R")
library(stringr)
library(dplyr)

mbt_counts <- count_frame_2 %>% filter(recording == "MBT")

##counts analysis
library(nlme)
lme_mbt <- lme(total.counts ~ strain * categories.allowed, 
                random = ~1| rat.id/categories.allowed,
                data = mbt_counts)
anova.lme(lme_mbt)

library(lsmeans)
mbt.lme.sum <- summary(lsmeans(lme_mbt, pairwise ~ strain|categories.allowed, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.sum[mbt.lme.sum$p.value < 0.05,])

#so we're going to use: flat, flat-z, flat-mz, short and trill
main.cats <- c("flat","flat-z","flat-mz","short","trill")

#by proportion only include flat, flat-z, flat-mz, short and trill
mbt.per <- mbt_counts %>% filter(categories.allowed %in% main.cats) %>%
  group_by(strain, rat.id, file.name) %>% 
  mutate(usv.per = total.counts/sum(total.counts))
mbt.per.lme <- lme(usv.per~strain*categories.allowed, 
                   random=~1|rat.id/categories.allowed, 
                   data=mbt.per)
anova.lme(mbt.per.lme)
mbt.per.sum <- summary(lsmeans(mbt.per.lme, pairwise~strain|categories.allowed, adjust="Tukey"))[["contrasts"]]
View(mbt.per.sum[mbt.per.sum$p.value<0.05,])

#frequency data

mbt_freqs <- data_freqs %>% filter(recording == "MBT") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label == "flat" | label == "flat-z" | label=="flat-mz" | 
           label == "short" | label=="trill")

#shapiro wilks test of normality
mbt_freqs_sw <- data_freqs %>% filter(recording == "MBT") %>% 
  select(strain, duration, m.freq, label)
stargazer(mbt_freqs_sw, type = "text")

shapiro.test(log(pup_freqs_sh[mbt_freqs_sw$strain == "SD" & mbt_freqs_sw$label == "flat",]$m.freq))


##frequency analysis
lme_mbt_freq <- lme(mean.freq ~ strain * label, 
                    random = ~1|rat.id/label, 
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
  filter(label == "flat" | label == "flat-z" | label=="flat-mz" | 
           label == "short" | label=="trill")


##duration analysis
lme_mbt_dur <- lme(mean.dur ~ strain * label,
                   random = ~1|rat.id/label,
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

wilcox.test(mbt_behavior_unmelted$retrieval ~ mbt_behavior_unmelted$strain)

wilcox.test(mbt_behavior_unmelted$corporal ~ mbt_behavior_unmelted$strain)

wilcox.test(mbt_behavior_unmelted$anogenital ~ mbt_behavior_unmelted$strain)

#behavioral latencies
mbt.beh.lats.data <- mbt_behavior %>% filter(behavior == "lat.retrieve" |
                                               behavior == "lat.group" |
                                               behavior == "lat.hover" |
                                               behavior == "lat.nurse")

wilcox.test(mbt_behavior_unmelted$lat.retrieve ~ mbt_behavior_unmelted$strain)

wilcox.test(mbt_behavior_unmelted$lat.group ~ mbt_behavior_unmelted$strain)

wilcox.test(mbt_behavior_unmelted$lat.hover ~ mbt_behavior_unmelted$strain)

wilcox.test(mbt_behavior_unmelted$lat.nurse ~ mbt_behavior_unmelted$strain)

#behavioral durations
mbt.beh.durs.data <- mbt_behavior %>% filter(behavior == "dur.hover" |
                                               behavior == "dur.nurse")

wilcox.test(mbt_behavior_unmelted$dur.hover ~ mbt_behavior_unmelted$strain)

wilcox.test(mbt_behavior_unmelted$dur.nurse ~ mbt_behavior_unmelted$strain)

#usvs during behaviors
beh_times_data <- read.csv("data/Exp1_MBT_beh_times.csv",stringsAsFactors = F)
names(beh_times_data) <- c("rat.id","file.name","strain","behavior","start","end","beh.dur")
beh_times_data$behavior <- gsub("HKP","LKP",beh_times_data$behavior)
breaks <- c(rbind(beh_times_data$start,beh_times_data$end))


merged.data <- merge(start.rows, beh_times_data, by = c("file.name","rat.id","strain"))

#USV counts per category by behavior

beh.total.usv.counts <- merged.data %>% 
  filter(start.time > start & start.time < end) %>%
  group_by(strain, behavior, label) %>%
  summarise(usv.count = length(label)) %>% filter(behavior == "Retrieval" |
                                                    behavior == "Hover" |
                                                    behavior == "LKP")
beh.label.expand <- expand.grid(file.name = unique(mbt_counts$file.name),
                                behavior = c("Retrieval","Hover","LKP"),
                                label = unique(mbt_counts$categories.allowed))
usv.beh.counts.filenames <- left_join(beh.label.expand, beh.total.usv.counts)
usv.beh.counts.join <- left_join(usv.beh.counts.filenames, file.name.key)
usv.beh.counts <- usv.beh.counts.join %>% 
  mutate(usv.count = ifelse(is.na(usv.count),0,usv.count))

lme.tot.usv.beh <- lme(usv.count ~ strain * behavior * label,
                       random = ~1|rat.id/behavior/label,
                       data = usv.beh.counts)
anova.lme(lme.tot.usv.beh)

lme.beh.tot.usv.sum <- summary(lsmeans(lme.tot.usv.beh, pairwise~strain*behavior|label, adjust = "Tukey"))[["contrasts"]]
View(lme.beh.tot.usv.sum[lme.beh.tot.usv.sum$p.value < 0.05,])

#usv percents by type by behavior
usv.beh.per <- usv.beh.counts %>% 
  filter(label=="flat" | label == "flat-z" | label=="flat-mz" |
           label == "short" | label == "trill") %>%
  group_by(file.name, rat.id, behavior) %>%
  mutate(percent = ifelse(usv.count>0, usv.count/sum(usv.count), 0))

lme.type.per <- lme(percent ~ strain * behavior * label,
                    random = ~1|rat.id/behavior/label,
                    data = usv.beh.per)
anova.lme(lme.type.per)
type.per.sum <- summary(lsmeans(lme.type.per, pairwise~strain *behavior|label, adjust="Tukey"))[["contrasts"]]
View(type.per.sum[type.per.sum$p.value < 0.05,])

#see if duration changes by behavior
dur.merged <- merge(data_durs, beh_times_data, by = c("file.name","rat.id","strain"))
beh.dur.data <- dur.merged %>% ungroup() %>%
  filter(start.time > start & start.time < end) %>%
  group_by(file.name, rat.id, strain, behavior, label) %>%
  summarise(m.dur = mean(duration)) %>%
  filter((behavior == "Retrieval" |behavior == "Hover" |behavior == "LKP") & 
           (label=="flat"|label=="flat-z"|label=="short"))
beh.dur.lme <- lme(m.dur ~ strain * behavior * label,
                   random = ~1|rat.id/behavior/label,
                   data = beh.dur.data)
anova.lme(beh.dur.lme)

#see if frequency changes by behavior
freq.merged <- merge(data_freqs, beh_times_data, by = c("file.name","rat.id","strain"))
beh_freq_data <- freq.merged %>% ungroup() %>%
  filter(start.time > start & start.time < end) %>%
  group_by(file.name, rat.id, strain, behavior, label) %>%
  summarise(mean.freq = mean(m.freq)) %>%
  filter((behavior == "Retrieval" |behavior == "Hover" |behavior == "LKP") & 
           (label=="flat"|label=="flat-z"|label=="short"))

beh.freq.lme <- lme(mean.freq ~ strain * behavior * label,
                    random = ~1|rat.id/behavior/label,
                    data = beh_freq_data)
anova.lme(beh.freq.lme)

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


###compare mbt to ma!!
mbt_15 <- start.rows %>% filter(recording=="MBT" & start.time<900) %>%
  group_by(strain, rat.id, label) %>%
  summarise(mbt.count = length(label))

ma_all <- momanesth_counts %>% ungroup() %>%
  select(label = categories.allowed,
          ma.count = total.counts,
          strain, rat.id)

mbt.ma.join <- left_join(ma_all, mbt_15) %>%
  mutate(mbt.count = ifelse(is.na(mbt.count), 0, mbt.count))


mbt.ma.data <- melt(mbt.ma.join, id.vars = c("strain","rat.id","label"),
                    variable.name = "recording", value.name = "usv.count")

mbt.ma.lme <- lme(usv.count~strain*label*recording,
                  random = ~1|rat.id/recording/label,
                  data = mbt.ma.data)
anova.lme(mbt.ma.lme)
mbt.ma.sum <- summary(lsmeans(mbt.ma.lme, pairwise~strain*recording|label, adjust="Tukey"))[["contrasts"]]
View(mbt.ma.sum[mbt.ma.sum$p.value < 0.05,])

#percentages
mbt.ma.per <- mbt.ma.data %>% 
  filter(label == "flat" | label=="flat-z"|label=="flat-mz"|
           label=="short"|label=="short-c"|label== "short-ur"|
           label=="trill") %>%
  group_by(strain, rat.id, recording) %>%
  mutate(usv.per = usv.count/sum(usv.count))
mbt.ma.per.lme <- lme(usv.per~strain*label*recording,
                  random = ~1|rat.id/recording/label,
                  data = mbt.ma.per)
anova.lme(mbt.ma.per.lme)
mbt.ma.per.sum <- summary(lsmeans(mbt.ma.per.lme, pairwise~strain*recording|label, adjust="Tukey"))[["contrasts"]]
View(mbt.ma.per.sum[mbt.ma.per.sum$p.value < 0.05,])


#totals by strain
total.strain.mbt.ma.data <- mbt.ma.data %>% 
  group_by(recording, strain, rat.id) %>% summarise(tot.count = sum(usv.count))

anova.lme(lme(tot.count~strain*recording, random=~1|rat.id, 
              data=total.strain.mbt.ma.data))
lme.tot.mbt.ma <- summary(lsmeans(lme(tot.count~strain*recording, random=~1|rat.id, 
                                      data=total.strain.mbt.ma.data), pairwise~recording * strain, adjust="Tukey"))[["contrasts"]]
View(lme.tot.mbt.ma[lme.tot.mbt.ma$p.value < 0.05,])






##distinguish between mom and pup
param.cutoffs <- data.frame(strain = table.freqs$strain,
                               label = table.freqs$label,
                            rat.id = table.freqs$rat.id,
                               freq.min = table.freqs$Min*1000,
                               freq.max = table.freqs$Max*1000,
                               dur.min = table.durs$Min/1000,
                               dur.max = table.durs$Max/1000)

mom_calls_1 <- data_freqs %>% filter(recording=="MBT") %>% 
  left_join(param.cutoffs) %>%
  filter(m.freq < freq.min | m.freq > freq.max | duration < dur.min |
           duration > dur.max) %>%
  select(-freq.min, -freq.max, -dur.min, -dur.max)
mom_calls_2 <- data_freqs %>% filter(recording=="MBT" & (label=="trill"|
                                                           label=="trill-c"))

mom_freqs_joined <- rbind.data.frame(mom_calls_1, mom_calls_2)  
#frequency data of maternal USVs
mom_freqs <- mom_freqs_joined %>% group_by(strain, rat.id, label) %>%
  summarise(mean.freq = mean(m.freq), sem.freq = sd(m.freq)/sqrt(length(m.freq)), 
            count = length(label))
mom.freqs.lme <- lme(mean.freq ~ strain * label,
                     random = ~1|rat.id/label,
                     data = mom_freqs)
anova.lme(mom.freqs.lme)

#duration data of maternal usvs
mom_durs_1 <- data_durs %>% filter(recording=="MBT") %>%
  anti_join(mom_freqs_joined, by=c("start.time"))
mom_durs_2 <- mom_durs_1 %>%
  left_join(param.cutoffs) %>%
  filter(duration < dur.min |duration > dur.max) %>%
  select(-freq.min, -freq.max, -dur.min, -dur.max, -unique.id)
mom_freqs_joined_fix <- mom_freqs_joined %>%
  select(file.name, label, strain, rat.id, recording, duration, start.time)

mom_durs <- rbind.data.frame(mom_durs_2, mom_freqs_joined_fix) 
mom_durs_data <- mom_durs %>% group_by(strain, rat.id, label) %>%
  summarise(mean.dur = mean(duration), count = length(duration))

mom.durs.lme <- lme(mean.dur ~ strain * label,
                    random = ~1|rat.id/label,
                    data=mom_durs_data)
anova.lme(mom.durs.lme)
mom.durs.sum <- summary(lsmeans(mom.durs.lme, pairwise~strain|label, adjust="Tukey"))[["contrasts"]]
View(mom.durs.sum[mom.durs.sum$p.value<0.05,])

#count data of maternal USVs
mom_calls_1_adjust <- mom_durs %>% select(-duration)

mom_calls_3 <- start.rows %>% filter(recording=="MBT") %>%
  anti_join(mom_durs, by = c("start.time")) %>%
  filter(label=="trill"|label=="trill-c") %>%
  select(-unique.id)

mom_counts_1 <- rbind.data.frame(mom_calls_1_adjust, mom_calls_3) 
mom_counts_2 <- mom_counts_1 %>% group_by(strain, label, rat.id) %>%
  summarise(usv.count= length(label))

grid <- expand.grid(rat.id = unique(mom_counts_1$rat.id),
                    label = unique(mom_counts_1$label))
mom_counts_3 <- left_join(grid, mom_counts_2)
mom_counts <- mom_counts_3 %>% mutate(strain = ifelse(str_detect(rat.id,"SD")==T,
                                                      "SD","WK"),
                                      usv.count = ifelse(is.na(usv.count),0,usv.count))

mom.counts.lme <- lme(usv.count ~ strain * label,
                      random = ~1|rat.id/label,
                      data = mom_counts)
anova.lme(mom.counts.lme)
mom.counts.sum <- summary(lsmeans(mom.counts.lme, pairwise~strain|label, adjust="Tukey"))[["contrasts"]]
View(mom.counts.sum[mom.counts.sum$p.value<0.05,])

#start times!!
mom_start_1 <- mom_counts_1


##rough maternal retrieval
r.ret <- data.frame(rat.id = c("WK77","SD78","SD76","WK95"),
                    no.rough = c(3,2,2,4))

pru <- left_join(r.ret, mbt.beh.counts.data[mbt.beh.counts.data$behavior=="retrieval",])
fart <- mbt.trills %>% select(label, usv.per,
                              strain, rat.id) %>%
  filter(label=="short")

r.ret.per <- pru %>% mutate(percent = no.rough/score) %>%
  left_join(fart)
r.ret.lme <- lme(percent ~ strain,
                 random=~1|rat.id,
                 data=r.ret.per)
anova(r.ret.lme)

##see if maternal trills are correlated with behaviors
all.trills <- mom_counts %>% filter(label=="trill"|label=="trill-c") %>% 
  group_by(rat.id, strain) %>% summarise(trill.count = sum(usv.count))
trill.group.corr <- left_join(all.trills, mbt_behavior_unmelted) %>%
  mutate(dur.retrieve = lat.group - lat.retrieve)

cor.test(trill.group.corr$trill.count, trill.group.corr$lat.group)

cor.test(trill.group.corr$trill.count, trill.group.corr$retrieval)

cor.test(trill.group.corr$trill.count, trill.group.corr$corporal)

cor.test(trill.group.corr$trill.count, trill.group.corr$anogenital)

cor.test(trill.group.corr$trill.count, trill.group.corr$nest.building)

cor.test(trill.group.corr$trill.count, trill.group.corr$lat.retrieve)

cor.test(trill.group.corr$trill.count, trill.group.corr$lat.hover)

cor.test(trill.group.corr$trill.count, trill.group.corr$lat.nurse)

cor.test(trill.group.corr$trill.count, trill.group.corr$dur.hover)

cor.test(trill.group.corr$trill.count, trill.group.corr$dur.nurse)

cor.test(trill.group.corr$trill.count, trill.group.corr$crossing)

cor.test(trill.group.corr$trill.count, trill.group.corr$rearing)

cor.test(trill.group.corr$trill.count, trill.group.corr$self.groom)

cor.test(trill.group.corr$trill.count, trill.group.corr$dur.retrieve)



