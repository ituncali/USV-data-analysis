##prepare count data
data_counts <- read.csv(file.choose(),stringsAsFactors = F)
source("C:/Users/ituncali/Documents/Master's Thesis/USV-Data-Management/src/count_total.R")
library(stringr)
library(dplyr)


momanesth_counts <- count_frame_2 %>% filter(recording == "MA")
momanesth_counts$strain <- as.factor(momanesth_counts$strain)

#kruskal-wallis
kruskal.test(total.counts ~ strain,
             data = momanesth_counts)

kruskal.test(total.counts ~ categories.allowed,
             data = momanesth_counts)

pairwise.wilcox.test(momanesth_counts$total.counts, momanesth_counts$categories.allowed,
                     p.adjust.method = "BH")




#momanesth counts lme
lme_momanesth <- lme(total.counts ~ strain * categories.allowed, 
                     random = ~1| rat.id,
                     data = momanesth_counts)
anova.lme(lme_momanesth)
lsmeans(lme_momanesth, pairwise ~ categories.allowed |strain, adjust = "Tukey")
ma_pairwise_summary <- (summary(lsmeans(lme_momanesth,pairwise~categories.allowed |strain, adjust = "Tukey")[["contrasts"]]))
View(ma_pairwise_summary[ma_pairwise_summary$p.value < 0.05,])


#compare momanesth counts to isolated pup counts
ma_5min <- start.rows %>% filter(recording=="MA" & start.time < 300)
ma_5min_frame <- ma_5min %>% group_by(label, rat.id) %>%
  summarise(total.counts = length(label))
expanded.ma.grid <- expand.grid(rat.id = unique(momanesth_counts$rat.id),
                                label = unique(momanesth_counts$categories.allowed))
ma_5min_countframe <- left_join(expanded.ma.grid, ma_5min_frame)  
ma_5min <- ma_5min_countframe %>% 
  mutate(total.counts = ifelse(is.na(total.counts),0,total.counts/8)) %>%
  left_join(file.name.key[file.name.key$recording=="MA",]) %>% 
  select(-file.name)
pup_to_join <- pup_counts %>% select(label = categories.allowed,
                                     total.counts, file.name, strain,
                                     rat.id, recording) %>%
  group_by(label, rat.id,strain) %>%
  summarise(total.counts = mean(total.counts)) %>%
  mutate(recording=c(rep("Pup",length(total.counts))))

ma_pup_counts <- rbind.data.frame(ma_5min, pup_to_join)

lme_ma_pup <- lme(total.counts ~ strain * label * recording,
                  random = ~1|rat.id/recording/label,
                  data=ma_pup_counts)
anova.lme(lme_ma_pup)
ma.pup.sum <- summary(lsmeans(lme_ma_pup, pairwise~recording|strain*label,
                              adjust="Tukey"))[["contrasts"]]
View(ma.pup.sum[ma.pup.sum$p.value<0.05,])

ma_pup_per <- ma_pup_counts %>% ungroup() %>%
  filter(label=="flat"|label=="flat-z"|
           label=="flat-mz"|label=="short"|
           label=="short-c"|label=="short-ur") %>%
  group_by(strain, rat.id, recording) %>% 
  mutate(total.filecounts = sum(total.counts),
           percent = ifelse(total.filecounts>0,total.counts/total.filecounts,0)) 
  

lme_mapup_per <- lme(percent ~ strain * label * recording,
                     random = ~1|rat.id/recording/label,
                     data=ma_pup_per)
anova.lme(lme_mapup_per)
ma.pup.per.sum <- summary(lsmeans(lme_mapup_per, pairwise~recording|strain*label,
                              adjust="Tukey"))[["contrasts"]]
View(ma.pup.per.sum[ma.pup.per.sum$p.value<0.05,])

##prepare freq data
data_freqs <- read.csv(file.choose(), stringsAsFactors = F)
m.freq <- apply(data_freqs[,c(1:50)],1,function(y) mean(y))
data_freqs <- mutate(data_freqs, m.freq = m.freq)

momanesth_freqs <- data_freqs %>% filter(recording == "MA") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label == "flat" | label == "flat-z" | label == "short" | label == "short-c" | label == "short-ur")

#lme for mom anesth freqs
lme_ma_freq <- lme(mean.freq ~ strain * label, 
                   random = ~1|rat.id, 
                   data = momanesth_freqs)
anova.lme(lme_ma_freq)
lsmeans(lme_ma_freq, pairwise ~ label, adjust = "Tukey")
ma_pairwise_freq <- (summary(lsmeans(lme_ma_freq,pairwise~label * strain, adjust = "Tukey")[["contrasts"]]))
View(ma_pairwise_freq[ma_pairwise_freq$p.value < 0.05,])

#compare freqs between isolated pups and MA
ma_freqs_comp <- momanesth_freqs %>% filter(label=="flat"|label=="flat-z"|
                                              label=="short")

pup_freqs_comp <- data_freqs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label == "flat" | label == "flat-z" |label == "short")

ma.pup.freqs <- rbind.data.frame(ma_freqs_comp, pup_freqs_comp)

lme.ma.pup.freqs <- lme(mean.freq ~ strain * recording * label,
                        random = ~1|rat.id,
                        data = ma.pup.freqs)

##prepare dur data
data_durs <- read.csv(file.choose(),stringsAsFactors = F)

momanesth_durs <- data_durs %>% filter(recording == "MA") %>% 
  group_by(strain, label, file.name,recording, rat.id) %>% 
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label == "flat" | label == "flat-z" | label == "short" | label == "short-c" | label == "short-ur")

#lme for mom anesth durs
lme_ma_dur <- lme(mean.dur ~ strain * label, 
                  random = ~1|rat.id, 
                  data = momanesth_durs)
anova.lme(lme_ma_dur)
lsmeans(lme_ma_dur, pairwise ~ label * strain, adjust = "Tukey")
ma_pairwise_dur <- (summary(lsmeans(lme_ma_dur,pairwise~label * strain, adjust = "Tukey")[["contrasts"]]))
View(ma_pairwise_freq[ma_pairwise_dur$p.value < 0.05,])

#compare pup iso and ma durs
ma_durs_comp <- momanesth_durs %>% filter(label=="flat"|label=="flat-z"|
                                            label=="short")
pup_durs_comp <- data_durs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label == "flat" | label == "flat-z" |label == "short")

ma.pups.durs <- rbind.data.frame(ma_durs_comp, pup_durs_comp)

lme.ma.pup.durs <- lme(mean.dur ~ strain*recording*label,
                       random=~1|rat.id,
                       data=ma.pups.durs)

##pup litter weights
litter.w <- read.csv(file.choose(),stringsAsFactors = F)

litter.w.lme <- lme(bw.pup ~ strain, random = ~1|rat.id, data = litter.w)
anova.lme(litter.w.lme)

#momanesth with litter weights
momanesth_counts_w <- left_join(momanesth_counts, litter.w)

momanesth.litter.w.lme.data <- momanesth_counts_w %>% group_by(strain, rat.id) %>%
  summarise(total.count.an = sum(total.counts), bw.pup = unique(bw.pup))

lme.ma.lme <- lme(total.count.an ~ strain * bw.pup, random = ~1|rat.id, data = momanesth.litter.w.lme.data)
anova.lme(lme.ma.lme)

library(ggpubr)
ggscatter(momanesth.litter.w.lme.data, x = "bw.pup", y = "total.count.an", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "pup weight (g)", ylab = "USVs (n)")
cor.test(momanesth.litter.w.lme.data$bw.pup, momanesth.litter.w.lme.data$total.count.an, method = "pearson")


#freqs corr w/ weight
momanesth_freqs_w <- left_join(momanesth_freqs, litter.w)
lme.ma.freqs <- lme(m.freq ~ strain * bw.pup, random = ~1|rat.id, data = momanesth_freqs_w)
anova.lme(lme.ma.freqs)

momanesth_avg_freqs <- data_freqs %>% filter(recording == "MA") %>% 
  group_by(strain, file.name, rat.id) %>% 
  summarise(mean.freq.an = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq))
momanesth_avg_freqs_w <- left_join(momanesth_avg_freqs, litter.w)


ggscatter(momanesth_freqs_w, x = "bw.pup", y = "mean.freq", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "pup weight (g)", ylab = "mean frequency (Hz)")
cor.test(momanesth_freqs_w$bw.pup, momanesth_freqs_w$m.freq, method = "pearson")


momanesth_durs_w <- left_join(momanesth_durs,litter.w)
lme.ma.durs <- lme(mean.dur ~ strain * bw.pup, random = ~1|rat.id, data = momanesth_durs_w)
anova.lme(lme.ma.durs)

ggscatter(momanesth_durs_w, x = "bw.pup", y = "mean.dur", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "pup weight (g)", ylab = "duration (s)")
cor.test(momanesth_durs_w$bw.pup, momanesth_durs_w$mean.dur, method = "pearson")

momanesth_avg_durs <- data_durs %>% filter(recording == "MA") %>% 
  group_by(strain, file.name, rat.id) %>% 
  summarise(mean.dur.an = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration))
momanesth_avg_durs_w <- left_join(momanesth_avg_durs,litter.w)



##momanesth before and after grouping
#want to compare 1st 3 min, 3rd 3 min and 5th 3 min...
ma.xtra <- start.rows %>% filter(recording == "MA")

three.times.data <- ma.xtra %>% 
  mutate(time.group = ifelse(start.time < 180, "1",
                             ifelse(start.time >= 420 & start.time < 540, "2",
                                    ifelse(start.time > 720, "3", NA)))) %>%
  filter(!is.na(time.group)) %>%
  group_by(strain, time.group, rat.id, file.name, label) %>%
  summarise(time.count = length(label))

to.join <- expand.grid(label=unique(momanesth_counts$categories.allowed), 
                       time.group=c("1","2","3"), 
                       file.name=unique(momanesth_counts$file.name))
to.join <- left_join(to.join, file.name.key)

pru <- to.join %>% left_join(three.times.data) %>%
  mutate(time.count = ifelse(is.na(time.count), 0, time.count))

three.times.lme <- lme(time.count ~ strain * time.group * label, 
                         random = ~1|rat.id/time.group/label, 
                         data = pru)
anova.lme(three.times.lme)

ma_summary <- (summary(lsmeans(three.times.lme,pairwise ~ (strain * time.group)|label, adjust = "Tukey")[["contrasts"]]))
View(ma_summary[ma_summary$p.value < 0.05,])

pru2 <- pru %>% filter(label=="flat"| label=="flat-z"|
                       label=="short"|label=="short-ur"|
                       label=="short-c")

three.times.profile <- pru2 %>% group_by(strain, rat.id, time.group) %>% 
  mutate(tot.count.in.time = sum(time.count)) %>% 
  group_by(strain, rat.id, time.group, label) %>% 
  mutate(rel.call.count = time.count/tot.count.in.time)
three.times.profile.lme <- lme(rel.call.count ~ strain * time.group * label, 
                               random = ~1|rat.id/time.group/label, 
                               data = three.times.profile)
anova.lme(three.times.profile.lme)
fart.sum <- summary(lsmeans(three.times.profile.lme, pairwise~strain|label))[["contrasts"]]
View(fart.sum[fart.sum$p.value<0.05,])

##momanesth before and during grouping
group.induced.data <- ma.xtra %>% 
  mutate(time.group = ifelse(start.time < 600 & start.time > 570, "b4",
                             ifelse(start.time >= 600 & start.time < 630, "during", NA))) %>%
  filter(!is.na(time.group)) %>%
  group_by(strain, time.group, rat.id, file.name, label) %>%
  summarise(time.count = length(label))
group.induced.data$file.name <- as.factor(as.character(group.induced.data$file.name))

to.join <- expand.grid(label=unique(momanesth_counts$categories.allowed), 
                       time.group=c("b4","during"), 
                       file.name=unique(momanesth_counts$file.name))
to.join <- left_join(to.join, file.name.key)

all.levels <- to.join %>% left_join(group.induced.data) %>%
  mutate(time.count = ifelse(is.na(time.count), 0, time.count))

all.levels$strain <- as.factor(all.levels$strain)
all.levels$time.group <- as.factor(all.levels$time.group)
all.levels$label <- as.factor(all.levels$label)


group.induced.lme <- lme(time.count ~ strain * time.group * label, 
                         random = ~1|rat.id/time.group/label, 
                         data = all.levels)
anova.lme(group.induced.lme)
ma_summary <- (summary(lsmeans(group.induced.lme,pairwise ~ time.group | strain*label, adjust = "Tukey")[["contrasts"]]))
View(ma_summary[ma_summary$p.value < 0.05,])

all.levs.profile <- all.levels %>% group_by(strain, rat.id, time.group) %>% 
  mutate(tot.count.in.time = sum(time.count)) %>% 
  group_by(strain, rat.id, time.group, label) %>% 
  mutate(rel.call.count = ifelse(tot.count.in.time >0, time.count/tot.count.in.time, 0))

group.profile.lme <- lme(rel.call.count ~ strain * time.group * label, 
                               random = ~1|rat.id/time.group/label, 
                               data = all.levs.profile)
anova.lme(group.profile.lme)
ma_summary <- (summary(lsmeans(group.profile.lme,pairwise ~ (time.group * label) | strain, adjust = "Tukey")[["contrasts"]]))
View(ma_summary[ma_summary$p.value < 0.05,])

