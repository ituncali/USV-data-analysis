##prepare count data
data_counts <- read.csv(file.choose(),stringsAsFactors = F)
source("C:/Users/ituncali/Documents/Master's Thesis/USV-Data-Management/src/count_total.R")
library(stringr)
library(dplyr)
allowed.categories <- c("flat", "flat-z", "flat-mz", "short", "short-su", "short-sd",
                        "short-ur", "short-dr", "short-c", "complex", "upward ramp",
                        "downward ramp", "step up", "step down", "multi-step", "multi-step-s",
                        "trill", "trill-c", "trill-f", "inverted-U", "unclear") 
count_list <- by(data = data_counts, INDICES = data_counts$file.name, FUN = function(x) count_total(x, allowed.categories))
count_frame <- do.call(rbind, count_list)
count_frame$file.name <- row.names(count_frame)
row.names(count_frame) <- NULL
count_frame$file.name <- str_replace(count_frame$file.name, pattern = ".[0-9]+$", replacement =  "" )
count_frame <- count_frame %>% group_by(file.name) %>%
  mutate(total.filecounts=sum(total.counts), rel.filecount=total.counts/total.filecounts)

library(reshape2)
file.name.key.unmelted <- read.csv(file.choose(),stringsAsFactors = F)

file.name.key <- melt(data=file.name.key.unmelted,id = c("strain","rat.id"),
                      variable.name="recording",value.name = "file.name")

library(dplyr)
count_frame_2 <- left_join(count_frame,file.name.key)

momanesth_counts <- count_frame_2 %>% filter(recording == "MA")
#momanesth counts lme
lme_momanesth <- lme(total.counts ~ strain * categories.allowed, 
                     random = ~1| rat.id,
                     data = momanesth_counts)
anova.lme(lme_momanesth)
lsmeans(lme_momanesth, pairwise ~ categories.allowed * strain, adjust = "Tukey")
ma_pairwise_summary <- (summary(lsmeans(lme_momanesth,pairwise~categories.allowed * strain, adjust = "Tukey")[["contrasts"]]))
View(ma_pairwise_summary[ma_pairwise_summary$p.value < 0.05,])


##prepare freq data
data_freqs <- read.csv(file.choose(), stringsAsFactors = F)
m.freq <- apply(data_freqs[,c(1:50)],1,function(y) mean(y))
data_freqs <- mutate(data_freqs, m.freq = m.freq)

momanesth_freqs <- data_freqs %>% filter(recording == "MA") %>% 
  group_by(strain, label, file.name, rat.id) %>% 
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


##prepare dur data
data_durs <- read.csv(file.choose(),stringsAsFactors = F)

momanesth_durs <- data_durs %>% filter(recording == "MA") %>% 
  group_by(strain, label, file.name, rat.id) %>% 
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

momanesth_freqs_w <- left_join(data_freqs[data_freqs$recording=="MA",], litter.w)
lme.ma.freqs <- lme(m.freq ~ strain * bw.pup, random = ~1|rat.id, data = momanesth_freqs_w)
anova.lme(lme.ma.freqs)

momanesth_durs_w <- left_join(data_durs[data_durs$recording=="MA",],litter.w)
lme.ma.durs <- lme(duration ~ strain * bw.pup, random = ~1|rat.id, data = momanesth_durs_w)
anova.lme(lme.ma.durs)

##momanesth time bins
#by 3 min bins
bins.ma <- cut(start_data[start_data$recording == "MA",]$start.time,5,include.lowest=T, labels = c("1","2","3","4","5"))

bin.data <- start_data %>% filter(recording == "MA") %>% 
  mutate(bin = bins.ma) %>% 
  group_by(bin,strain,file.name) %>% 
  summarise(count = sum(no.counts))

bin.ma.lme <- lme(count ~ strain * bin, random = ~1|file.name, data = bin.data)
anova.lme(bin.ma.lme)


##momanesth before and after grouping
grouping.data <- start_data %>% filter(recording == "MA") %>%
  mutate(time.group = ifelse(start.time < 660, "b4","aft")) %>%
  group_by(strain, time.group, rat.id, file.name) %>%
  summarise(m.count = sum(no.counts))

grouping.lme <- lme(m.count ~ strain * time.group, random = ~1|rat.id, data = grouping.data)
anova.lme(grouping.lme)


