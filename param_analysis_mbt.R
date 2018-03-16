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

mbt_counts <- count_frame_2 %>% filter(recording == "MBT")

##counts analysis
library(nlme)
lme_mbt <- lme(total.counts ~ strain * categories.allowed, 
                random = ~1| rat.id,
                data = mbt_counts)
anova.lme(lme_mbt)

library(lsmeans)
mbt.lme.sum <- summary(lsmeans(lme_mbt, pairwise ~ strain * categories.allowed, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.sum[mbt.lme.sum$p.value < 0.05,])

#frequency data
data_freqs <- read.csv(file.choose(), stringsAsFactors = F)
m.freq <- apply(data_freqs[,c(1:50)],1,function(y) mean(y))
data_freqs <- mutate(data_freqs, m.freq = m.freq)
mbt_freqs <- data_freqs %>% filter(recording == "MBT") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short" | label == "short-c" | label == "complex" | label == "trill" | label == "trill-c")


##frequency analysis
lme_mbt_freq <- lme(mean.freq ~ strain * label, 
                    random = ~1|rat.id, 
                    data = mbt_freqs)
anova.lme(lme_mbt_freq)
mbt.lme.freq.sum <- summary(lsmeans(lme_mbt_freq, pairwise ~ strain * label, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.freq.sum[mbt.lme.freq.sum$p.value < 0.05,])


#duration data
mbt_durs <- data_durs %>% filter(recording == "MBT") %>% group_by(strain, label,recording,file.name,rat.id) %>%
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short" | label == "short-c" | label == "complex" | label == "trill" | label == "trill-c")


##duration analysis
lme_mbt_dur <- lme(mean.dur ~ strain * label,
                   random = ~1|rat.id,
                   data = mbt_durs)
anova.lme(lme_mbt_dur)
mbt.lme.dur.sum <- summary(lsmeans(lme_mbt_dur, pairwise~strain*label, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.dur.sum[mbt.lme.dur.sum$p.value < 0.05,])


##maternal behavior data and analysis
mbt_behavior_unmelted <- read.csv(file.choose(),stringsAsFactors = F)
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
mbt.lme.beh.lats.sum <- summary(lsmeans(lme_mbt_beh_lats, pairwise~strain * behavior, adjust = "Tukey"))[["contrasts"]]
View(mbt.lme.beh.lats.sum[mbt.lme.beh.lats.sum$p.value < 0.05,])

#behavioral durations
mbt.beh.durs.data <- mbt_behavior %>% filter(behavior == "dur.hover" |
                                               behavior == "dur.nurse")

lme_mbt_beh_durs <- lme(score ~ strain * behavior,
                        random = ~1|rat.id,
                        data = mbt.beh.durs.data)
anova.lme(lme_mbt_beh_durs)


#usvs during behaviors
beh_times_data <- read.csv(file.choose(),stringsAsFactors = F)
names(beh_times_data) <- c("rat.id","file.name","strain","behavior","start","end","beh.dur")
breaks <- c(rbind(beh_times_data$start,beh_times_data$end))

allowed.categories <- c("flat", "flat-z", "flat-mz", "short", "short-su", "short-sd",
                        "short-ur", "short-dr", "short-c", "complex", "upward ramp",
                        "downward ramp", "step up", "step down", "multi-step", "multi-step-s",
                        "trill", "trill-c", "trill-f", "inverted-U", "unclear")

categories.allowed.search <- paste0(allowed.categories, "(?!-)")
library(stringr)
start_data <- cbind(data_counts, 
                    no.counts = mapply(function(q) sum(str_count(q,categories.allowed.search)),data_counts$label))

merged.data <- merge(start_data, beh_times_data, by = c("file.name","rat.id","strain"))

beh.usv.counts.data <- merged.data %>% group_by(file.name, strain, start.time, behavior, start, end) %>%
  mutate(usv.counts = ifelse(start.time > start && start.time < end, no.counts, 0))

#total usv counts by behavior
beh.total.usv.counts <- beh.usv.counts.data %>% group_by(strain, file.name, rat.id, behavior) %>%
  summarise(tot.usv = sum(usv.counts)) %>% filter(behavior == "Retrieval" |
                                                    behavior == "Hover" |
                                                    behavior == "LKP")

lme.tot.usv.beh <- lme(tot.usv ~ strain * behavior,
                       random = ~1|rat.id,
                       data = beh.total.usv.counts)
anova.lme(lme.tot.usv.beh)

lme.beh.tot.usv.sum <- summary(lsmeans(lme.tot.usv.beh, pairwise~strain * behavior, adjust = "Tukey"))[["contrasts"]]
View(lme.beh.tot.usv.sum[lme.beh.tot.usv.sum$p.value < 0.05,])

#usv counts by type by behavior
beh.usv.type.counts <- beh.usv.counts.data %>% group_by(strain, file.name, rat.id, behavior, label) %>%
  summarise(usv.count = sum(usv.counts)) %>% filter((behavior == "Retrieval" |behavior == "Hover" |
                                                       behavior == "LKP") &
                                                      (label == "flat" | label == "flat-mz" |
                                                      label == "flat-z" |
                                                      label == "complex" |
                                                      label == "short" |
                                                      label == "short-c" |
                                                      label == "trill" |
                                                      label == "trill-c"))
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

