mbt_counts_all <- count_frame_2 %>% filter(recording == "MBT")
mbt_counts_all$strain <- as.factor(mbt_counts_all$strain)

mbt.per.1 <- mbt_counts_all %>% 
       group_by(strain, categories.allowed) %>% 
       summarise(usv.count = sum(total.counts)) %>%
       group_by(strain) %>% mutate(percent=usv.count/sum(usv.count)*100) %>%
  filter(percent>2)

label.to.keep <- as.vector(unique(mbt.per.1$categories.allowed))

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


##USV counts during maternal behaviors
beh_times_data <- read.csv("data/Exp1_MBT_beh_times.csv",stringsAsFactors = F)
names(beh_times_data) <- c("rat.id","file.name","strain","behavior","start","end","beh.dur")
beh_times_data$behavior <- gsub("HKP","LKP",beh_times_data$behavior)
breaks <- c(rbind(beh_times_data$start,beh_times_data$end))
merged.data <- merge(start.rows, beh_times_data, by = c("file.name","rat.id","strain"))

beh.total.usv.counts <- merged.data %>% 
  filter(start.time > start & start.time < end) %>%
  group_by(strain, behavior, rat.id, label) %>%
  summarise(usv.count = length(label)) %>% filter(behavior == "Retrieval" |
                                                    behavior == "Hover" |
                                                    behavior == "LKP")
beh.total.pers <- beh.total.usv.counts %>% 
  group_by(strain, behavior, label) %>%
  summarise(usv.count = sum(usv.count)) %>%
  group_by(strain, behavior) %>%
  mutate(total = sum(usv.count)) %>%
  ungroup() %>%
  mutate(percent = usv.count/total*100) %>%
  filter(percent>2)
label.to.keep.behs <- as.vector(unique(beh.total.pers$label))

beh.label.expand <- expand.grid(file.name = unique(mbt_counts$file.name),
                                behavior = c("Retrieval","Hover","LKP"),
                                label = label.to.keep.behs)
usv.beh.counts.filenames <- left_join(beh.label.expand, file.name.key)
usv.beh.counts.join <- left_join(usv.beh.counts.filenames, beh.total.usv.counts)
usv.beh.counts <- usv.beh.counts.join %>% 
  mutate(usv.count = ifelse(is.na(usv.count),0,usv.count))

usv.beh.pers <- usv.beh.counts %>%
  group_by(strain, behavior, rat.id, file.name) %>%
  mutate(total = sum(usv.count)) %>%
  ungroup() %>%
  mutate(percent = usv.count/total)

##now adding chisquare test!!
chisqr.data <- usv.beh.counts %>% group_by(behavior, label, strain) %>% 
  summarise(usv.sum = sum(usv.count)) %>%
  group_by(label, strain) %>%
  mutate(pers=usv.sum/sum(usv.sum))
#sd flats
sdflats <- chisqr.data %>% filter(strain=="SD" & label=="flat")
chisq.test(sdflats$pers)

wkflats <- chisqr.data %>% filter(strain=="WK" & label=="flat") %>% 
  ungroup() %>% mutate(pers = usv.sum/sum(usv.sum))

##usv frequencies and durations by behavior
#see if duration changes by behavior
dur.merged <- merge(data_durs, beh_times_data, by = c("file.name","rat.id","strain"))
beh.dur.data <- dur.merged %>% ungroup() %>%
  filter(start.time > start & start.time < end) %>%
  group_by(file.name, rat.id, strain, behavior, label) %>%
  summarise(m.dur = mean(duration)) %>%
  filter((behavior == "Retrieval" |behavior == "Hover" |behavior == "LKP") & 
           (label %in% label.to.keep.behs))


#see if frequency changes by behavior
freq.merged <- merge(data_freqs, beh_times_data, by = c("file.name","rat.id","strain"))
beh_freq_data <- freq.merged %>% ungroup() %>%
  filter(start.time > start & start.time < end) %>%
  group_by(file.name, rat.id, strain, behavior, label) %>%
  summarise(mean.freq = mean(m.freq)) %>%
  filter((behavior == "Retrieval" |behavior == "Hover" |behavior == "LKP") & 
           (label %in% label.to.keep.behs))


##compare first 15 min of mbt to 15 min mom anesth recording
mbt_15 <- start.rows %>% filter(recording=="MBT" & start.time<900) %>%
  group_by(strain, rat.id, label) %>%
  summarise(mbt.count = length(label))

ma_all <- momanesth_counts %>% ungroup() %>%
  select(label = categories.allowed,
         ma.count = total.counts,
         strain, rat.id)

mbt.ma.join <- left_join(mbt_15, ma_all) %>%
  mutate(ma.count = ifelse(is.na(ma.count), 0, ma.count))

mbt.ma.data <- melt(mbt.ma.join, id.vars = c("strain","rat.id","label"),
                    variable.name = "recording", value.name = "usv.count")

mbt.ma.pers <- mbt.ma.data %>% group_by(strain, label, recording) %>%
  summarise(usv.strain = sum(usv.count)) %>%
  group_by(strain, recording) %>%
  mutate(total = sum(usv.strain)) %>%
  ungroup() %>%
  mutate(percent = usv.strain/total*100) %>%
  filter(percent>2)
label.to.keep.mbtma <- as.vector(unique(mbt.ma.pers$label))

mbt.ma.percents <- mbt.ma.data %>%
  filter(label %in% label.to.keep.mbtma) %>%
  group_by(strain, rat.id, recording) %>%
  mutate(total = sum(usv.count)) %>%
  ungroup() %>%
  mutate(percent = usv.count/total)








