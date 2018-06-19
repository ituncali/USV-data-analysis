momanesth_counts_all <- count_frame_2 %>% filter(recording == "MA")

ma.per <- momanesth_counts_all %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(usv.count = sum(total.counts)) %>%
  group_by(strain) %>% mutate(percent=usv.count/sum(usv.count)*100) %>%
  filter(percent>2)

label.to.keep.ma <- as.vector(unique(ma.per$categories.allowed))

momanesth_counts <- momanesth_counts_all %>% filter(categories.allowed %in% label.to.keep.ma)


ma.per.2 <- momanesth_counts %>% 
  group_by(strain, rat.id, file.name) %>%
  mutate(usv.per=total.counts/sum(total.counts)) 


##frequencies
momanesth_freqs <- data_freqs %>% filter(recording == "MA") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label %in% label.to.keep.ma)


##durations
momanesth_durs <- data_durs %>% filter(recording == "MA") %>% 
  group_by(strain, label, file.name,recording, rat.id) %>% 
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label %in% label.to.keep.ma)


##compare isolation recordings to first 5 min 
#COUNTS
mapup.compare <- ma_pup_counts %>% group_by(strain, label,recording) %>% 
  summarise(total.usv=sum(total.counts)) %>% 
  group_by(strain, recording) %>% 
  mutate(total.all = sum(total.usv)) %>% ungroup() %>% 
  mutate(percent = total.usv/total.all*100) %>%
  filter(percent >= 2)
label.to.keep.mapup <- as.vector(unique(mapup.compare$label))

ma_pup_counts_2 <- ma_pup_counts %>%
  filter(label %in% label.to.keep.mapup)

ma_pup_per_2 <- ma_pup_counts_2 %>%
  group_by(strain, recording, rat.id) %>%
  mutate(total.filecount = sum(total.counts)) %>%
  ungroup() %>%
  mutate(percent = total.counts/total.filecount)
#FREQUENCIES
ma_freqs_comp <- momanesth_freqs %>% filter(label %in% label.to.keep.mapup)

pup_freqs_comp <- data_freqs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label %in% label.to.keep.mapup)

ma.pup.freqs <- rbind.data.frame(ma_freqs_comp, pup_freqs_comp)
ma.pup.freqs <- ma.pup.freqs %>% ungroup() %>%
  filter(!label=="short-c")
#DURATIONS
ma_durs_comp <- momanesth_durs %>% filter(label %in% label.to.keep.mapup)

pup_durs_comp <- data_durs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label %in% label.to.keep.mapup)

ma.pups.durs <- rbind.data.frame(ma_durs_comp, pup_durs_comp)
ma.pups.durs <- ma.pups.durs %>% ungroup() %>%
  filter(!(label=="short-c"|label=="short-ur"))

##three time points
ma.xtra <- start.rows %>% filter(recording == "MA")

three.times.data <- ma.xtra %>% 
  mutate(time.group = ifelse(start.time < 180, "1",
                             ifelse(start.time >= 420 & start.time < 540, "2",
                                    ifelse(start.time > 720, "3", NA)))) %>%
  filter(!is.na(time.group)) %>%
  group_by(strain, time.group, rat.id, file.name, label) %>%
  summarise(time.count = length(label))

three.per <- three.times.data %>% group_by(strain, time.group, label) %>% 
  summarise(strain.total=sum(time.count)) %>% 
  group_by(strain, time.group) %>% 
  mutate(total=sum(strain.total)) %>% ungroup() %>% 
  mutate(percent = strain.total/total*100) %>%
  filter(percent>2)

label.to.keep.bme <- as.vector(unique(three.per$label))

to.join <- expand.grid(label=label.to.keep.bme, 
                       time.group=c("1","2","3"), 
                       file.name=unique(momanesth_counts$file.name))
to.join <- left_join(to.join, file.name.key)

pru <- to.join %>% left_join(three.times.data) %>%
  mutate(time.count = ifelse(is.na(time.count), 0, time.count)) %>%
  droplevels()

three.times.profile <- pru %>% group_by(strain, rat.id, time.group) %>% 
  mutate(tot.count.in.time = sum(time.count)) %>% 
  group_by(strain, rat.id, time.group, label) %>% 
  mutate(percent = time.count/tot.count.in.time)


##before vs. during grouping
group.induced.data <- start.rows %>% filter(recording == "MA") %>% 
  mutate(time.group = ifelse(start.time < 600 & start.time > 570, "b4",
                             ifelse(start.time >= 600 & start.time < 630, "during", NA))) %>%
  filter(!is.na(time.group)) %>%
  group_by(strain, time.group, rat.id, file.name, label) %>%
  summarise(time.count = length(label))

group.per <- group.induced.data %>% group_by(strain, time.group, label) %>% 
  summarise(strain.total=sum(time.count)) %>% 
  group_by(strain, time.group) %>% 
  mutate(total=sum(strain.total)) %>% ungroup() %>% 
  mutate(percent = strain.total/total*100) %>%
  filter(percent>2)
label.to.keep.group <- as.vector(unique(group.per$label))

to.join <- expand.grid(label=label.to.keep.group, 
                       time.group=c("b4","during"), 
                       file.name=unique(momanesth_counts$file.name))
to.join <- left_join(to.join, file.name.key)

all.levels <- to.join %>% left_join(group.induced.data) %>%
  mutate(time.count = ifelse(is.na(time.count), 0, time.count)) %>%
  droplevels()

all.levs.profile <- all.levels %>% group_by(strain, rat.id, time.group) %>% 
  mutate(tot.count.in.time = sum(time.count)) %>% 
  group_by(strain, rat.id, time.group, label) %>% 
  mutate(rel.call.count = ifelse(tot.count.in.time >0, time.count/tot.count.in.time, 0))



