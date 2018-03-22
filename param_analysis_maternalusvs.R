#prepare count data
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
#need to add the files that had no counts
missing.data <- data.frame(categories.allowed = c(rep(allowed.categories,10)), 
                     total.counts = c(rep(0,210)),
                     file.name = c(rep("T0000094",21),rep("T0000139",21),rep("T0000122",21),
                                   rep("T0000146",21),rep("T0000147",21),rep("T0000183",21),
                                   rep("T0000188",21),rep("T0000193",21),rep("T0000227",21),
                                   rep("T0000301",21)),
                     total.filecounts = c(rep(0,210)),
                     rel.filecount = c(rep(0,210)))

library(dplyr)
count_frame <- rbind.data.frame(count_frame, missing.data)
count_frame_2 <- left_join(count_frame,file.name.key)

maternal_counts <- count_frame_2 %>% filter(recording == "MomAlone" | recording == "PupsSep")

#maternal counts lme
lme_maternal <- lme(total.counts ~ strain * categories.allowed * recording, 
                     random = ~1| rat.id,
                     data = maternal_counts)
anova.lme(lme_maternal)
lsmeans(lme_maternal, pairwise ~ categories.allowed * strain * recording, adjust = "Tukey")
maternal_pairwise_summary <- (summary(lsmeans(lme_maternal,pairwise ~ strain * recording * categories.allowed, adjust = "Tukey")[["contrasts"]]))
View(maternal_pairwise_summary[maternal_pairwise_summary$p.value < 0.05,])


#maternal freqs
data_freqs <- read.csv(file.choose(), stringsAsFactors = F)
m.freq <- apply(data_freqs[,c(1:50)],1,function(y) mean(y))
data_freqs <- mutate(data_freqs, m.freq = m.freq)

#need to combine both recordings because cant compare them because not enough animals
#vocalized in both recordings
maternal_freqs <- data_freqs %>% filter(recording == "MomAlone" | recording == "PupsSep") %>% 
  group_by(strain, label, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label == "flat" | label == "short" | label == "trill-c")

#lme for maternal freqs
lme_maternal_freq <- lme(mean.freq ~ strain * label, 
                   random = ~1|rat.id, 
                   data = maternal_freqs)
anova.lme(lme_maternal_freq)
lsmeans(lme_maternal_freq, pairwise ~ label, adjust = "Tukey")
maternal_pairwise_freq <- (summary(lsmeans(lme_maternal_freq,pairwise~label * strain, adjust = "Tukey")[["contrasts"]]))
View(maternal_pairwise_freq[maternal_pairwise_freq$p.value < 0.05,])

#maternal durs
maternal_durs <- data_durs %>% filter(recording == "MomAlone" | recording == "PupsSep") %>% 
  group_by(strain, label, rat.id) %>% 
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label == "flat" | label == "short" | label == "trill-c")

#lme for maternal freqs
lme_maternal_dur <- lme(mean.dur ~ strain * label, 
                         random = ~1|rat.id, 
                         data = maternal_durs)
anova.lme(lme_maternal_dur)
lsmeans(lme_maternal_dur, pairwise ~ label, adjust = "Tukey")
maternal_pairwise_dur <- (summary(lsmeans(lme_maternal_dur,pairwise~label, adjust = "Tukey")[["contrasts"]]))
View(maternal_pairwise_dur[maternal_pairwise_dur$p.value < 0.05,])


###maternal usvs in the half hour recording
maternal.usvs <- to.adapt %>% filter(label=="trill" | label == "trill-c" |
                                       label == "trill-f" | duration > 0.4) #%>% select(mean.freqs)



freqs_call_strain.M <- maternal.usvs %>% group_by(strain,label) %>% 
  summarise(mean = mean(mean.freqs),
            SEM = sd(mean.freqs)/sqrt(length(mean.freqs)), 
            count = length(mean.freqs))

dur_call_strain.M <- maternal.usvs %>% group_by(strain,label) %>% 
  summarise(mean = mean(duration), 
            SEM = sd(duration)/sqrt(length(duration)), 
            count = length(duration))

freqs_strain.M <- to.adapt%>%group_by(strain)%>%summarise(mean(mean.freqs))
dur_strain.M <- to.adapt%>%group_by(strain)%>%summarise(mean(duration))


##22kHz USVs in these recordings!!
neg.data <- read.csv(file.choose(),stringsAsFactors = F)
neg.data <- left_join(neg.data, file.name.key)

neg.total.counts <- neg.data %>% 
  mutate(neg.col = ifelse(peak.freq.mean. < 30000, 1, 0),
         pos.col = ifelse(peak.freq.mean. >= 30000, 1, 0)) %>%
  group_by(strain, rat.id, recording) %>%
  summarise(neg.count = sum(neg.col), pos.count = sum(pos.col), tot.count = length(peak.freq.mean.))

neg.percent <- neg.total.counts %>% group_by(strain,recording) %>%
  summarise(count.neg = sum(neg.count),
            count.tot = sum(tot.count),
            percent.neg = sum(neg.count)/sum(tot.count)*100)
maternal.negs <- neg.data %>% filter(recording == "MomAlone" | recording == "PupsSep") %>%
  group_by(strain, rat.id, recording) %>%
  filter(peak.freq.mean. < 28000) %>%
  summarise(neg.count = length(peak.freq.mean.))

grid.arrange(tableGrob(neg.percent, rows = NULL))






