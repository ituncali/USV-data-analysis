#start with Raw_CleanUp.r
#Choose the folder in Master's Thesis/USV-Data-Management/data/exp1
#this will have every raw excel file
#run until the creation of which_are_wk

##instead of adding additional columns based on vectors of file names, just read in
##a .csv file containing all file names and their respective strain, rat.id and recording

library(reshape2)
file.name.key.unmelted <- read.csv(file.choose(),stringsAsFactors = F)

file.name.key <- melt(data=file.name.key.unmelted,id = c("strain","rat.id"),
                      variable.name="recording",value.name = "file.name")

BIG <- left_join(BIG,file.name.key)


write.csv(BIG,"C:/Users/ituncali/Documents/Master's Thesis/Data/Experiment 1/Exp1_All_Raw.csv",
          row.names = F)


####Now to make the sub-setted data!!!


#1. all usvs w/ labels & start times
BIG_1 <- select(BIG, unique.id:label,start.time,strain,rat.id,recording)

write.csv(BIG_1,"C:/Users/ituncali/Documents/Master's Thesis/Data/Experiment 1/Exp1_All_Labels_Starttimes.csv",
          row.names = F)

#2. non-overlapping usvs w/ labels & durations

allowed.categories <- c("flat", "flat-z", "flat-mz", "short", "short-su", "short-sd",
                        "short-ur", "short-dr", "short-c", "complex", "upward ramp",
                        "downward ramp", "step up", "step down", "multi-step", "multi-step-s",
                        "trill", "trill-c", "trill-f", "inverted-u", "unclear") 

BIG_non_overlap <- filter(BIG, label %in% allowed.categories)

BIG_2 <- select(BIG_non_overlap, unique.id:duration,strain,rat.id,recording)

write.csv(BIG_2,"C:/Users/ituncali/Documents/Master's Thesis/Data/Experiment 1/Exp1_Non_Overlap_Durations.csv",
          row.names = F)

#3. filtered usvs w/ labels & freqs
#need to go into Interpolation_CleanUp.r for this
#need to start with non-overlapping calls:

write.csv(BIG_non_overlap, "C:/Users/ituncali/Documents/Master's Thesis/USV-Interpolation-CleanUp/data/Exp1_Non_Overlap_all_vars.csv",
          row.names=F)

QQ <- read.csv("data/QQ.csv",stringsAsFactors = F)
filtered_50 <- read.csv("data/filtered_50.csv", stringsAsFactors = F)

to.adapt <- data.frame(QQ, file.name = filtered_50$file.name, 
                           label = filtered_50$label, 
                           strain = filtered_50$strain, 
                           rat.id = filtered_50$rat.id, 
                           recording = filtered_50$recording, 
                           duration=filtered_50$duration, 
                           start.time = filtered_50$start.time)

write.csv(to.adapt, "C:/Users/ituncali/Documents/Master's Thesis/Data/Experiment 1/Exp1_Filtered_All_Vars.csv",
          row.names = F)

##to get a table with counts from all three files
plot_frame_2 <- count_frame_2 %>% group_by(strain,categories.allowed) %>% 
  summarise(sum(total.counts))
names(plot_frame_2) <- c("strain","label","total_raw")

Strain_recovery <- read.csv("C:/Users/ituncali/Documents/Master's Thesis/USV-Interpolation-CleanUp/Strain_call_recovery.csv",stringsAsFactors = F)
names(Strain_recovery) <- c("label","strain","total_filtered","total_nonoverlap","percent_filtered_from_nonoverlap")
data.file.counts <- left_join(plot_frame_2,Strain_recovery)

#right now, percent_filtered is the percent out of non-overlapping calls, we want 
#it to be the percent out of the total raw calls

data.file.counts <- mutate(data.file.counts, percent_filtered_from_raw = 100*(total_filtered/total_raw))

write.csv(data.file.counts,"C:/Users/ituncali/Documents/Master's Thesis/Data/Experiment 1/Exp1_Nonoverlap_Filtered_counts.csv",
          row.names = F)

#to get SD and WK columns
dcast(Strain_recovery, label~strain)

grid.arrange(tableGrob(dcast(Strain_recovery, label~strain), rows = NULL, cols = c("USV Category", "SD (%)","WKY (%)"), 
                       theme = ttheme_default(base_size = 6)))

#to get counts from different recordings

missing.data <- data.frame(categories.allowed = c(rep(allowed.categories,10)),
total.counts = c(rep(0,210)),
file.name = c(rep("T0000094",21),rep("T0000139",21),rep("T0000122",21),
rep("T0000146",21),rep("T0000147",21),rep("T0000183",21),
rep("T0000188",21),rep("T0000193",21),rep("T0000227",21),
rep("T0000301",21)),
total.filecounts = c(rep(0,210)),
rel.filecount = c(rep(0,210)))

count_frame <- rbind.data.frame(count_frame, missing.data)

wk95.f <- data.frame(categories.allowed = allowed.categories,
total.counts = c(rep(0,21)),
file.name = c(rep("T0000305",21)),
total.filecounts = c(rep(0,21)),
rel.filecount = c(rep(0,21)))

count_frame <- rbind.data.frame(count_frame, wk95.f)
count_frame_2 <- left_join(count_frame,file.name.key)
count_frame_2 <- count_frame_2 %>% select(label = categories.allowed, file.name, total.counts)

freq_frame <- data_freqs %>% group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(freq.count = length(m.freq)) %>%
  ungroup() %>%
  select(file.name, label, freq.count)

dur_frame <- data_durs %>% group_by(strain, label, recording, file.name, rat.id) %>%
  summarise(dur.count = length(duration)) %>%
  ungroup() %>%
  select(file.name, label, dur.count)

join1 <- left_join(count_frame_2, dur_frame)
join2 <- left_join(join1, freq_frame)
counts_by_file <- left_join(join2, file.name.key)
counts_by_file <- counts_by_file %>% mutate(dur.count = ifelse(is.na(dur.count), 0, dur.count),
                                            freq.count = ifelse(is.na(freq.count),0, freq.count))

per_by_file <- counts_by_file %>% mutate(percent = ifelse(dur.count > 0, freq.count/dur.count, 0)) %>%
  select(label, file.name, strain, rat.id, recording, percent)

per.lme <- lme(percent ~ label * strain * recording, random = ~1|rat.id, data = per_by_file)
anova.lme(per.lme)
per.sum <- summary(lsmeans(per.lme, pairwise ~ strain|label, adjust = "Tukey"))[["contrasts"]]
View(per.sum[per.sum$p.value < 0.05,])

per.lab.sum <- summary(lsmeans(per.lme, pairwise ~ strain*label|recording, adjust = "Tukey"))[["contrasts"]]
View(per.lab.sum[per.lab.sum$p.value < 0.05,])

#without recording in account
per.strain.label <- counts_by_file %>% group_by(label, strain, rat.id) %>% 
  summarise(sum.dur = sum(dur.count), sum.freq = sum(freq.count)) %>% 
  mutate(percent = ifelse(sum.dur >0, sum.freq/sum.dur, 0))
per.sl.lme <- lme(percent ~ strain * label, random = ~1|rat.id, data = pru)
anova.lme(per.sl.lme)


per.sl.sum <- summary(lsmeans(per.sl.lme, pairwise ~ strain|label, adjust = "Tukey"))[["contrasts"]]
View(per.sl.sum[per.sl.sum$p.value < 0.05,])




