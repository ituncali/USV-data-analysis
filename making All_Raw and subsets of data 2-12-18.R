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

to.adapt <- data.frame(QQ, file.name = filtered_50$file.name, 
                           label = filtered_50$label, 
                           strain = filtered_50$strain, 
                           rat.id = filtered_50$rat.id, 
                           recording = filtered_50$recording, 
                           duration=filtered_50$duration, 
                           start.time = filtered_50$start.time)

write.csv(to.adapt, "C:/Users/ituncali/Documents/Master's Thesis/Data/Experiment 1/csv files/Exp1_Filtered_All_Vars.csv",
          row.names = F)

##to get a table with counts from all three files
plot_frame_2 <- plot_frame_1 %>% group_by(strain,categories.allowed) %>% 
  summarise(sum(total.counts))
names(plot_frame_2) <- c("strain","label","total_raw")

call_per <- read.csv(file.choose(),stringsAsFactors = F)
names(call_per) <- c("label","strain","total_filtered","total_nonoverlap","percent_filtered_from_nonoverlap")
data.file.counts <- left_join(plot_frame_2,call_per)

#right now, percent_filtered is the percent out of non-overlapping calls, we want 
#it to be the percent out of the total raw calls

data.file.counts <- mutate(data.file.counts, percent_filtered_from_raw = 100*(total_filtered/total_raw))

write.csv(data.file.counts,"C:/Users/ituncali/Documents/Master's Thesis/Data/Experiment 1/Exp1_Nonoverlap_Filtered_counts.csv",
          row.names = F)

#to get SD and WK columns
dcast(Strain_recovery, label~strain)

#to get counts from different recordings

count_frame_2 <- left_join(count_frame, file.name.key)




