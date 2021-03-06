library(dplyr)
library(reshape2)
source("C:/Users/ituncali/Documents/Master's Thesis/USV-Data-Management/src/count_total.R")
library(stringr)

##counts
data_counts <- read.csv("data/Exp1_All_Labels_Starttimes.csv",stringsAsFactors = F)
allowed.categories <- c("flat", "flat-z", "flat-mz", "short", "short-su", "short-sd",
                        "short-ur", "short-dr", "short-c", "complex", "upward ramp",
                        "downward ramp", "step up", "step down", "multi-step", "multi-step-s",
                        "trill", "trill-c", "trill-f", "inverted-u", "unclear") 
count_list <- by(data = data_counts, INDICES = data_counts$file.name, FUN = function(x) count_total(x, allowed.categories))
count_frame <- do.call(rbind, count_list)
count_frame$file.name <- row.names(count_frame)
row.names(count_frame) <- NULL
count_frame$file.name <- str_replace(count_frame$file.name, pattern = ".[0-9]+$", replacement =  "" )
count_frame <- count_frame %>% group_by(file.name) %>%
  mutate(total.filecounts=sum(total.counts), rel.filecount=total.counts/total.filecounts)
#need to add missing data
wk95.f <- data.frame(categories.allowed = allowed.categories, 
                     total.counts = c(rep(0,21)),
                     file.name = c(rep("T0000305",21)),
                     total.filecounts = c(rep(0,21)),
                     rel.filecount = c(rep(0,21)))
ma.missing.data <- data.frame(categories.allowed = c(rep(allowed.categories,10)),
                           total.counts = c(rep(0,210)),
                           file.name = c(rep("T0000094",21),rep("T0000139",21),rep("T0000122",21),
                                         rep("T0000146",21),rep("T0000147",21),rep("T0000183",21),
                                         rep("T0000188",21),rep("T0000193",21),rep("T0000227",21),
                                         rep("T0000301",21)),
                           total.filecounts = c(rep(0,210)),
                           rel.filecount = c(rep(0,210)))
count_frame_added <- rbind.data.frame(count_frame, wk95.f, ma.missing.data)
##file.name.key
file.name.key.unmelted <- read.csv("data/Exp1_file.name_key.csv", stringsAsFactors = F)
file.name.key <- melt(data=file.name.key.unmelted,id = c("strain","rat.id"),
                      variable.name="recording",value.name = "file.name")

count_frame_2 <- left_join(count_frame_added, file.name.key)

exp1.counts.lme.data <- count_frame_2 %>% group_by(categories.allowed, rat.id) %>% 
  summarise(tot.count =sum(total.counts))
exp1.counts.lme <- lme(tot.count ~ categories.allowed, 
                       random = ~1|rat.id, 
                       data = exp1.counts.lme.data)
anova.lme(exp1.counts.lme)
exp1.counts.sum <- summary(lsmeans(exp1.counts.lme, pairwise ~ categories.allowed, adjust = "Tukey"))[["contrasts"]]

exp1.pers.lme.data <- count_frame_2 %>% group_by(categories.allowed, rat.id) %>% 
  summarise(tot.count =sum(total.counts)) %>% ungroup() %>% mutate(percent = tot.count/sum(tot.count)*100)
exp1.pers.lme <- lme(percent ~ categories.allowed, 
                       random = ~1|rat.id, 
                       data = exp1.pers.lme.data)
anova.lme(exp1.pers.lme)
exp1.pers.sum <- summary(lsmeans(exp1.pers.lme, pairwise ~ categories.allowed, adjust = "Tukey"))[["contrasts"]]

#I'm getting rid of a bunch of calls that made up <3% of the data all together
exp1.counts <- count_frame_2 %>% group_by(categories.allowed) %>% summarise(tot.count =sum(total.counts))
exp1.counts.pers <- exp1.counts %>% ungroup() %>% mutate(percent = tot.count/sum(tot.count) *100)

sum(exp1.counts.pers[order(exp1.counts.pers$percent, decreasing=T)[14:21],]$percent)
calls.to.get.rid.of <- as.vector(exp1.counts.pers[order(exp1.counts.pers$percent, decreasing = T)[14:21],]$categories.allowed)

count_frame_3 <- count_frame_2 %>% group_by(categories.allowed) %>%
  filter(!(categories.allowed == calls.to.get.rid.of[1]|categories.allowed == calls.to.get.rid.of[2]|
             categories.allowed == calls.to.get.rid.of[3]|categories.allowed == calls.to.get.rid.of[4]|
             categories.allowed == calls.to.get.rid.of[5]|categories.allowed == calls.to.get.rid.of[6]|
             categories.allowed == calls.to.get.rid.of[7]|categories.allowed == calls.to.get.rid.of[8]))
count_frame_3$categories.allowed <- droplevels(count_frame_3$categories.allowed)

#now see how many of the remaining calls are emitted per animal per recording
counts.props <- count_frame_3 %>% group_by(categories.allowed, strain, rat.id) %>% 
  summarise(number = paste0(sum(ifelse(total.counts==0,0,1)),"/6"), 
            proportion = sum(ifelse(total.counts==0,0,1))/6) 



##durations
data_dur_readin <- read.csv("data/Exp1_Non_Overlap_Durations.csv", stringsAsFactors = F)
data_durs <- data_dur_readin %>% 
  mutate(label = ifelse(duration < 0.012 & label !="short", "short",
                        ifelse(duration > 0.012 & label == "short", "flat",label))) %>%
  filter(duration < 0.8 & duration > 0.002)


##frequencies
data_read_in <- read.csv("data/Exp1_Filtered_All_Vars.csv", stringsAsFactors = F)
m.freq <- apply(data_read_in[,c(1:50)],1,function(y) mean(y))
data_freqs_readin <- mutate(data_read_in, m.freq = m.freq)
data_freqs <- data_freqs_readin %>% 
  mutate(label = ifelse(duration < 0.012 & label !="short", "short",
                        ifelse(duration > 0.012 & label == "short", "flat",label))) %>%
  filter(duration < 0.8 & duration > 0.002)



##start time counts
categories.allowed.search <- paste0(allowed.categories, "(?!-)")
start_data <- cbind(data_counts, 
                    no.counts = mapply(function(q) sum(str_count(q,categories.allowed.search)),data_counts$label))
#rows.to.rep <- start_data %>% filter(no.counts > 1)
#needed to add extra regex at beginning to account for multi-step-s
xtra.list <- mapply(function(q) str_extract_all(q, pattern = "[a-z]+(-[a-z]+)(-[a-z]+)|[a-z]+(-[a-z]+)|[a-z]+"),
                    start_data$label)
xtra.notlist <- unlist(xtra.list)
xtra.notlist.fix <- str_replace_all(xtra.notlist, 
                                    pattern = c("downward" = "downward ramp",
                                                "upward" = "upward ramp",
                                                "up$" = "step up",
                                                "down$" = "step down"#,
                                                #"inverted$" = "inverted-u",
                                                #"^s$" = "multi-step-s"
                                    ))                           
#first remove the rows that do not contain calls
xtra.rows <- unlist(mapply(function(q) str_subset(q, categories.allowed.search), xtra.notlist.fix))
start.rows <- cbind.data.frame(label = xtra.rows,
                                    rat.id = unlist(mapply(function(q,x) c(rep(q,x)),
                                                           q = start_data$rat.id, x = start_data$no.counts)),
                                    strain = unlist(mapply(function(q,x) c(rep(q,x)),
                                                           q = start_data$strain, x = start_data$no.counts)),
                                    start.time = unlist(mapply(function(q,x) c(rep(q,x)),
                                                               q = start_data$start.time, x = start_data$no.counts)),
                                    file.name = unlist(mapply(function(q,x) c(rep(q,x)),
                                                              q = start_data$file.name, x = start_data$no.counts)),
                                    unique.id = unlist(mapply(function(q,x) c(rep(q,x)),
                                                              q = start_data$unique.id, x = start_data$no.counts)),
                                    recording = unlist(mapply(function(q,x) c(rep(q,x)),
                                                              q = start_data$recording, x = start_data$no.counts)))
#now get rid of calls we said aren't important


