##counts
library(dplyr)
library(reshape2)
source("C:/Users/ituncali/Documents/Master's Thesis/USV-Data-Management/src/count_total.R")
library(stringr)

data_counts <- read.csv("data/Exp2_All_Labels_Starttimes.csv",stringsAsFactors = F)
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

file.name.key <- read.csv("data/Exp2_file_name_key.csv", stringsAsFactors = F)

count_frame_2 <- left_join(count_frame, file.name.key)

count_frame_3 <- count_frame_2 %>% filter(!(file.name=="T0000399"|
                                              file.name=="T0000410"|
                                              file.name=="T0000291"))

##durations
data_dur_readin <- read.csv("data/Exp2_Non_Overlap_Durations.csv", stringsAsFactors = F)
data_durs_2 <- data_dur_readin %>% 
  mutate(label = ifelse(duration < 0.012 & label !="short", "short",
                        ifelse(duration > 0.012 & label == "short", "flat",label))) %>%
  filter(duration < 0.8 & duration > 0.002)
data_durs <- data_durs_2 %>% filter(!(file.name=="T0000399"|
                                       file.name=="T0000410"|
                                       file.name=="T0000291"))


##frequencies
data_read_in <- read.csv("data/Exp2_Filtered_All_Vars.csv", stringsAsFactors = F)
m.freq <- apply(data_read_in[,c(1:50)],1,function(y) mean(y))
data_freqs_readin <- mutate(data_read_in, m.freq = m.freq)
data_freqs_2 <- data_freqs_readin %>% 
  mutate(label = ifelse(duration < 0.012 & label !="short", "short",
                        ifelse(duration > 0.012 & label == "short", "flat",label))) %>%
  filter(duration < 0.8 & duration > 0.002)
data_freqs <- data_freqs_2 %>% filter(!(file.name=="T0000399"|
                                          file.name=="T0000410"|
                                          file.name=="T0000291"))


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
start.rows.almost <- cbind.data.frame(label = xtra.rows,
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
start.rows <- start.rows.almost %>% filter(!(file.name=="T0000399"|
                                               file.name=="T0000410"|
                                               file.name=="T0000291"))
