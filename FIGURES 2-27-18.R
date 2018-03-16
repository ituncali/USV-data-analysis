###count figures
count_data <- read.csv(file.choose(),stringsAsFactors = F)
library(stringr)
source("C:/Users/ituncali/Documents/Master's Thesis/USV-Data-Management/src/count_total.R")
allowed.categories <- c("flat", "flat-z", "flat-mz", "short", "short-su", "short-sd",
"short-ur", "short-dr", "short-c", "complex", "upward ramp",
"downward ramp", "step up", "step down", "multi-step", "multi-step-s",
"trill", "trill-c", "trill-f", "inverted-U", "unclear")
count_list <- by(data = count_data, INDICES = count_data$file.name, FUN = function(x) count_total(x, allowed.categories))
count_frame <- do.call(rbind, count_list)
count_frame$file.name <- row.names(count_frame)
row.names(count_frame) <- NULL
count_frame$file.name <- str_replace(count_frame$file.name, pattern = ".[0-9]+$", replacement =  "" )
count_frame <- count_frame %>% group_by(file.name) %>%
mutate(total.filecounts=sum(total.counts), rel.filecount=total.counts/total.filecounts)
library(dplyr)
file.name.key.unmelted <- read.csv(file.choose(),stringsAsFactors = F)
library(reshape2)
file.name.key <- melt(data=file.name.key.unmelted,id = c("strain","rat.id"),
                      variable.name="recording",value.name = "file.name")
count_frame_2 <- left_join(count_frame, file.name.key)
#pie chart

pup_counts %>% group_by(categories.allowed) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>%
  mutate(per = round(total.count/sum(total.count) *100, digits = 0),
         pie.lab = ifelse(per > 2, paste0(categories.allowed, " ", per, "%"), NA)) %>%
  ggplot(aes(x="", y=total.count, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5)) +
  coord_polar("y", start=0) + 
  theme_void() +
  theme(legend.position = "none")

#total pup counts bar plot
pup_counts %>% group_by(strain, categories.allowed) %>%
  filter(total.counts > 0) %>%
  ggplot(aes(x=strain, y=total.counts, fill=categories.allowed))+
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Total USV counts per Strain") +
  theme(legend.position = "bottom", legend.title = element_blank())

#mean pup counts bar plot
pup_counts %>% group_by(strain, categories.allowed) %>%
  summarise(m.count = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>%
  filter(m.count > 0) %>%
  ggplot(aes(x=strain, y=m.count, fill=categories.allowed))+
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = m.count - sem, ymax = m.count + sem), width = .1, size = .5) +
  theme_minimal() +
  ggtitle("Mean USV counts per Strain")

#scatter plot
plot(x = pup_freqs_2[pup_freqs_2$label=="flat",]$mean.freq, 
     y = pup_freqs_2[pup_freqs_2$label == "flat",]$mean.dur, 
     xlab = "Frequency", ylab = "Duration", main = "Pup Flat Isolation Calls",
     col = ifelse(pup_freqs_2[pup_freqs_2$label=="flat",]$strain=="SD","blue","red"),
     legend(x="topright", legend = c("SD","WKY"), col=c("blue","red"), pch=1))

data_freqs %>% filter((recording == "Mpupiso" | recording == "Fpupiso") & (label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short")) %>%
  with(plot(x = duration * 1000, 
     y = m.freq/1000, 
     xlab = "Duration (ms)", ylab = "Frequency (kHz)",
     col = ifelse(strain == "SD" & (label == "flat" | label == "flat-z" | label == "flat-mz"), "blue", 
                  ifelse(strain == "SD" & (label == "short"), "yellow",
                         ifelse(strain == "WK" & (label == "short"), "red", "purple")))))
legend(x="topright",legend=c("SD Flat-type","WKY Flat-type", "SD Short", "WKY Short"),col=c("blue","purple","yellow", "red"),pch=1)

data_freqs %>% filter((recording == "Mpupiso" | recording == "Fpupiso") & (label == "flat-mz")) %>%
  with(plot(x = duration, 
            y = m.freq, 
            xlab = "Duration", ylab = "Frequency", main = "Pup Flat-mz Isolation Calls",
            col = ifelse(strain == "SD", "red",
                         "blue")))
legend(x="topright",legend=c("SD","WKY"),col=c("red","blue"),pch=1)

momanesth_counts %>% group_by(strain, categories.allowed) %>%
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>%
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .1, size = .5, position = position_dodge(.9)) +
  theme_bw()


##start time figures
#need to add a column to count_data that has a count for each row
categories.allowed.search <- paste0(allowed.categories, "(?!-)")
start_data <- cbind(count_data, 
    no.counts = mapply(function(q) sum(str_count(q,categories.allowed.search)),count_data$label))
#want to divide times into 6 5-minute bins
bins <- cut(start_data[start_data$recording == "MBT",]$start.time,6,include.lowest=T, labels = c("0-5","5-10","10-15","15-20","20-25","25"))
line.data <- start_data %>% filter(recording == "MBT") %>% mutate(bin = bins) %>% 
  group_by(bin,strain,file.name) %>% 
  summarise(count = sum(no.counts))%>% 
  group_by(bin,strain) %>% 
  summarise(mean = mean(count), sem = sd(count)/sqrt(length(count)),count = length(count))


ggplot(line.data, aes(x=bin, y=mean, group = strain, colour=strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5) +
  geom_line(size=1) +
  geom_point(size=3) + 
  ylab("mean no. calls") +
  xlab("time bin (m)") +
  ggtitle("Number of Calls Emitted Per 5-min Time Bin") +
  theme_bw() +
  theme(legend.justification=c(2,0),
        legend.position=c(1,.5))

base <- pup_starts %>% group_by(strain) %>% 
  summarise(mean.start = mean(lat.to.start),sem.start = sd(lat.to.start)/sqrt(length(lat.to.start)),
            mean.stop = mean(lat.to.stop), sem.stop = sd(lat.to.stop)/sqrt(length(lat.to.stop)))
melt1 <-  melt(data = base[,c("strain","mean.start","mean.stop")], id = "strain", variable.name = "start.or.stop", value.name = "m.time")
melt2 <- melt(data = base[,c("strain","sem.start","sem.stop")],id = "strain", value.name = "sem")  
melt1$sem <- melt2$sem

melt1 %>%
  ggplot(aes(x=start.or.stop, y = m.time, group = strain, fill = strain)) +
  geom_errorbar(aes(ymin = m.time - sem, ymax = m.time + sem), width=.1, size=.5, position=position_dodge(0.9)) +
  geom_bar(stat = "identity", position = position_dodge())

pup_starts %>% group_by(strain) %>% filter(rate != Inf) %>%
  summarise(mean.rate = mean(rate), sem.rate = sd(rate)/sqrt(length(rate))) %>%
  ggplot(aes(x=strain, y = mean.rate)) +
  geom_errorbar(aes(ymin = mean.rate - sem.rate, ymax = mean.rate + sem.rate), width=.1, size=.5, position=position_dodge(0.9)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw() +
  ggtitle("Rate of Isolation USV Emission per Second by Strain") +
  ylab("Total USVs/second")

#hist with error bars for pup isolation recordings
categories.allowed.search <- paste0(allowed.categories, "(?!-)")
start_data <- cbind(data_counts, 
                    no.counts = mapply(function(q) sum(str_count(q,categories.allowed.search)),data_counts$label))
#want to divide times into 5 1-minute bins
bins <- cut(start_data[start_data$recording == "Mpupiso" | start_data$recording == "Fpupiso",]$start.time,5,include.lowest=T, labels = c("1","2","3","4","5"))
line.data <- start_data %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  mutate(bin = bins) %>% 
  group_by(bin,strain,file.name) %>% 
  summarise(count = sum(no.counts))%>% 
  group_by(bin,strain) %>% 
  summarise(mean = mean(count), sem = sd(count)/sqrt(length(count)),count = length(count))

ggplot(line.data, aes(x=bin, y=mean, group = strain, colour=strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5) +
  geom_line(size=1) +
  geom_point(size=3) + 
  ylab("Ultrasonic Vocalizations (n)") +
  xlab("Test Duration (min)") +
  theme_bw() +
  theme(legend.justification=c(2,0),
        legend.position=c(1,.6))

###frequency figures
freq_data <- read.csv(file.choose(),stringsAsFactors = F)

m.freq <- apply(freq_data[,c(1:50)],1,function(y) mean(y))

freq_data <- mutate(freq_data, m.freq = m.freq)

sd.freq <- apply(freq_data[,c(1:50)],1,function(y) sd(y))

freq_data <- mutate(freq_data, sd.freq = sd.freq)

bw.freq <- apply(freq_data[,c(1:50)],1,function(y) max(y)-min(y))

freq_data <- mutate(freq_data, bw.freq = bw.freq)

freq_data$calltype <- ifelse(freq_data$label=="flat" | freq_data$label=="flat-z" | freq_data$label=="flat-mz","flat-type",
                             ifelse(freq_data$label=="short" | freq_data$label=="short-c" | freq_data$label=="short-ur" | freq_data$label=="short-dr" | freq_data$label=="short-su" | freq_data$label=="short-sd", "short-type","FM-type"))

#strain-wise
freq_data %>% filter(recording == "MBT") %>% group_by(strain,calltype) %>%
  summarise(mean = mean(m.freq),sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  ggplot(aes(x=calltype, y=mean/1000, fill=strain)) + 
  geom_errorbar(aes(ymin=(mean/1000)-(sem/1000), ymax=(mean/1000)+(sem/1000)), width=.2,position = position_dodge(.9)) +
  geom_bar(stat = "identity",position=position_dodge()) +
  ylab("mean frequency (kHz)") +
  xlab("USV category") +
  ggtitle("Strain call frequencies") +
  theme_minimal() +
  scale_fill_brewer(palette="Paired")

#pup isolation calls
#plot
#boxplot of pup freqs
pup_freqs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.freq/1000, colour = strain)) +
  theme_bw() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq/1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Frequency (kHz)")

#histogram of pup freqs
data_freqs %>% 
  filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  ggplot(aes(x = m.freq/1000, y=(..count..)/sum(..count..) *100, group = strain, fill = strain)) + geom_histogram(position = "dodge", colour = "black") +
  xlab("Frequency (kHz)") +
  ylab("Percent (%)") +
  geom_vline(xintercept = 75, size = 1, colour = "#FF3721",
                    linetype = "dashed")

#point graph of call freqs by individual
crap <- data_freqs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>%
filter(label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short")
crap$rat.id.fix <- ifelse(crap$recording=="Fpupiso",paste(crap$rat.id,"F"),paste(pup_freqs$rat.id,"M"))

crap %>% group_by(strain, label, rat.id.fix) %>%
ggplot(aes(x = rat.id.fix, y = m.freq, colour = label)) +
theme_bw() +
geom_point(aes(y = m.freq, group = rat.id.fix), size = 1, position = position_jitterdodge(.5)) +
ggtitle("frequencies of 4 major USV types") +
theme(axis.text.x = element_text(angle = 90, hjust = 1))


pup_counts %>% 
  filter(categories.allowed == "flat" | categories.allowed == "flat-z" |categories.allowed == "flat-mz" | categories.allowed == "short") %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity", position = position_dodge(.8)) +
  theme_bw() +
  theme(legend.justification = c(0,1), legend.position = c(.75, .75)) +
  xlab("USV Category") +
  ylab("Ultrasonic Vocalizations (n)")

#boxplot of pup durs
pup_durs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.dur * 1000, colour = strain)) +
  theme_bw() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.dur * 1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Duration (ms)")

#histogram of pup durs
data_durs %>% 
  filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  ggplot(aes(x = duration*1000, y=(..count..)/sum(..count..) *100, group = strain, fill = strain)) + geom_histogram(position = "dodge", colour = "black") +
  xlab("Duration (ms)") +
  ylab("Percent (%)")


#line bar of pup durs
data_durs %>% filter((recording == "Mpupiso" | recording == "Fpupiso") & (label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short")) %>%
  group_by(strain, label) %>%
  summarise(mean = mean(duration), sem = sd(duration)/sqrt(length(duration))) %>%
  ggplot(aes(x = label, y = mean, group = strain, colour = strain)) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1, size=.5) +
  geom_line(size=1) +
  geom_point(size=3) +
  theme_bw() +
  ggtitle("durations of 4 major USV types") +
  theme(legend.justification=c(2,0),
        legend.position=c(.5,.25))

##mom anesth recordings

momanesth_freqs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.freq, colour = strain)) +
  theme_bw() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq, group = strain),position=position_jitterdodge(0.2)) +
  ggtitle("frequencies of 5 major USV types") +
  theme(legend.justification=c(2,0),
        legend.position=c(.9,.15))

momanesth_freqs %>% group_by(strain, label, rat.id) %>%
  ggplot(aes(x = rat.id, y = mean.freq, colour = label)) +
  theme_bw() +
  geom_point(aes(y = mean.freq, group = rat.id), size = 2) +
  geom_errorbar(aes(ymin = mean.freq - sem, ymax = mean.freq + sem), width = .2, size = .1) +
  ggtitle("frequencies of 5 major USV types") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


momanesth_durs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.dur, colour = strain)) +
  theme_bw() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.dur, group = strain),position=position_jitterdodge(0.2)) +
  ggtitle("durations of 5 major USV types") +
  theme(legend.justification=c(2,0),
        legend.position=c(.9,.15))

momanesth_durs %>% group_by(strain, label, rat.id) %>%
  ggplot(aes(x = rat.id, y = mean.dur, colour = label)) +
  theme_bw() +
  geom_point(aes(y = mean.dur, group = rat.id), size = 2) +
  geom_errorbar(aes(ymin = mean.dur - sem, ymax = mean.dur + sem), width = .2, size = .1) +
  ggtitle("durations of 5 major USV types") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##mom calls

View(freq_data %>% filter(recording == "MomAlone" | recording == "PupsSep") %>% 
       group_by(strain, label) %>% 
       summarise(mean = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
                 count = length(m.freq)))

mom_freqs <- freq_data %>% 
  filter(recording == "MomAlone" | recording == "PupsSep") %>% 
  select(label, strain, rat.id, recording, duration, start.time, m.freq, sd.freq, bw.freq)
mom.freq.aov <- aov(m.freq ~ strain*label,data = mom_freqs)
summary(mom.freq.aov)
TukeyHSD(mom.freq.aov)

mom.dur.aov <- aov(duration ~ strain*label, data = mom_durs)
summary(mom.dur.aov)

tukey.mbt.counts.df <- as.data.frame(tukey.mbt.counts$`strain:categories.allowed`)
names <- rownames(tukey.mbt.counts.df)
rownames(tukey.mbt.counts.df) <- NULL
tukey.mbt.counts.df <- cbind(names,tukey.mbt.counts.df)


#mom anesth

freq_data %>% filter(recording == "MA" & strain == "SD") %>% with(hist(m.freq))
freq_data %>% filter(recording == "MA" & strain == "WK") %>% with(hist(m.freq))

freq_data %>% filter(recording == "MBT") %>% with(hist(m.freq))

freq_data %>% filter(recording == "PupsSep") %>% with(hist(m.freq))


freq_data %>% filter(recording == "PupsSep") %>% with(hist(m.freq))

##duration figures
dur_data <- read.csv(file.choose(),stringsAsFactors = F)

freq_data %>% filter(recording == "MBT") %>% group_by(strain,calltype) %>%
  summarise(mean = mean(duration),sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  ggplot(aes(x=calltype, y=mean, fill=strain)) + 
  geom_errorbar(aes(ymin=(mean)-(sem), ymax=(mean)+(sem)), width=.2,position = position_dodge(.9)) +
  geom_bar(stat = "identity",position=position_dodge()) +
  ylab("mean duration (s)") +
  xlab("USV category") +
  ggtitle("Strain call frequencies") +
  theme_minimal() +
  scale_fill_brewer(palette="Paired")

##time figures
bins <- cut(start_data[start_data$recording == "MA",]$start.time,5,include.lowest=T, labels = c("1","2","3","4","5"))
line.data <- start_data %>% filter(recording == "MA") %>% 
  mutate(bin = bins) %>% 
  group_by(bin,strain,file.name) %>% 
  summarise(count = sum(no.counts))%>% 
  group_by(bin,strain) %>% 
  summarise(mean = mean(count), sem = sd(count)/sqrt(length(count)),count = length(count))

ggplot(line.data, aes(x=bin, y=mean, group = strain, colour=strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5) +
  geom_line(size=1) +
  geom_point(size=3) + 
  ylab("mean no. calls/3min") +
  xlab("time bin (3m)") +
  ggtitle("Number of Calls Emitted Per 3-min Time Bin") +
  theme_bw() +
  theme(legend.justification=c(2,0),
        legend.position=c(1,.5))

##mom alone and pups separated
maternal_counts %>% group_by(strain, recording, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = recording, y = mean, group = categories.allowed, colour = categories.allowed)) + 
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .1, size = .5, position = position_dodge(.5)) + 
  geom_line(size=1,position = position_dodge(.5)) + geom_point(size=3, position = position_dodge(.5)) + facet_wrap(~strain) +
  scale_colour_manual(values=c("#BF383E", "#73BF5C", "#1F5694", "#13355C", "#3594FF", "#FF4540", "#FFC19E", "#94705C", "#FF71B7", "#FF1FE0", "#E51CCA", "#BF5BA8", "#FFB3FA", "#5C0B51", "#BFBE3F", "#949339", "#66BF9A", "#447F67", "#224033", "#5C5B1E", "#375C2C"))

maternal_counts %>% group_by(strain, recording, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) + 
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .1, size = .5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity",position = "dodge") +
  scale_fill_grey()



