library(cowplot)

#usvs figure
plot_grid(plot1, plot2, plot3, labels = "AUTO", ncol = 2, align = 'w')

#maternal behaviors figure
plot_grid(plot5, empty, labels = "AUTO", ncol=1)

#bar graph
plot1 <- mbt_counts %>% 
  #filter(categories.allowed == "flat" | categories.allowed == "flat-z" | 
  #         categories.allowed=="flat-mz" |categories.allowed == "short"|
  #         categories.allowed=="trill") %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity", position = position_dodge(.8), colour="black",
           width=.7) +
  theme_classic() +
  theme(legend.justification = c(0,0), legend.position = c(.5, .5), 
        legend.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size = 17), legend.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, size = 17), 
        axis.text.y = element_text(size = 17),
        axis.title=element_text(size=14)) +
  scale_fill_grey(labels = c("SD (n=8)","WKY (n=7)")) +
  xlab("USV Category") +
  ylab("USVs / 30min (n)")

#pie chart
plot2 <- ggplot(mbt.pie.data, aes(x="", y=per, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values=as.vector(colour.key$label.colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 6) +
  coord_polar("y", start=0) + 
  facet_wrap(~strain,labeller = as_labeller(c("SD"="SD","WK"="WKY"))) + 
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        aspect.ratio = 1, 
        legend.title = element_blank(),
        legend.key.size = unit(.1,"cm"),
        legend.background = element_rect(colour = "transparent", fill = "transparent"),
        strip.text=element_text(size=20))

mbt.pie.data <- count_frame_2 %>% filter(recording=="MBT") %>% 
  group_by(strain, categories.allowed) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>% group_by(strain) %>%
  mutate(per = total.count/sum(total.count) *100,
         pie.lab = ifelse(per > 2.5, paste0(categories.allowed, " ", round(per), "%"), NA))
label_colours <- read.csv("data/usv_label_colours.csv", stringsAsFactors = F)
names(label_colours) <- c("categories.allowed", "colour")
colour.key <- left_join(mbt.pie.data, label_colours)
colour.key <- colour.key %>% group_by(categories.allowed) %>%
  summarise(label.colour = unique(colour))


#scatter plot freqs and durs
library(gridExtra)
plot3 <- grid.arrange(dur.hist, empty, try.inset, freq.hist, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.line = element_blank())

base.plot.mbt <- ggplot(plot3.data.mbt, aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 2, aes(colour = strain), alpha = 0.3) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_classic() +
  scale_colour_grey() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=14))

mbt.freq.hist <- data_freqs %>% filter(recording == "MBT")
freq.hist <- ggplot(mbt.freq.hist, aes(x = m.freq/1000, fill = strain)) + 
  geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100), position = "identity", colour = "black", alpha=0.5) +
  ylab("(%)")+
  scale_fill_grey() +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x=element_text(size=20),
        axis.title.x=element_text(size=14)) +
  coord_flip()

mbt.dur.hist <- data_durs %>% filter(recording == "MBT" & duration < 0.4)
dur.hist <- ggplot(mbt.dur.hist, aes(x = duration * 1000, fill = strain)) + 
  geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100), 
                 position = "identity", colour = "black", alpha=0.5) +
  ylab("(%)")+
  scale_fill_grey() +
  theme_classic()+
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=14))


try.inset <- base.plot.mbt + annotation_custom(grob = inset.plot.mbt, xmin = 165, xmax = 610, ymin = 52, ymax = 95)

plot3.data.mbt <- data_freqs %>% filter((recording == "MBT") & 
                                          (label %in% label.to.keep))
label_colours <- read.csv("data/usv_label_colours.csv", stringsAsFactors = F)
scatter.colour.key <- label_colours %>% filter(label %in% label.to.keep)

inset.plot.mbt <- ggplotGrob(ggplot(plot3.data.mbt, aes(x = duration * 1000, y = m.freq/1000)) +
                           geom_point(size = 2, aes(colour = label), alpha = 0.5) +
                           scale_colour_manual(values=scatter.colour.key$colour) +
                           xlab("Duration (ms)") +
                           ylab("Frequeny (kHz)") +
                           theme_classic() +
                             theme(legend.background = element_rect(colour = "transparent", fill = "transparent"), 
                                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                   legend.position = c(.7,.85), legend.title = element_blank(), 
                                   legend.direction = "horizontal",legend.text = element_text(size = 6),
                                   plot.background = element_rect(colour = "transparent", fill = "transparent"),
                                   axis.title = element_text(size = 6),
                                   axis.text = element_text(size = 5)))




##duration bars
plot4.data <- mbt_durs %>% group_by(strain, label) %>%
  summarise(m.dur = mean(mean.dur)*1000, sem = sd(mean.dur)/sqrt(length(mean.dur))*1000)

plot4 <- ggplot(plot4.data, aes(x=label, y=m.dur, fill=strain)) +
  geom_bar(stat="identity", position=position_dodge(.8)) +
  geom_errorbar(aes(ymin=m.dur-sem, ymax=m.dur+sem),size=.5, width=.1,
                position=position_dodge(.8)) +
  xlab("USV Category") +
  ylab("Duration (ms)") +
  theme_classic() +
  theme(legend.position="none",axis.text.x = element_text(size = 14)) +
  scale_fill_grey()
  


#behavioral counts
library(gridExtra)
plot5 <- grid.arrange(behavior.count, empty, ncol=2, widths=c(5, 3), heights=c(4))


behavior.count <- ggplot(mbt.beh.counts.data, aes(x = behavior, y = score, colour = strain)) +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(size = 2, aes(y = score, group = strain), position = position_jitterdodge(0.2)) +
  scale_colour_grey(start=0, end=.6, labels = c("SD","WKY")) +
  theme_classic() +
  theme(legend.position = c(.2,.7), legend.text = element_text(size=17),
        legend.title = element_blank(), axis.text=element_text(size=17),
        axis.title=element_text(size=14)) +
  scale_x_discrete(labels = c(retrieval = "Retrieval",
                              mouthing = "Mouthing",
                              corporal = "Corporal 
Licking",
                              anogenital = "Anogenital 
Licking",
                              nest.building = "Nest Build")) +
  xlab("Behavior") +
  ylab("# Behaviors / 30min")



#behavioral latencies

ggplot(mbt.beh.lats.data, aes(x = behavior, y = score, colour = strain)) +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(size = 2, aes(y = score, group = strain), position = position_jitterdodge(0.2))

#behavioral durations

ggplot(mbt.beh.durs.data, aes(x = behavior, y = score, colour = strain)) +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(size = 2, aes(y = score, group = strain), position = position_jitterdodge(0.2))

#usv emissions during behaviors

plot6 <- ggplot(beh.total.usv.counts, aes(x = behavior, y = usv.count, fill = strain)) +
  geom_bar(stat = "identity",position = position_dodge(0.8)) +
  scale_colour_grey(start = 0, end = .6, labels = c("SD","WKY")) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("Behavior") +
  ylab("Latency (s)") +
  facet_wrap(~label)



##pie chart usv profiles during behaviors
plot7.data <- beh.type.per %>% group_by(strain, behavior, label) %>%
  summarise(tot.count = sum(usv.counts)) %>%
  group_by(strain, behavior) %>%
  mutate(per = tot.count/sum(tot.count)*100, 
         pie.lab = ifelse(per > 2.5, paste0(label, " ", round(per),
                                            "%"), NA),
         behavior_f = factor(behavior, levels = c("Retrieval","Hover","LKP")),
         strain_f = factor(strain, levels = c("SD","WK")))

label_colours <- read.csv("data/usv_label_colours.csv", stringsAsFactors = F)

colour.key <- label_colours %>% filter(label == "flat"| label == "flat-z"|
                                       label == "short")


plot7 <- ggplot(plot7.data, aes(x="", y=per, fill=label))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values= as.vector(colour.key$colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar("y", start=0) + 
  facet_wrap(~strain_f*behavior_f) + 
  theme_void() +
  theme(legend.position = "right", aspect.ratio = 1, legend.title = element_blank(),
        legend.key.size = unit(.1,"cm"), legend.background = element_rect(colour = "transparent", fill = "transparent"))




#behavioral latency histogram
ggplot(mbt.beh.lats.data, aes(x = score, group = strain, fill = strain)) +
  geom_histogram(aes(y = ..count..), position = "dodge",
                 binwidth = 100)


##summary statistics table
mbt_freqs_sw <- data_freqs %>% filter(recording == "MBT") %>% 
  select(strain, duration, m.freq, label)
mbt.freq.table <- mbt_freqs_sw %>% select(-duration) %>% 
  mutate(m.freq = m.freq/1000) %>%
  group_by(strain,label) %>% 
  summarise_all(funs(length, mean, sd, median, min, max))
names(mbt.freq.table) <- c("Strain","USV Category","N","Mean","Standard Dev","Median","Min","Max")

mbt_durs_sw <- data_durs %>% filter(recording == "MBT")
mbt.dur.table <- mbt_durs_sw %>% select(label, duration, strain) %>%
  mutate(duration = duration*1000) %>%
  group_by(strain, label) %>%
  summarise_all(funs(length, mean, sd, median, min, max))
names(mbt.dur.table) <- c("Strain","USV Category","N","Mean","Standard Dev","Median","Min","Max")

mbt.dur.freq.table <- merge(mbt.freq.table, mbt.dur.table, by = c("Strain","USV Category"))

##summary stats of pup parameters
table.durs <- data_durs %>% filter((recording=="MA"|recording=="Mpupiso"|
                                     recording=="Fpupiso") & (label=="flat"|
                                      label=="flat-z"|
                                      label=="short")) %>% 
  group_by(strain, label, rat.id) %>% 
  summarise(Min=round(min(duration)*1000,3), Mean=round(mean(duration)*1000,3),
            SEM = round(sd(duration)/sqrt(length(duration))*1000,3),
            fiveninetyfive=paste0(round(quantile(duration,0.05,type=3)*1000,3),", ",
                                  round(quantile(duration,0.95,type=3)*1000,3)),
            Max=round(max(duration)*1000,3))

table.freqs <- data_freqs %>% filter((recording=="MA"|recording=="Mpupiso"|
                                       recording=="Fpupiso") & (label=="flat"|
                                        label=="flat-z"|
                                          label=="short")) %>%
  group_by(strain, label, rat.id) %>% 
  summarise(Min=round(min(m.freq)/1000,3),Mean=round(mean(m.freq)/1000,3),
            SEM = round(sd(m.freq)/sqrt(length(m.freq))/1000,3),
            fiveninetyfive=paste0(round(quantile(m.freq,0.05,type=3)/1000,3),", ",
                                  round(quantile(m.freq,0.95,type=3)/1000,3)),
            Max=round(max(m.freq)/1000,3))

#summary stats of behavioral latencies and durations
mbt.lats.sum <- mbt_behavior %>% filter(behavior == "lat.retrieve" |
                                               behavior == "lat.group" |
                                               behavior == "lat.hover" |
                                               behavior == "lat.nurse") %>% 
  group_by(strain, behavior) %>% 
  summarise(Min = round(min(score)/60,2), Median = round(median(score)/60,2), 
            SIQrange = paste0(round((median(score) - 
                                ((quantile(score,.75,type=3)-quantile(score,.25,type=3))/2))/60,2),
                              ", ",
            round((median(score) + ((quantile(score,.75,type=3)-quantile(score,.25,type=3))/2))/60,2)),
            Max = round(max(score)/60,2))

mbt.beh.durs.sum <- with.added %>%
  group_by(strain, behavior) %>% 
  summarise(Min = round(min(score)/60,2), Median = round(median(score)/60,2), 
            SIQrange = paste0(round((median(score) - ((quantile(score,.75,type=3)-quantile(score,.25,type=3))/2))/60,2),
                              ", ",
                              round((median(score) + ((quantile(score,.75,type=3)-quantile(score,.25,type=3))/2))/60,2)),
            Max = round(max(score)/60,2))
to.add <- beh_times_data %>% group_by(strain, behavior, file.name,rat.id) %>% 
  summarise(score = sum(beh.dur)) %>%
  filter(behavior=="Retrieval")
with.added <- rbind.data.frame(mbt_behavior[mbt_behavior$behavior == "dur.hover"|
                                              mbt_behavior$behavior == "dur.nurse",],
                               to.add)



###mom usvs
ggplot(mom_counts_1, aes(x = start.time/60, fill = label)) +
  geom_histogram(colour="black", binwidth=1, position="identity") +
  scale_fill_manual(values=c("#000000","#000000","#000000","#FF0000","#FF0000")) +
  theme_classic() +
  facet_wrap(~rat.id)+
  theme(legend.justification = c(0,0), 
        legend.position = c(0.85, 0),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        axis.text=element_text(size=17),
        axis.title=element_text(size=17),
        legend.text=element_text(size=17)) +
  xlab("Start Time (min)") +
  ylab("#USVs / 30min (n)")

##hists totals
ggplot(mom_counts_1, aes(x = start.time, fill = label)) +
  geom_histogram(colour="black", binwidth=60, position="identity") +
  scale_fill_manual(values=c("#000000","#000000","#000000","#FF0000","#FF0000")) +
  theme_classic() +
  facet_wrap(~strain)+
  theme(legend.justification = c(0,0), 
        legend.position = c(0.85, .45),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  xlab("Start Time (s)") +
  ylab("USVs (n)")

##mom counts
bar.moms <- mom_counts %>% group_by(strain, label) %>%
  summarise(m.count = mean(usv.count),sem=sd(usv.count)/sqrt(length(usv.count)))

mom.plot <- ggplot(bar.moms, aes(x = label, y=m.count, fill = strain)) +
  geom_bar(colour="black", stat="identity",position=position_dodge(.6),
           width=.5) +
    geom_errorbar(aes(ymin=m.count-sem,ymax=m.count+sem),
                  position=position_dodge(.6),
                  size=.1, width=.2) +
  scale_fill_grey(labels=c("SD","WKY")) +
  theme_classic() +
  theme(legend.justification = c(0,0), 
        legend.position = c(0.85, .45),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        axis.text=element_text(size=17),
        axis.title=element_text(size=14)) +
  xlab("USV Category") +
  ylab("#USVs / 30min (n)")
