library(cowplot)

#usvs figure
plot_grid(plot1, plot2, plot3, plot4, labels = "AUTO", ncol = 2, align = 'w')

#maternal behaviors figure
plot_grid(plot5, plot6, labels = "AUTO")

#bar graph
plot1 <- mbt_counts %>% 
  filter(categories.allowed == "flat" | categories.allowed == "flat-z" | categories.allowed == "short") %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity", position = position_dodge(.8)) +
  theme_classic() +
  theme(legend.justification = c(0,0), legend.position = c(.5, .5), 
        legend.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size = 14), legend.title = element_blank(),
        axis.text.x = element_text(size = 14)) +
  scale_fill_grey(labels = c("SD","WKY")) +
  xlab("USV Category") +
  ylab("Ultrasonic Vocalizations (n)")

#pie chart
plot2 <- ggplot(mbt.pie.data, aes(x="", y=per, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values=as.vector(colour.key$label.colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar("y", start=0) + 
  facet_wrap(~strain,labeller = as_labeller(c("SD"="SD","WK"="WKY"))) + 
  theme_void() +
  theme(legend.position = "right", aspect.ratio = 1, legend.title = element_blank(),
        legend.key.size = unit(.1,"cm"), legend.background = element_rect(colour = "transparent", fill = "transparent"))

mbt.pie.data <- mbt_counts %>% group_by(strain, categories.allowed) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>%
  mutate(per = total.count/sum(total.count) *100,
         pie.lab = ifelse(per > 2.5, paste0(categories.allowed, " ", round(per), "%"), NA))
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
        panel.grid.minor = element_blank())

mbt.freq.hist <- data_freqs %>% filter(recording == "MBT")
freq.hist <- ggplot(mbt.freq.hist, aes(x = m.freq/1000, fill = strain)) + 
  geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100), position = "identity", colour = "black", alpha=0.5) +
  ylab("(%)")+
  scale_fill_grey() +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_blank()) +
  coord_flip()

mbt.dur.hist <- data_durs %>% filter(recording == "MBT" & duration < 0.4)
dur.hist <- ggplot(mbt.dur.hist, aes(x = duration * 1000, fill = strain)) + 
  geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100), position = "identity", colour = "black", alpha=0.5) +
  ylab("(%)")+
  scale_fill_grey() +
  theme_classic()+
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank())


try.inset <- base.plot.mbt + annotation_custom(grob = inset.plot.mbt, xmin = 165, xmax = 410, ymin = 55, ymax = 95)

plot3.data.mbt <- data_freqs %>% filter((recording == "MBT") & 
                                          (label == "flat" | label == "flat-z" | 
                                             label == "short") &
                                          duration < 0.4)


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

label_colours <- read.csv("data/usv_label_colours.csv", stringsAsFactors = F)
scatter.colour.key <- label_colours %>% filter(label == "flat"|
                                              label == "flat-z"|
                                              label == "short")


##call categories by time
data.set <- start.rows %>% filter(recording == "MBT" & (label == "flat" | label == "flat-z"|
                                                         label == "short"))
bins <- cut(data.set$start.time,6,include.lowest=T, labels = c("1","2","3","4","5","6"))

line.data <- data.set %>%  
  mutate(bin = bins) %>% 
  group_by(bin,strain,file.name, label) %>% 
  summarise(count = length(label))%>% 
  group_by(bin,strain, label) %>% 
  summarise(m.count = mean(count),tot.count = length(count), sem = ifelse(tot.count > 1, sd(count)/sqrt(length(count)),0))
line.data <- rbind.data.frame(line.data, data.frame(bin = c("1","6"),
                                                    strain = c("WK","WK"),
                                                    label = c("flat-z","flat-z"),
                                                    m.count = c(0,0),
                                                    tot.count = c(0,0),
                                                    sem = c(0,0)))

plot4 <- ggplot(line.data, aes(x=bin, y=m.count, group = label, colour=label)) + 
  geom_errorbar(aes(ymin=m.count-sem, ymax=m.count+sem), width=.1,size=.5) +
  geom_line(size=1) +
  geom_point(size=3) + 
  ylab("USVs (n)") +
  xlab("Time in test (5 min intervals)") +
  theme_classic() +
  theme(legend.position = c(.75, .75)) +
  scale_colour_manual(values=scatter.colour.key$colour) +
  facet_wrap(~strain,labeller = as_labeller(c("SD"="SD","WK"="WKY")))


#boxplot freqs
mbt_freqs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.freq/1000, colour = strain)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq/1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Frequency (kHz)") +
  scale_colour_grey()


#boxplot durs
mbt_durs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.dur * 1000, colour = strain)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.dur * 1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Duration (ms)") +
  scale_colour_grey()



#behavioral counts
plot5 <- ggplot(mbt.beh.counts.data, aes(x = behavior, y = score, colour = strain)) +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(size = 2, aes(y = score, group = strain), position = position_jitterdodge(0.2)) +
  scale_colour_grey() +
  theme_classic() +
  theme(legend.position = "none")



#behavioral latencies

ggplot(mbt.beh.lats.data, aes(x = behavior, y = score, colour = strain)) +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(size = 2, aes(y = score, group = strain), position = position_jitterdodge(0.2))

#behavioral durations

ggplot(mbt.beh.durs.data, aes(x = behavior, y = score, colour = strain)) +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(size = 2, aes(y = score, group = strain), position = position_jitterdodge(0.2))

#usv emissions during behaviors
ggplot(beh.total.usv.counts, aes(x = behavior, y = tot.usv, colour = strain)) +
  geom_boxplot(fill = "white", position = position_dodge(0.8)) +
  geom_point(size = 2, aes(y = tot.usv, group = strain), position = position_jitterdodge(0.2))

plot6.data <- beh.usv.counts.data %>% filter(label == "flat" |
                                               label == "flat-z" |
                                               label == "short")
plot6 <- ggplot(plot6.data, aes(x = label, y = usv.count, colour = strain)) +
  geom_boxplot(fill = "white", position = position_dodge(0.8)) +
  geom_point(size = 2, aes(y = usv.count, group = strain), position = position_jitterdodge(0.2)) +
  scale_colour_grey() +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14)) +
  facet_wrap(~behavior, labeller = as_labeller(c("Retrieval"="Grouping","Hover"="Active Caregiving",
                                                 "LKP" = "Nursing"))) +
  xlab("USV Category") +
  ylab("Ultrasonic Vocalizations (n)")

plot6.line.data <- beh.usv.counts.data %>% filter(label == "flat" |label == "flat-z" |
                                                   label == "short") %>%
  group_by(strain, label, behavior) %>%
  summarise(m.count = mean(usv.count), sem = sd(usv.count)/sqrt(length(usv.count))) %>%
  mutate(behavior.key = ifelse(behavior == "Retrieval", "1",
                               ifelse(behavior == "Hover", "2", "3")))


plot6 <- ggplot(plot6.line.data, aes(x = behavior.key, y = m.count, colour = label, group = label)) +
  geom_errorbar(aes(ymin=m.count-sem, ymax=m.count+sem), width=.1,size=.5) +
  geom_line(size=1) +
  geom_point(size=3) + 
  ylab("USVs (n)") +
  xlab("Behavior") +
  theme_classic() +
  theme(legend.position = c(.75, .75)) +
  scale_colour_manual(values=scatter.colour.key$colour) +
  facet_wrap(~strain,labeller = as_labeller(c("SD"="SD","WK"="WKY"))) +
  scale_x_discrete(labels = c("1"="Grouping","2"="Active Caregiving",
                              "3" = "Nursing"))

  



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



