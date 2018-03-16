library(cowplot)


plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, labels = "AUTO", ncol = 2, align = 'w')


#bar graph
plot1 <- mbt_counts %>% 
  filter(categories.allowed == "flat" | categories.allowed == "flat-z" |categories.allowed == "flat-mz" | categories.allowed == "short" | categories.allowed == "short-c" | categories.allowed == "complex" | categories.allowed == "trill" | categories.allowed == "trill-c") %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity", position = position_dodge(.8)) +
  theme_classic() +
  theme(legend.justification = c(0,0), legend.position = c(.5, .5), legend.background = element_rect(colour = "transparent", fill = "transparent")) +
  scale_fill_grey() +
  xlab("USV Category") +
  ylab("Ultrasonic Vocalizations (n)")

#pie chart
plot2 <- ggplot(mbt.pie.data, aes(x="", y=per, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values=c("#BF383E", "#73BF5C", "#1F5694", "#13355C", "#3594FF", "#FF4540", "#FFC19E", "#94705C", "#FF71B7", "#FF1FE0", "#E51CCA", "#BF5BA8", "#FFB3FA", "#5C0B51", "#BFBE3F", "#949339", "#66BF9A", "#447F67", "#224033", "#5C5B1E", "#375C2C")) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 2) +
  coord_polar("y", start=0) + facet_wrap(~strain) + theme_void() +
  theme(legend.position = "none", aspect.ratio = 1)

mbt.pie.data <- mbt_counts %>% group_by(strain, categories.allowed) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>%
  mutate(per = total.count/sum(total.count) *100,
         pie.lab = ifelse(per > 2.5, paste0(categories.allowed, " ", round(per), "%"), NA))


#scatter plot freqs and durs
base.plot.mbt + annotation_custom(grob = inset.plot.mbt, xmin = 250, xmax = 950, ymin = 50, ymax = 90)

plot4.data.mbt <- data_freqs %>% filter((recording == "MBT") & (label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short" | label == "short-c" | label == "complex" | label == "trill" | label == "trill-c"))


inset.plot.mbt <- ggplotGrob(ggplot(plot4.data.mbt, aes(x = duration * 1000, y = m.freq/1000)) +
                           geom_point(size = 2, aes(colour = label), alpha = 0.5) +
                           scale_colour_manual(values=c("#BF383E", "#1F5694", "#13355C", "#3594FF", "#FF71B7","#FF1FE0", "#66BF9A", "#447F67")) +
                           xlab("Duration (ms)") +
                           ylab("Frequeny (kHz)") +
                           theme_bw() +
                           theme(legend.background = element_rect(colour = "transparent", 
                                                                  fill = "transparent"), 
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank(), 
                                 legend.position = c(.65,.85), legend.title = element_blank(), 
                                 legend.direction = "horizontal"))


base.plot.mbt <- ggplot(plot4.data.mbt, aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 2, aes(colour = strain), alpha = 0.3) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_bw() +
  scale_colour_grey() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())



#boxplot freqs
plot7 <- mbt_freqs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.freq/1000, colour = strain)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq/1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Frequency (kHz)") +
  scale_colour_grey()


#boxplot durs
plot8 <- mbt_durs %>% group_by(strain, label) %>%
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

plot6 <- ggplot(beh.usv.type.counts[beh.usv.type.counts$behavior == "Hover"|beh.usv.type.counts$behavior=="Retrieval",], aes(x = label, y = usv.count, colour = strain)) +
  geom_boxplot(fill = "white", position = position_dodge(0.8)) +
  geom_point(size = 2, aes(y = usv.count, group = strain), position = position_jitterdodge(0.2)) +
  facet_wrap(~behavior) +
  scale_colour_grey() +
  theme_classic() +
  theme(legend.position = "none")



#behavioral latency histogram
ggplot(mbt.beh.lats.data, aes(x = score, group = strain, fill = strain)) +
  geom_histogram(aes(y = ..count..), position = "dodge",
                 binwidth = 100)

#frequency histogram
mbt.freq.hist <- data_freqs %>% filter(recording == "MBT")
plot3 <- ggplot(mbt.freq.hist, aes(x = m.freq/1000, fill = strain)) + 
  geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100), position = "identity", colour = "black", alpha=0.5) +
  xlab("Frequency (kHz)") +
  ylab("Percent (%)")+
  scale_fill_grey() +
  theme_classic() +
  theme(legend.position = "none")
 

#duration histogram
mbt.dur.hist <- data_durs %>% filter(recording == "MBT" & duration < 0.4)
plot4 <- ggplot(mbt.dur.hist, aes(x = duration * 1000, fill = strain)) + 
  geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100), position = "identity", colour = "black", alpha=0.5) +
  xlab("Duration (ms)") +
  ylab("Percent (%)")+
  scale_fill_grey() +
  theme_classic()+
  theme(legend.position = "none")





