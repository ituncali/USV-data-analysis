library(cowplot)

plot_grid(plot1, plot3, plot4, labels = "AUTO", ncol = 2, align = 'w')


#boxplot
plot1 <- maternal_counts %>% 
  filter(categories.allowed == "flat" | categories.allowed == "short" |categories.allowed == "trill-c") %>% 
  group_by(strain, recording, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = recording, colour = recording)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity", position = position_dodge(.8), fill = "white") +
  theme_classic() +
  theme(legend.justification = c(0,0), legend.position = c(.24, .7), 
        legend.background = element_rect(colour = "transparent", fill = "transparent"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank()) +
  xlab("USV Category") +
  ylab("Ultrasonic Vocalizations (n)") +
  facet_wrap(~strain)

#mom alone and pups separated line graph
plot2 <- maternal_counts %>% group_by(strain, recording, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>%
  mutate(point.lab = ifelse(recording == "PupsSep" & mean > 5, paste0(categories.allowed), NA)) %>%
  ggplot(aes(x = recording, y = mean, group = categories.allowed, colour = categories.allowed)) + 
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = .1, size = .5, position = position_dodge(.5)) + 
  geom_line(size=1,position = position_dodge(.5)) + geom_point(size=3, position = position_dodge(.5)) + facet_wrap(~strain) +
  scale_colour_manual(values=c("#BF383E", "#73BF5C", "#1F5694", "#13355C", "#3594FF", "#FF4540", 
                               "#FFC19E", "#94705C", "#FF71B7", "#FF1FE0", "#E51CCA", "#BF5BA8", 
                               "#FFB3FA", "#5C0B51", "#BFBE3F", "#949339", "#66BF9A", "#447F67", 
                               "#224033", "#5C5B1E", "#375C2C")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(aes(label=point.lab),hjust=0, vjust=-2, colour = "black")


#boxplot of freqs
plot3 <- maternal_freqs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.freq/1000, colour = strain)) + 
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq/1000, group = strain), size=2, position = position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Frequency (kHz)") +
  scale_colour_grey(start = 0, end = .6) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14))


#boxplot of durs
plot4 <- maternal_durs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.dur * 1000, colour = strain)) + 
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.dur * 1000, group = strain), size=2, position = position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Duration (ms)") +
  scale_colour_grey(start = 0, end = .6) + 
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14))

  