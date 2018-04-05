####msx3 figures



##bar plot of calls

bp.data <- count_frame_3 %>% filter(categories.allowed=="flat"|
                                      categories.allowed=="short"|
                                      categories.allowed=="trill"|
                                      categories.allowed=="trill-c") %>%
  group_by(strain, recording, categories.allowed) %>%
  summarise(m.count = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts)))

ggplot(bp.data, aes(x = categories.allowed, y = m.count, fill = strain)) +
  geom_bar(stat="identity", position=position_dodge(.8)) +
  geom_errorbar(aes(ymin = m.count-sem, ymax=m.count+sem),
                position = position_dodge(.8)) +
  facet_wrap(~recording)+
  scale_fill_grey()



#maternal behavior scores
library(reshape2)
bp.data <- melt(data=msx3_beh_scores,id = c("strain","rat.id", "file.name","recording"),
                     variable.name="behavior",value.name = "score")
#active scores
bp.active <- bp.data %>% filter(behavior=="retrieval"|behavior == "mouthing"|
                                  behavior=="corporal"|behavior=="anogenital"|
                                  behavior=="nest.building")

ggplot(bp.active[bp.active$strain=="WK",], aes(x=behavior, y=score, colour = recording))+
  theme_classic() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = score, group = recording),position=position_jitterdodge(0.2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14)) +
  xlab("Behavior") +
  ylab("# Behavior / 30min")

#latencies
bp.lat <- bp.data %>% filter(behavior=="lat.retrieve"|behavior == "lat.group"|
                   behavior=="lat.hover"|behavior=="lat.nurse")

ggplot(bp.lat[bp.lat$strain=="WK",], aes(x=behavior, y=score, colour = recording))+
  theme_classic() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = score, group = recording),position=position_jitterdodge(0.2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14)) +
  xlab("Behavior") +
  ylab("Latency (s)")

#durations
bp.dur <- bp.data %>% filter(behavior=="dur.hover"|behavior == "dur.nurse")

ggplot(bp.dur[bp.dur$strain=="WK",], aes(x=behavior, y=score, colour = recording))+
  theme_classic() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = score, group = recording),position=position_jitterdodge(0.2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14)) +
  xlab("Behavior") +
  ylab("Duration (s)")


#usvs emitted during behaviors
ubplot.data <- usv.beh.counts %>% filter(label=="flat"|label=="short"|
                                         label=="trill"|label=="trill-c"|
                                         label=="flat-mz") %>%
  group_by(strain, label, behavior, recording) %>%
  summarise(m.count = mean(usv.count), sem = sd(usv.count)/sqrt(length(usv.count)))

ggplot(ubplot.data, aes(x=label, y=m.count, fill=strain)) +
  geom_bar(stat="identity", position=position_dodge(.8)) +
  geom_errorbar(aes(ymin = m.count-sem, ymax=m.count+sem),
                position = position_dodge(.8)) +
  facet_wrap(~recording*behavior)+
  scale_fill_grey()

