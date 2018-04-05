####msx3 figures



##bar plot of calls

bp.data <- count_frame_3 %>% filter(categories.allowed=="flat"|
                                      categories.allowed=="short"|
                                      categories.allowed=="trill"|
                                      categories.allowed=="trill-c") %>%
  group_by(strain, recording, categories.allowed) %>%
  summarise(m.count = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>%
  mutate(group = ifelse(strain=="SD" & recording=="VEH", "1",
                        ifelse(strain=="SD" & recording=="MSX3", "2",
                               ifelse(strain=="WK" & recording=="VEH","3","4"))))

ggplot(bp.data, aes(x = categories.allowed, y = m.count, fill = group)) +
  geom_bar(stat="identity", position=position_dodge(.8)) +
  geom_errorbar(aes(ymin = m.count-sem, ymax=m.count+sem),
                position = position_dodge(.8)) +
  scale_fill_grey(labels=c("SD VEH","SD MSX-3","WKY VEH","WKY MSX-3")) +
  theme_classic() +
  theme(legend.justification = c(0,0), legend.position = c(.6,.6)) 


#pie charts of calls
plot2 <- ggplot(mbt.pie.data, aes(x="", y=per, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values=as.vector(colour.key$label.colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar("y", start=0) + 
  facet_wrap(~strain * recording) + 
  theme_void() +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        aspect.ratio = 1, 
        legend.title = element_blank(),
        legend.key.size = unit(.1,"cm"),
        legend.background = element_rect(colour = "transparent", fill = "transparent"))

mbt.pie.data <- count_frame_3 %>% group_by(strain, categories.allowed, recording) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>% group_by(strain, recording) %>%
  mutate(per = total.count/sum(total.count) *100,
         pie.lab = ifelse(per > 2.5, paste0(categories.allowed, " ", round(per), "%"), NA))
label_colours <- read.csv("data/usv_label_colours.csv", stringsAsFactors = F)
names(label_colours) <- c("categories.allowed", "colour")
colour.key <- left_join(mbt.pie.data, label_colours)
colour.key <- colour.key %>% group_by(categories.allowed) %>%
  summarise(label.colour = unique(colour))



#maternal behavior scores
library(reshape2)
bp.data <- melt(data=msx3_beh_scores,id = c("strain","rat.id", "file.name","recording"),
                     variable.name="behavior",value.name = "score")
#active scores
bp.active <- bp.data %>% filter(behavior=="retrieval"|behavior == "mouthing"|
                                  behavior=="corporal"|behavior=="anogenital"|
                                  behavior=="nest.building") %>%
  group_by(strain, behavior, recording) %>%
  summarise(med.score = median(score)) %>%
  mutate(group = ifelse(strain=="SD" & recording == "VEH", "1",
                        ifelse(strain=="SD" & recording == "MSX3","2",
                               ifelse(strain=="WK" & recording=="VEH","3","4"))))

ggplot(bp.active, aes(x=behavior, y=med.score, fill = group))+
  theme_classic() +
  geom_bar(stat = "identity", position=position_dodge(0.8)) +
 # geom_point(aes(y = score, group = recording),position=position_jitterdodge(0.2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14),
        legend.title=element_blank()) +
  xlab("Behavior") +
  ylab("# Behavior / 30min") +
  scale_fill_grey(labels = c("SD VEH","SD MSX-3","WKY VEH","WKY MSX-3"))

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

ggplot(ubplot.data, aes(x=behavior, y=m.count, fill=label)) +
  geom_bar(stat="identity") +
  #geom_errorbar(aes(ymin = m.count-sem, ymax=m.count+sem),position = position_dodge(.8)) +
  facet_wrap(~strain*recording)


#pie charts of USVs emitted during behaviors
plot7.data <- usv.beh.counts %>% group_by(strain, behavior, label) %>%
  summarise(tot.count = sum(usv.count)) %>%
  group_by(strain, behavior) %>%
  mutate(per = tot.count/sum(tot.count)*100, 
         pie.lab = ifelse(per > 2.5, paste0(label, " ", round(per),
                                            "%"), NA),
         behavior_f = factor(behavior, levels = c("Retrieval","Hover","LKP")),
         strain_f = factor(strain, levels = c("SD","WK")))

label_colours <- read.csv("data/usv_label_colours.csv", stringsAsFactors = F)

colour.key <- label_colours %>% filter(label=="flat"|label=="short"|
                                         label=="trill"|label=="trill-c"|
                                         label=="flat-mz")


plot7 <- ggplot(plot7.data, aes(x="", y=per, fill=label))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values= as.vector(colour.key$colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar("y", start=0) + 
  facet_wrap(~strain_f*behavior_f) + 
  theme_void() +
  theme(legend.position = "right", aspect.ratio = 1, legend.title = element_blank(),
        legend.key.size = unit(.1,"cm"), legend.background = element_rect(colour = "transparent", fill = "transparent"))



