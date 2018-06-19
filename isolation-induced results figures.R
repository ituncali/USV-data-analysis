library(gridExtra)
library(cowplot)
###pup isolation

#create layout
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2)

plot_grid(plot1, plot2, plot3, plot4, labels = "AUTO", ncol = 2, align = 'w')

#bar graph of major 4 call type counts by strain
plot1 <- pup_counts %>% 
  #filter(categories.allowed == "flat" | categories.allowed == "flat-z" |categories.allowed == "flat-mz" | categories.allowed == "short") %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity", position = position_dodge(.8), colour="black",
           width=.7) +
  theme_classic() +
  theme(legend.justification = c(0,0), legend.position = c(.6, .60), 
        legend.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size = 17), legend.title = element_blank(),
        axis.text = element_text(size = 17),
        axis.title = element_text(size=14)) +
  scale_fill_grey(labels = c("SD (n=8)","WKY (n=7)")) +
  xlab("USV Category") +
  ylab("# USVs / 5min")


#pie chart
label_colours <- read.csv("data/usv_label_colours.csv", stringsAsFactors = F)
names(label_colours) <- c("categories.allowed","colour")
plot2.data <- count_frame_2 %>% filter(recording=="Fpupiso"|recording=="Mpupiso") %>% 
  group_by(strain, categories.allowed) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>%
  mutate(per = round(total.count/sum(total.count) *100, digits = 0),
         pie.lab = ifelse(per > 2, paste0(categories.allowed, " ", per, "%"), NA))
colour.key <- left_join(plot2.data, label_colours)
colour.key <- colour.key %>% group_by(categories.allowed) %>%
  summarise(label.colour = unique(colour))

plot2 <- ggplot(plot2.data, aes(x="", y=per, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity", alpha = 0.5, colour = "black") +
  scale_fill_manual(values= as.vector(colour.key$label.colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = .5), size = 6) +
  coord_polar("y", start=0) + 
  facet_wrap(~strain, labeller = as_labeller(c("SD"="SD","WK"="WKY"))) + 
  theme_void() +
  theme(legend.position = "bottom", aspect.ratio = 1, legend.title = element_blank(),
        legend.key.size = unit(.1,"cm"), legend.background = element_rect(colour = "transparent", 
                                                                          fill = "transparent"),
        strip.text=element_text(size=24))

#line plot of number of calls per 1 min time bin
allowed.categories <- c("flat", "flat-z", "flat-mz", "short", "short-su", "short-sd",
                        "short-ur", "short-dr", "short-c", "complex", "upward ramp",
                        "downward ramp", "step up", "step down", "multi-step", "multi-step-s",
                        "trill", "trill-c", "trill-f", "inverted-u", "unclear")

categories.allowed.search <- paste0(allowed.categories, "(?!-)")
library(stringr)
start_data <- cbind(data_counts, 
                    no.counts = mapply(function(q) sum(str_count(q,categories.allowed.search)),data_counts$label))

bins <- cut(start_data[start_data$recording == "Mpupiso" | start_data$recording == "Fpupiso",]$start.time,5,include.lowest=T, labels = c("1","2","3","4","5"))
line.data <- start_data %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  mutate(bin = bins) %>% 
  group_by(bin,strain,file.name) %>% 
  summarise(count = sum(no.counts))%>% 
  group_by(bin,strain) %>% 
  summarise(mean = mean(count), sem = sd(count)/sqrt(length(count)),count = length(count))

plot3 <- ggplot(line.data, aes(x=bin, y=mean, group = strain, colour=strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5) +
  geom_line(size=1) +
  geom_point(size=3) + 
  ylab("# USVs / 1min") +
  xlab("Test Duration (min)") +
  theme_classic() +
  theme(legend.position="none",
        axis.text = element_text(size = 17),
        axis.title = element_text(size=14)) +
  scale_colour_grey(start = 0, end = .7)



#scatter plot duration vs. frequency
plot4 <- grid.arrange(dur.hist, empty, main.plot, freq.hist, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.line = element_blank())

pup.freq.hist <- data_freqs %>% filter(recording == "Mpupiso"|recording=="Fpupiso")
freq.hist <- ggplot(pup.freq.hist, aes(x = m.freq/1000, fill = strain)) + 
  geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100), position = "identity", colour = "black", alpha=0.5) +
  ylab("(%)")+
  scale_fill_grey() +
  theme_classic() +
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.text.x=element_text(size=20),
        axis.title.x=element_text(size=14)) +
  coord_flip()

pup.dur.hist <- data_durs %>% filter((recording == "Mpupiso"|recording=="Fpupiso") & duration < 0.4)
dur.hist <- ggplot(pup.dur.hist, aes(x = duration * 1000, fill = strain)) + 
  geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100), 
                 position = "identity", colour = "black", alpha=0.5) +
  ylab("(%)")+
  scale_fill_grey() +
  theme_classic()+
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=14))


main.plot <- base.plot + annotation_custom(grob = inset.plot, xmin = 155, xmax = 340, ymin = 50, ymax = 95)
  
plot4.data <- data_freqs %>% filter((recording == "Mpupiso" | recording == "Fpupiso") & (label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short"))
  #mutate(label2 = str_replace(label, pattern = "-.*", replacement = "")),
scatter.colour.key <- colour.key %>% filter(categories.allowed == "flat"|
                                              categories.allowed == "flat-mz"|
                                              categories.allowed == "flat-z"|
                                              categories.allowed == "short")
   
inset.plot <- ggplotGrob(ggplot(plot4.data, aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 1, aes(colour = label), alpha = 0.5) +
  scale_colour_manual(values=as.vector(scatter.colour.key$label.colour)) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_classic() +
  theme(legend.background = element_rect(colour = "transparent", fill = "transparent"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.5,1), legend.title = element_blank(), 
        legend.direction = "horizontal", legend.text = element_text(size = 6),
        plot.background = element_rect(colour = "transparent", fill = "transparent"),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5)))

base.plot <- data_freqs %>% filter((recording == "Mpupiso" | recording == "Fpupiso") & (label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short")) %>%
  ggplot(aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 2, aes(colour = strain), alpha = 0.3) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_classic() +
  scale_colour_grey() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text=element_text(size=20),
        axis.title=element_text(size=14))



#boxplot of pup freqs
plot5 <- pup_freqs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.freq/1000, colour = strain)) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14)) +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq/1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Frequency (kHz)") +
  scale_colour_grey(start = 0, end = .6)


#boxplot of pup durs
plot6 <- pup_durs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.dur * 1000, colour = strain)) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 14)) +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.dur * 1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Duration (ms)") +
  scale_colour_grey(start = 0, end = .6)

###extra graphs

#histogram of pup freqs
data_freqs %>% 
  filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  ggplot(aes(x = m.freq/1000, fill = strain)) + 
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                 position = "identity", colour = "black", alpha = 0.5) +
  xlab("Frequency (kHz)") +
  ylab("Percent (%)") +
  scale_fill_grey() +
  theme_classic()

data_durs %>% 
  filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  ggplot(aes(x = duration*1000, fill = strain)) + 
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                 position = "identity", colour = "black", alpha = 0.5) +
  xlab("Duration (ms)") +
  ylab("Percent (%)") +
  scale_fill_grey() +
  theme_classic()


##horizontal duration bar
dur.bar <- pup_durs %>% group_by(strain, label, rat.id) %>% 
  summarise(m.litter.dur = mean(mean.dur)) %>% 
  group_by(strain, label) %>% 
  summarise(m.d = mean(m.litter.dur)*1000, sem = 1000*sd(m.litter.dur)/sqrt(length(m.litter.dur))) %>% 
  ggplot(aes(x=label, y=m.d, fill=strain)) +
  geom_errorbar(aes(ymin=m.d-sem, ymax=m.d+sem),position=position_dodge(.8),
                width=.1,size=.5)+
  geom_bar(stat="identity", position=position_dodge(.8), width=.7, colour="black")+
  #coord_flip() +
  theme_classic() +
  theme(legend.justification = c(0,0),
        legend.position = c(.8,.6),
        axis.text=element_text(size=20),
        axis.title=element_text(size=14)) +
  scale_fill_grey() +
  ylab("Duration (ms)")+
  xlab("USV Category")


#frequency bars
freq.bar <- pup_freqs %>% #group_by(strain, label, rat.id) %>% 
  #summarise(m.litter.freq = mean(mean.freq)) %>%
  group_by(strain, label) %>%
  summarise(m.f = mean(mean.freq)/1000, sem=(sd(mean.freq)/sqrt(length(mean.freq)))/1000) %>%
  ggplot(aes(x=label, y=m.f, fill=strain)) +
  geom_errorbar(aes(ymin=m.f-sem, ymax=m.f+sem), position=position_dodge(.8),
                width=.1, size=.5)+
  geom_bar(stat="identity",position=position_dodge(.8), width=.7, colour="black")+
  theme_classic() +
  theme(legend.position="none",
        axis.text=element_text(size=20),
        axis.title=element_text(size=14))+
  scale_fill_grey() +
  ylab("Peak Frequency (kHz)")+
  xlab("USV Category")



