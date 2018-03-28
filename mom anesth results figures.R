library(cowplot)


plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, labels = "AUTO", ncol = 2, align = 'w')

#bar graph of major 4 call type counts by strain
plot1 <- momanesth_counts %>% 
  filter(categories.allowed == "flat" | categories.allowed == "flat-z" |categories.allowed == "short" | categories.allowed == "short-c" | categories.allowed == "short-ur") %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity", position = position_dodge(.8)) +
  theme_classic() +
  theme(legend.justification = c(0,0), legend.position = c(.75, .70),
        legend.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size = 14), legend.title = element_blank(),
        axis.text.x = element_text(size = 14)) +
  scale_fill_grey(labels = c("SD","WKY")) +
  xlab("USV Category") +
  ylab("USVs (n)")


#pie chart
label_colours <- read.csv(file.choose(), stringsAsFactors = F)
names(label_colours) <- c("categories.allowed","colour")

plot2 <- ggplot(ma.pie.data, aes(x="", y=per, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values= as.vector(colour.key$label.colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar("y", start=0) + 
  facet_wrap(~strain,labeller = as_labeller(c("SD"="SD","WK"="WKY"))) + 
  theme_void() +
  theme(legend.position = "right", aspect.ratio = 1, legend.title = element_blank(),
        legend.key.size = unit(.1,"cm"), legend.background = element_rect(colour = "transparent", fill = "transparent"))

ma.pie.data <- momanesth_counts %>% group_by(strain, categories.allowed) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>%
  mutate(per = total.count/sum(total.count) *100,
         pie.lab = ifelse(per > 2.5, paste0(categories.allowed, " ", round(per), "%"), NA))
colour.key <- left_join(ma.pie.data, label_colours)
colour.key <- colour.key %>% group_by(categories.allowed) %>%
  summarise(label.colour = unique(colour))


#line data counts by time bins
bins <- cut(start_data[start_data$recording == "MA",]$start.time,5,include.lowest=T, labels = c("1","2","3","4","5"))
line.data <- start_data %>% filter(recording == "MA") %>% 
  mutate(bin = bins) %>% 
  group_by(bin,strain,file.name) %>% 
  summarise(count = sum(no.counts))%>% 
  group_by(bin,strain) %>% 
  summarise(mean = mean(count), sem = sd(count)/sqrt(length(count)),count = length(count))

line.data <- ma.xtra %>% 
  mutate(time.group = ifelse(start.time < 180, "1",
                             ifelse(start.time >= 360 & start.time < 540, "2",
                                    ifelse(start.time > 720, "3", NA)))) %>%
  filter(!is.na(time.group)) %>%
  group_by(strain, time.group, rat.id, file.name) %>%
  summarise(time.count = length(label)) %>%
  group_by(strain, time.group) %>%
  summarise(m.time.group = mean(time.count), sem = sd(time.count)/sqrt(length(time.count))) %>%
  mutate(time.group = ifelse(time.group == "1", "1-3",
                             ifelse(time.group == "2", "7-9", "12-15")))
line.data$time.group <- as.factor(line.data$time.group)

plot3 <- ggplot(line.data, aes(x=time.group, y=m.time.group, group = strain, colour=strain)) + 
  geom_errorbar(aes(ymin=m.time.group-sem, ymax=m.time.group+sem), width=.1,size=.5) +
  geom_line(size=1) +
  geom_point(size=3) + 
  ylab("USVs (n)") +
  xlab("Time in test (3 min intervals)") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_colour_grey(start = 0, end = .7) +
  scale_x_discrete(labels = c("1"="start","2"="middle","3"="end"))

#scatter plot duration x frequency
plot4 <- base.plot + annotation_custom(grob = inset.plot, xmin = 100, xmax = 400, ymin = 55, ymax = 95)

plot4.data <- data_freqs %>% filter((recording == "MA") & (label == "flat" | label == "flat-z" | label == "short" | label == "short-c" | label == "short-ur"))
scatter.colour.key <- colour.key %>% filter(categories.allowed == "flat"|
                                              categories.allowed == "flat-z"|
                                              categories.allowed == "short"|
                                              categories.allowed == "short-c"|
                                              categories.allowed == "short-ur")

inset.plot <- ggplotGrob(ggplot(plot4.data, aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 2, aes(colour = label), alpha = 0.5) +
  scale_colour_manual(values=as.vector(scatter.colour.key$label.colour)) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_classic() +
  theme(legend.background = element_rect(colour = "transparent", fill = "transparent"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = c(.5,1), legend.title = element_blank(), 
        legend.direction = "horizontal",legend.text = element_text(size = 6),
        plot.background = element_rect(colour = "transparent", fill = "transparent"),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 5)))

  
base.plot <- ggplot(plot4.data, aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 2, aes(colour = strain), alpha = 0.3) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_classic() +
  scale_colour_grey() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#frequency boxplot
plot5.data <- data_freqs %>% filter(recording == "MA") %>% 
  group_by(strain, label, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label == "flat" | label == "flat-z" | label == "short" | label == "short-c" | label == "short-ur")


plot5 <- ggplot(plot5.data, aes(x = label, y = mean.freq/1000, colour = strain)) +
  theme_classic() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq/ 1000, group = strain),position=position_jitterdodge(0.2)) +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14)) +
  scale_colour_grey(start = 0, end = .6) +
  xlab("USV Category") +
  ylab("Frequency (kHz)")


#duration boxplot
plot6.data <- data_durs %>% filter(recording == "MA") %>% 
  group_by(strain, label, file.name, rat.id) %>% 
  summarise(mean.dur = mean(duration), sem = sd(duration)/sqrt(length(duration)),
            count = length(duration)) %>%
  filter(label == "flat" | label == "flat-z" | label == "short" | label == "short-c" | label == "short-ur")


plot6 <- ggplot(plot6.data, aes(x = label, y = mean.dur * 1000, colour = strain)) +
  theme_classic() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.dur*1000, group = strain),position=position_jitterdodge(0.2)) +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14)) +
  scale_colour_grey(start = 0, end = .6) +
  xlab("USV Category") +
  ylab("Duration (ms)")


##extra plots

#histogram of pup call start times

#by category by .5 min bins
ma.xtra <- start.rows %>% filter(recording == "MA")
bins <- cut(ma.xtra$start.time,30,include.lowest=T, labels = as.character(c(seq(from=0,to=14.5,by=.5))))
hist.flat.data <- ma.xtra %>% mutate(bin = bins) %>%
  group_by(bin, strain, file.name) %>%
  summarise(no.counts = length(label)) %>%
  group_by(bin,strain) %>% 
  summarise(count = mean(no.counts), sem = sd(no.counts)/sqrt(length(no.counts)))

ggplot(hist.flat.data, aes(x = bin, y = count, fill = strain)) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = count - sem, ymax = count + sem), position = position_dodge(.8)) +
  xlab("bin") +
  ylab("count") +
  scale_fill_grey() +
  theme_classic()

#by 1 min bins
bins <- cut(start_data[start_data$recording == "MA",]$start.time,15,include.lowest=T, labels = as.character(c(1:15)))
hist.start.data <- start_data %>% filter(recording == "MA") %>% 
  mutate(bin = bins) %>% 
  group_by(bin, strain, file.name) %>%
  summarise(no.counts = sum(no.counts)) %>%
  group_by(bin,strain) %>% 
  summarise(count = mean(no.counts), sem = sd(no.counts)/sqrt(length(no.counts)))

ggplot(hist.start.data, aes(x = bin, y = count, fill = strain)) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = count - sem, ymax = count + sem), position = position_dodge(.8)) +
  xlab("bin") +
  ylab("count") +
  scale_fill_grey() +
  theme_classic()

#by 5 min bins
bins <- cut(start_data[start_data$recording == "MA",]$start.time,3,include.lowest=T, labels = c("1st 5", "2nd 5", "3rd 5"))
hist.start.data <- start_data %>% filter(recording == "MA") %>% 
  mutate(bin = bins) %>% 
  group_by(bin, strain, file.name) %>%
  summarise(no.counts = sum(no.counts)) %>%
  group_by(bin,strain) %>% 
  summarise(count = mean(no.counts), sem = sd(no.counts)/sqrt(length(no.counts)))

ggplot(hist.start.data, aes(x = bin, y = count, fill = strain)) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = count - sem, ymax = count + sem), position = position_dodge(.8)) +
  xlab("bin") +
  ylab("count") +
  scale_fill_grey() +
  theme_classic()

#line plot of before and during grouping
line.group.data <- all.levels %>% filter(label == "flat" | label == "short") %>%
  group_by(strain, label, time.group) %>%
  summarise(mcount = mean(time.count), sem = sd(time.count)/sqrt(length(time.count)))

ggplot(line.group.data, aes(x = time.group, y = mcount, group = strain, colour = strain)) +
  geom_point(size=3) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin = mcount - sem, ymax = mcount + sem),width=.1,size=.5) +
  facet_wrap(~label)


#pie chart before and during grouping
ggplot(grouping.pie.data, aes(x="", y=per, fill=label))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values= as.vector(ck$label.colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 2) +
  coord_polar("y", start=0) + 
  facet_wrap(~strain * time.group) + 
  theme_void() +
  theme(legend.position = "none", aspect.ratio = 1)

grouping.pie.data <- all.levels %>%
  group_by(strain, label, time.group) %>%
  summarise(tot.label.count = sum(time.count)) %>%
  group_by(strain, time.group) %>%
  mutate(tot.time = sum(tot.label.count)) %>%
  group_by(strain, label, time.group) %>%
  filter(tot.time > 1) %>%
  mutate(per = tot.label.count/sum(tot.time) *100,
         pie.lab = ifelse(per > 2.5, paste0(label, " ", round(per), "%"), NA))

lc <- label_colours
names(lc) <- c("label","colour")
ck <- left_join(grouping.pie.data, lc)
ck <- ck %>% group_by(label) %>%
  summarise(label.colour = unique(colour))


#pie of three times... by strain
three.times.pie.data <- pru %>%
  group_by(strain, label, time.group) %>%
  summarise(strain.time.count = sum(time.count)) %>%
  group_by(strain, time.group) %>%
  mutate(tot.time = sum(strain.time.count)) %>%
  group_by(strain, label, time.group) %>%
  mutate(per = strain.time.count/tot.time *100,
         pie.lab = ifelse(per > 2.5, paste0(label, " ", round(per), "%"), NA))

three.times.pie.data %>%
ggplot(aes(x="", y=per, fill=label))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values= as.vector(ck$label.colour)) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 2) +
  coord_polar("y", start=0) + 
  facet_wrap(~time.group * strain) + 
  theme_void() +
  theme(legend.position = "none", aspect.ratio = 1)

lc <- label_colours
names(lc) <- c("label","colour")
ck <- left_join(three.times.pie.data, lc)
ck <- ck %>% group_by(label) %>%
  summarise(label.colour = unique(colour))