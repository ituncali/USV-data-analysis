library(gridExtra)
library(cowplot)
###pup isolation

#create layout
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2)

plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, labels = "AUTO", ncol = 2, align = 'w')

#bar graph of major 4 call type counts by strain
plot1 <- pup_counts %>% 
  filter(categories.allowed == "flat" | categories.allowed == "flat-z" |categories.allowed == "flat-mz" | categories.allowed == "short") %>% 
  group_by(strain, categories.allowed) %>% 
  summarise(mean = mean(total.counts), sem = sd(total.counts)/sqrt(length(total.counts))) %>% 
  ggplot(aes(x = categories.allowed, y = mean, group = strain, fill = strain)) + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.1,size=.5, position = position_dodge(.8)) + 
  geom_bar(stat = "identity", position = position_dodge(.8)) +
  theme_bw() +
  theme(legend.justification = c(0,0), legend.position = c(.85, .76), legend.background = element_rect(colour = "transparent", fill = "transparent")) +
  scale_fill_grey() +
  xlab("USV Category") +
  ylab("Ultrasonic Vocalizations (n)")


#pie chart
plot2 <- pup_counts %>% group_by(strain, categories.allowed) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>%
  mutate(per = round(total.count/sum(total.count) *100, digits = 0),
         pie.lab = ifelse(per > 2, paste0(categories.allowed, " ", per, "%"), NA)) %>%
  ggplot(aes(x="", y=per, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values=c("#BF383E", "#73BF5C", "#1F5694", "#13355C", "#3594FF", "#FFC19E", "#FF71B7", "#FF1FE0", "#E51CCA", "#BF5BA8", "#FFB3FA", "#5C0B51", "#BFBE3F", "#949339", "#5C5B1E", "#375C2C")) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 2) +
  coord_polar("y", start=0) + facet_wrap(~strain) + theme_void() +
  theme(legend.position = "none", aspect.ratio = 1)

#line plot of number of calls per 1 min time bin
allowed.categories <- c("flat", "flat-z", "flat-mz", "short", "short-su", "short-sd",
                        "short-ur", "short-dr", "short-c", "complex", "upward ramp",
                        "downward ramp", "step up", "step down", "multi-step", "multi-step-s",
                        "trill", "trill-c", "trill-f", "inverted-U", "unclear")

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
  ylab("Ultrasonic Vocalizations (n)") +
  xlab("Test Duration (min)") +
  theme_bw() +
  theme(legend.justification=c(2,0),
        legend.position=c(1,.6)) +
  scale_colour_grey()



#scatter plot duration vs. frequency
plot4 <- base.plot + annotation_custom(grob = inset.plot, xmin = 140, xmax = 340, ymin = 65, ymax = 90)
  
plot4.data <- data_freqs %>% filter((recording == "Mpupiso" | recording == "Fpupiso") & (label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short"))
  #mutate(label2 = str_replace(label, pattern = "-.*", replacement = "")),
   
inset.plot <- ggplotGrob(ggplot(plot4.data, aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 2, aes(colour = label), alpha = 0.5) +
  scale_colour_manual(values=c("#1F5694", "#3594FF", "#13355C", "#FF71B7")) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_bw() +
  theme(legend.background = element_rect(colour = "transparent", fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.65,.85), legend.title = element_blank(), legend.direction = "horizontal"))

base.plot <- data_freqs %>% filter((recording == "Mpupiso" | recording == "Fpupiso") & (label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short")) %>%
  ggplot(aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 2, aes(colour = strain), alpha = 0.3) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_bw() +
  scale_colour_grey() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())



#boxplot of pup freqs
plot5 <- pup_freqs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.freq/1000, colour = strain)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq/1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Frequency (kHz)") +
  scale_colour_grey()


#boxplot of pup durs
plot6 <- pup_durs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.dur * 1000, colour = strain)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.dur * 1000, group = strain),position=position_jitterdodge(0.2)) +
  xlab("USV Category") +
  ylab("Duration (ms)") +
  scale_colour_grey()

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

