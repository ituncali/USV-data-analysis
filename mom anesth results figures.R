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
  theme_bw() +
  theme(legend.justification = c(0,1), legend.position = c(.75, .75),legend.background = element_rect(colour = "transparent", fill = "transparent")) +
  scale_fill_grey() +
  xlab("USV Category") +
  ylab("Ultrasonic Vocalizations (n)")


#pie chart
plot2 <- ggplot(ma.pie.data, aes(x="", y=per, fill=categories.allowed))+
  geom_bar(width = 1, stat = "identity", alpha = .5, colour = "black") +
  scale_fill_manual(values=c("#BF383E", "#73BF5C", "#1F5694", "#13355C", "#3594FF", "#FFC19E", "#94705C", "#FF71B7", "#FF1FE0", "#E51CCA", "#BF5BA8", "#FFB3FA", "#5C0B51", "#BFBE3F", "#949339", "#5C5B1E", "#375C2C")) +
  geom_text(aes(label = pie.lab), position = position_stack(vjust = 0.5), size = 2) +
  coord_polar("y", start=0) + facet_wrap(~strain) + theme_void() +
  theme(legend.position = "none", aspect.ratio = 1)
  
ma.pie.data <- momanesth_counts %>% group_by(strain, categories.allowed) %>%
  summarise(total.count = sum(total.counts)) %>%
  filter(total.count > 1) %>%
  mutate(per = total.count/sum(total.count) *100,
         pie.lab = ifelse(per > 2.5, paste0(categories.allowed, " ", round(per), "%"), NA))



#line data counts by time bins
bins <- cut(start_data[start_data$recording == "MA",]$start.time,5,include.lowest=T, labels = c("1","2","3","4","5"))
line.data <- start_data %>% filter(recording == "MA") %>% 
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
  xlab("Test Duration (3min)") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_grey()

#scatter plot duration x frequency
plot4 <- base.plot + annotation_custom(grob = inset.plot, xmin = 100, xmax = 400, ymin = 55, ymax = 90)

plot4.data <- data_freqs %>% filter((recording == "MA") & (label == "flat" | label == "flat-z" | label == "short" | label == "short-c" | label == "short-ur"))

inset.plot <- ggplotGrob(ggplot(plot4.data, aes(x = duration * 1000, y = m.freq/1000)) +
                            geom_point(size = 2, aes(colour = label), alpha = 0.5) +
                            scale_colour_manual(values=c("#1F5694", "#3594FF", "#FF71B7", "#FF1FE0", "#5C0B51")) +
                            xlab("Duration (ms)") +
                            ylab("Frequeny (kHz)") +
                            theme_bw() +
                            theme(legend.background = element_rect(colour = "transparent", fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.65,.85), legend.title = element_blank(), legend.direction = "horizontal"))

  
base.plot <- ggplot(plot4.data, aes(x = duration * 1000, y = m.freq/1000)) +
  geom_point(size = 2, aes(colour = strain), alpha = 0.3) +
  xlab("Duration (ms)") +
  ylab("Frequeny (kHz)") +
  theme_bw() +
  scale_colour_grey() +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#frequency boxplot
plot5 <- momanesth_freqs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.freq/1000, colour = strain)) +
  theme_bw() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.freq/ 1000, group = strain),position=position_jitterdodge(0.2)) +
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_grey() +
  xlab("USV Category") +
  ylab("Frequency (kHz)")


#duration boxplot
plot6 <- momanesth_durs %>% group_by(strain, label) %>%
  ggplot(aes(x = label, y = mean.dur * 1000, colour = strain)) +
  theme_bw() +
  geom_boxplot(fill = "white", position=position_dodge(0.8)) +
  geom_point(aes(y = mean.dur*1000, group = strain),position=position_jitterdodge(0.2)) +
  theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_grey() +
  xlab("USV Category") +
  ylab("Duration (ms)")


##extra plots

#histogram of pup call start times
#data
start_data <- cbind(data_counts, 
                    no.counts = mapply(function(q) sum(str_count(q,categories.allowed.search)),data_counts$label))
rows.to.rep <- start_data %>% filter(no.counts > 1)
#total is sum(rows.to.rep$no.counts) = 3086 so this should be the same after replicating rows
xtra.list <- mapply(function(q) str_extract_all(q, pattern = regex("[a-z]+(-[a-z]+)|[a-z]+")),
                    rows.to.rep$label)
xtra.notlist <- unlist(xtra.list)
xtra.notlist.fix <- str_replace_all(xtra.notlist, 
                                    pattern = c("downward" = "downward ramp",
                                                "upward" = "upward ramp",
                                                "up$" = "step up",
                                                "down$" = "step down",
                                                "inverted$" = "inverted-u"))                           
#first remove the rows that do not contain calls
xtra.rows <- unlist(mapply(function(q) str_subset(q, categories.allowed.search), xtra.notlist.fix))
xtra.rows.added <- cbind(label = xtra.rows,
                         rat.id = unlist(mapply(function(q,x) c(rep(q,x)),
                                                q = rows.to.rep$rat.id, x = rows.to.rep$no.counts)),
                         strain = unlist(mapply(function(q,x) c(rep(q,x)),
                                                q = rows.to.rep$strain, x = rows.to.rep$no.counts)),
                         start.time = unlist(mapply(function(q,x) c(rep(q,x)),
                                                    q = rows.to.rep$start.time, x = rows.to.rep$no.counts)),
                         file.name = unlist(mapply(function(q,x) c(rep(q,x)),
                                                   q = rows.to.rep$file.name, x = rows.to.rep$no.counts)),
                         unique.id = unlist(mapply(function(q,x) c(rep(q,x)),
                                                   q = rows.to.rep$unique.id, x = rows.to.rep$no.counts)),
                         recording = unlist(mapply(function(q,x) c(rep(q,x)),
                                                   q = rows.to.rep$recording, x = rows.to.rep$no.counts)))
start.category.data <- rbind(start_data[start_data$no.counts < 2,-8],xtra.rows.added)

#by category by 1 min bins
ggplot(start.category.data[start.category.data$label == "flat",], aes(x = start.time, fill = strain)) + 
  geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                 position = "dodge", colour = "black", alpha = 0.5) +
  scale_fill_grey() +
  xlab("Start Time of Flat USVs") +
  ylab("Percent %") +
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
