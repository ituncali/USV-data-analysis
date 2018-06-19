beh_times_data <- read.csv("data/Exp2_MBT_beh_times_2.csv",stringsAsFactors = F)

wkveh_startrows <- start.rows %>% filter(strain == "WK" & recording=="VEH")

merged.data <- merge(wkveh_startrows, beh_times_data, by = c("rat.id"))

beh.total.usv.counts <- merged.data %>% 
  filter(start.time > start & start.time < end) %>%
  group_by(behavior, rat.id, label) %>%
  summarise(usv.count = length(label))

rets <- c("ret.1","ret.2","ret.3","ret.4","ret.5","ret.6","ret.7","ret.8")

retrieval.total.counts <- beh.total.usv.counts %>%
  group_by(rat.id, label) %>%
  filter(behavior %in% rets) %>%
  summarise(usv.count = sum(usv.count)) %>%
  ungroup() %>%
  mutate(behavior = rep("retrieval",length(retrieval.total.counts$rat.id)))

beh.usv.counts.plot.data <- rbind.data.frame(beh.total.usv.counts, retrieval.total.counts) %>%
  ungroup() %>%
  filter(!(behavior %in% rets))
  
  
beh.usv.counts.plot.data %>%
  group_by(behavior, label) %>%
  summarise(m.count = mean(usv.count), 
            sd.count = ifelse(length(usv.count)>1, 
                              sd(usv.count)/sqrt(length(usv.count)), 0)) %>%
ggplot(aes(x = label, y = m.count, fill=label)) +
  geom_bar(stat = "identity", position=position_dodge(.8)) +
  geom_errorbar(aes(ymin=m.count-sd.count, ymax=m.count+sd.count),
                position=position_dodge(.8)) +
  #scale_colour_grey(start = 0, end = .6, labels = c("SD","WKY")) +
  theme_classic() +
  #theme(legend.position = "none") +
  xlab("Behavior") +
  ylab("# USVs") +
  facet_wrap(~behavior)
