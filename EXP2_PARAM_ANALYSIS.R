###analyzing msx-3 data

#count_frame_2
#check for UsV count outliers
source("http://goo.gl/UUyEzD")
outlierKD(ol.check.data[ol.check.data$strain=="WK" & ol.check.data$recording=="MSX3",], tot.count)


#remove outliers
count_frame_3 <- count_frame_2 %>% filter(!(file.name=="T0000399"|
                                              file.name=="T0000410"|
                                              file.name=="T0000291"))

library(nlme)
msx3_lme <- lme(total.counts~strain * categories.allowed * recording,
                random = ~1|rat.id/categories.allowed,
                data = count_frame_3)
anova.lme(msx3_lme)
library(lsmeans)
msx3.count.sum <- summary(lsmeans(msx3_lme, pairwise ~ strain|recording, adjust="Dunnet"))[["contrasts"]]
View(msx3.count.sum[msx3.count.sum$p.value<0.05,])
#diff in categories.allowed and strain x categories allowed
#first get total percents 
msx3_tot_call_percents <- count_frame_3 %>% group_by(categories.allowed) %>%
  summarise(tot.sum = sum(total.counts)) %>% ungroup() %>%
  mutate(percent = tot.sum / sum(tot.sum))
big.calls.frame <- msx3_tot_call_percents %>% filter(percent > 0.02) %>% 
  select(categories.allowed)
big.calls <- as.vector(big.calls.frame$categories.allowed)
#calls more than 2%: complex, flat-mz, trill-c, flat-z, trill, short, flat

#let's get percents
prop.data <- count_frame_2 %>% filter(categories.allowed %in% big.calls)


msx3_per_lme <- lme(rel.filecount ~ strain * categories.allowed * recording,
                    random = ~1|rat.id/categories.allowed,
                    data = prop.data)
anova(msx3_per_lme)
msx3.per.sum <- summary(lsmeans(msx3_per_lme, pairwise ~ strain|categories.allowed, adjust="Tukey"))[["contrasts"]]
View(msx3.per.sum[msx3.per.sum$p.value<0.05,])


#now time for behavior
msx3_beh_in <- read.csv("data/Exp2_MB_Scores_mariana.csv",stringsAsFactors = F)
#remove outliers
msx3_beh_scores_active <- msx3_beh_in %>% filter(!(file.name=="T0000399"|
                                              file.name=="T0000410"|
                                              file.name=="T0000291")) %>%
  mutate(group = ifelse(strain=="SD" & recording == "VEH", "1",
                        ifelse(strain=="SD" & recording == "MSX3","2",
                               ifelse(strain=="WK" & recording=="VEH","3","4"))))
msx3_beh_scores_active$group <- as.factor(msx3_beh_scores_active$group)

msx3_it_in <- read.csv("data/Exp2_MB_Scores.csv",stringsAsFactors = F)
msx3_beh_scores_latdur <- msx3_it_in %>% filter(!(file.name=="T0000399"|
                                                     file.name=="T0000410"|
                                                     file.name=="T0000291")) %>%
  mutate(group = ifelse(strain=="SD" & recording == "VEH", "1",
                        ifelse(strain=="SD" & recording == "MSX3","2",
                               ifelse(strain=="WK" & recording=="VEH","3","4"))))
msx3_beh_scores_latdur$group <- as.factor(msx3_beh_scores_latdur$group)


kruskal.test(retrieval ~ group, data=msx3_beh_scores_active)
kruskal.test(mouthing ~ group, data=msx3_beh_scores_active)
kruskal.test(corporal ~ group, data=msx3_beh_scores_active) #significant
kruskal.test(anogenital ~ group, data=msx3_beh_scores_active) #significant
kruskal.test(nest.building ~ group, data=msx3_beh_scores_active)

kruskal.test(lat.retrieve ~ group, data=msx3_beh_scores_active)

kruskal.test(lat.group ~ group, data=msx3_beh_scores_active) #signif

kruskal.test(dur.nurse ~ group, data=msx3_beh_scores_latdur)

#start with WKY
#active scores
msx3_beh_wk_a <- msx3_beh_scores_active %>% filter(strain=="WK")

wilcox.test(msx3_beh_wk_a$retrieval ~ msx3_beh_wk_a$recording)

wilcox.test(msx3_beh_wk_a$corporal ~ msx3_beh_wk_a$recording) #not signif

wilcox.test(msx3_beh_wk_a$anogenital ~ msx3_beh_wk_a$recording) #not signif

wilcox.test(msx3_beh_wk_a$nest.building ~ msx3_beh_wk_a$recording)

#latencies

wilcox.test(msx3_beh_wk_a$lat.retrieve ~ msx3_beh_wk_a$recording)

wilcox.test(msx3_beh_wk_a$lat.group ~ msx3_beh_wk_a$recording) #not signif

wilcox.test(msx3_beh_wk_ld$lat.hover ~ msx3_beh_wk_ld$recording)

wilcox.test(msx3_beh_wk_ld$lat.nurse ~ msx3_beh_wk_ld$recording)

#durations
msx3_beh_wk_ld <- msx3_beh_scores_latdur %>% filter(strain=="WK") 
wilcox.test(msx3_beh_wk_ld$dur.hover ~ msx3_beh_wk_ld$recording)

wilcox.test(msx3_beh_wk_ld$dur.nurse ~ msx3_beh_wk_ld$recording)

#non maternal
wilcox.test(msx3_beh_wk$crossing ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$rearing ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$self.groom ~ msx3_beh_wk$recording)


#try SD
msx3_beh_sd <- msx3_beh_scores %>% filter(strain=="SD")

wilcox.test(msx3_beh_sd$dur.nurse ~ msx3_beh_sd$recording)

#compare behaviors between strains within treatments
#vehicle
msx3_beh_veh <- msx3_beh_scores_active %>% filter(recording=="VEH")

wilcox.test(msx3_beh_veh$retrieval ~ msx3_beh_veh$strain)

wilcox.test(msx3_beh_veh$mouthing ~ msx3_beh_veh$strain)

wilcox.test(msx3_beh_veh$corporal ~ msx3_beh_veh$strain) #signif

wilcox.test(msx3_beh_veh$anogenital ~ msx3_beh_veh$strain) #signif

wilcox.test(msx3_beh_veh$nest.building ~ msx3_beh_veh$strain)

wilcox.test(msx3_beh_veh$lat.group ~ msx3_beh_veh$strain) #signif

#msx3
msx3_beh_msx3 <- msx3_beh_scores_active %>% filter(recording=="MSX3")

wilcox.test(msx3_beh_msx3$retrieval ~ msx3_beh_msx3$strain)

wilcox.test(msx3_beh_msx3$corporal ~ msx3_beh_msx3$strain)

wilcox.test(msx3_beh_msx3$anogenital ~ msx3_beh_msx3$strain) #signif

#compare MSX3 WK to VEH SD
msx3_beh_last <- msx3_beh_scores_active %>% filter((recording=="VEH"&strain=="SD")|
                                                     (recording=="MSX3"&strain=="WK"))
wilcox.test(msx3_beh_last$retrieval ~ msx3_beh_last$group)
wilcox.test(msx3_beh_last$corporal ~ msx3_beh_last$group)
wilcox.test(msx3_beh_last$anogenital ~ msx3_beh_last$group)
wilcox.test(msx3_beh_last$lat.group ~ msx3_beh_last$group)

#try ANOVA
#retrieval
retrieval.lme <- lme(retrieval ~ strain * recording,
               random = ~1|rat.id,
               data = msx3_beh_scores)
anova.lme(retrieval.lme)
retrieval.sum <- summary(lsmeans(retrieval.lme, pairwise ~ strain * recording, adjust = "Tukey"))[["contrasts"]]
View(retrieval.sum)

#corporal
corporal.lme <- lme(corporal ~ strain * recording,
                    random = ~1|rat.id,
                    data = msx3_beh_scores)
anova.lme(corporal.lme)
corporal.sum <- summary(lsmeans(corporal.lme, pairwise ~ strain, adjust = "Tukey"))[["contrasts"]]
View(corporal.sum)

#anogenital
anogen.lme <- lme(anogenital ~ strain * recording,
                    random = ~1|rat.id,
                    data = msx3_beh_scores)
anova.lme(anogen.lme)
anogen.sum <- summary(lsmeans(anogen.lme, pairwise ~ strain, adjust = "Tukey"))[["contrasts"]]
View(anogen.sum)

#nest building
nest.lme <- lme(nest.building ~ strain * recording,
                  random = ~1|rat.id,
                  data = msx3_beh_scores)
anova.lme(nest.lme)


#latency to group
l.group.lme <- lme(lat.group ~ strain * recording,
                  random = ~1|rat.id,
                  data = msx3_beh_scores)
anova.lme(l.group.lme)
l.group.sum <- summary(lsmeans(l.group.lme, pairwise ~ strain, adjust = "Tukey"))[["contrasts"]]
View(l.group.sum)

#latency to nurse
l.nurse.lme <- lme(lat.nurse ~ strain * recording,
                   random = ~1|rat.id,
                   data = msx3_beh_scores)
anova.lme(l.nurse.lme)


#duration hover over
dur.hover.lme <- lme(dur.hover ~ strain * recording,
                   random = ~1|rat.id,
                   data = msx3_beh_scores)
anova.lme(dur.hover.lme)


#duration nursing
dur.nurse.lme <- lme(dur.nurse ~ strain * recording,
                   random = ~1|rat.id,
                   data = msx3_beh_scores)
anova.lme(dur.nurse.lme)
dur.nurse.sum <- summary(lsmeans(dur.nurse.lme, pairwise ~ recording, adjust = "Tukey"))[["contrasts"]]
View(dur.nurse.sum)

###USVs emitted during behaviors
beh_times_data <- read.csv("data/Exp2_MBT_beh_times.csv",stringsAsFactors = F)
beh_times_data$behavior <- gsub("HKP","LKP",beh_times_data$behavior)


merged.data <- merge(start.rows, beh_times_data, by = c("file.name","rat.id","strain"))

beh.total.usv.counts <- merged.data %>% 
  filter(start.time > start & start.time < end) %>%
  group_by(behavior, label, file.name) %>%
  summarise(usv.count = length(label)) %>% filter(behavior == "Retrieval" |
                                                    behavior == "Hover" |
                                                    behavior == "LKP")
beh.label.expand <- expand.grid(file.name = unique(count_frame_3$file.name),
                                behavior = c("Retrieval","Hover","LKP"),
                                label = unique(count_frame_3$categories.allowed))
usv.beh.counts.filenames <- left_join(beh.label.expand, beh.total.usv.counts)
usv.beh.counts.join <- left_join(usv.beh.counts.filenames, file.name.key)
usv.beh.counts <- usv.beh.counts.join %>% 
  mutate(usv.count = ifelse(is.na(usv.count),0,usv.count))

##first try total counts per behavior
tot.usv.beh <- usv.beh.counts %>% group_by(strain, rat.id, file.name, behavior, recording) %>%
  summarise(sum.count = sum(usv.count))

tot.usv.beh.lme <- lme(sum.count ~ strain * behavior * recording,
                       random = ~1|rat.id/behavior,
                       data=tot.usv.beh)
anova.lme(tot.usv.beh.lme)



##determin maternal USVs
data_freqs_1 <- read.csv("data/Exp1_Filtered_All_Vars.csv", stringsAsFactors = F)
m.freq <- apply(data_freqs_1[,c(1:50)],1,function(y) mean(y))
data_freqs_readin <- mutate(data_freqs_1, m.freq = m.freq)
data_freqs_1 <- data_freqs_readin %>%
mutate(label = ifelse(duration < 0.012 & label !="short", "short",
ifelse(duration > 0.012 & label == "short", "flat",label))) %>%
filter(duration < 0.8 & duration > 0.002)

data_durs_1 <- read.csv("data/Exp1_Non_Overlap_Durations.csv", stringsAsFactors = F)
data_durs_1 <- data_durs_1 %>%
mutate(label = ifelse(duration < 0.012 & label !="short", "short",
ifelse(duration > 0.012 & label == "short", "flat",label))) %>%
filter(duration < 0.8 & duration > 0.002)

table.freqs <- data_freqs_1 %>% filter((recording=="MA"|recording=="Mpupiso"|
recording=="Fpupiso") & (label=="flat"|
label=="flat-z"|
label=="short")) %>%
group_by(strain, label) %>%
summarise(Min=round(min(m.freq)/1000,3),Mean=round(mean(m.freq)/1000,3),
SEM = round(sd(m.freq)/sqrt(length(m.freq))/1000,3),
five=round(quantile(m.freq,0.05,type=3)/1000,3),ninetyfive=
round(quantile(m.freq,0.95,type=3)/1000,3),
Max=round(max(m.freq)/1000,3))

table.durs <- data_durs_1 %>% filter((recording=="MA"|recording=="Mpupiso"|
recording=="Fpupiso") & (label=="flat"|
label=="flat-z"|
label=="short")) %>%
group_by(strain, label) %>%
summarise(Min=round(min(duration)*1000,3), Mean=round(mean(duration)*1000,3),
SEM = round(sd(duration)/sqrt(length(duration))*1000,3),
five=round(quantile(duration,0.05,type=3)*1000,3),ninetyfive=
round(quantile(duration,0.95,type=3)*1000,3),
Max=round(max(duration)*1000,3))




##distinguish between mom and pup
param.cutoffs <- data.frame(strain = table.freqs$strain,
                            label = table.freqs$label,
                            freq.min = table.freqs$five*1000,
                            freq.max = table.freqs$ninetyfive*1000,
                            dur.min = table.durs$five/1000,
                            dur.max = table.durs$ninetyfive/1000)

mom_calls_1 <- data_freqs %>%  
  left_join(param.cutoffs) %>%
  filter(m.freq < freq.min | m.freq > freq.max | duration < dur.min |
           duration > dur.max) %>%
  select(-freq.min, -freq.max, -dur.min, -dur.max)
mom_calls_2 <- data_freqs %>% filter(label=="trill"|label=="trill-c")

mom_freqs_joined <- rbind.data.frame(mom_calls_1, mom_calls_2)  
#frequency data of maternal USVs
mom_freqs_almost <- mom_freqs_joined %>% group_by(strain, rat.id, label) %>%
  summarise(mean.freq = mean(m.freq), sem.freq = sd(m.freq)/sqrt(length(m.freq)), 
            count = length(label))
mom_freqs <- left_join(mom_freqs_almost, file.name.key)

mom.freqs.lme <- lme(mean.freq ~ strain * label,
                     random = ~1|rat.id/label,
                     data = mom_freqs)
anova.lme(mom.freqs.lme)

#duration data of maternal usvs
mom_durs_1 <- data_durs %>% 
  anti_join(mom_freqs_joined, by=c("unique.id"))
mom_durs_2 <- mom_durs_1 %>%
  left_join(param.cutoffs) %>%
  filter((duration < dur.min |duration > dur.max)) %>%
  select(-freq.min, -freq.max, -dur.min, -dur.max, -unique.id)
mom_freqs_joined_fix <- mom_freqs_joined %>%
  select(file.name, label, strain, rat.id, recording, duration, start.time)

mom_durs_almost <- rbind.data.frame(mom_durs_2, mom_freqs_joined_fix) 
mom_durs_almost2 <- mom_durs_almost %>% group_by(strain, rat.id, label) %>%
  summarise(mean.dur = mean(duration), count = length(duration))
mom_durs <- left_join(mom_durs_almost2, file.name.key)

mom.durs.lme <- lme(mean.dur ~ strain * label,
                    random = ~1|rat.id/label,
                    data=mom_durs_data)
anova.lme(mom.durs.lme)
mom.durs.sum <- summary(lsmeans(mom.durs.lme, pairwise~strain|label, adjust="Tukey"))[["contrasts"]]
View(mom.durs.sum[mom.durs.sum$p.value<0.05,])

#count data of maternal USVs
mom_calls_1_adjust <- mom_durs_almost %>% select(-duration)

mom_calls_3 <- start.rows %>% 
  anti_join(mom_durs_almost, by = c("start.time")) 

mom_calls_4 <- mom_calls_3 %>%
  filter(label=="trill"|label=="trill-c") %>%
  select(-unique.id)

mom_counts_1 <- rbind.data.frame(mom_calls_1_adjust, mom_calls_4) 
mom_counts_2 <- mom_counts_1 %>% group_by(strain, label, rat.id) %>%
  summarise(usv.count= length(label))

grid <- expand.grid(rat.id = unique(count_frame_3$rat.id),
                    label = c("flat","flat-z","short","trill","trill-c"))
mom_counts_3 <- left_join(grid, mom_counts_2)
mom_counts.almost <- mom_counts_3 %>% mutate(strain = ifelse(str_detect(rat.id,"SD")==T,
                                                      "SD","WK"),
                                      usv.count = ifelse(is.na(usv.count),0,usv.count))
mom_counts <- left_join(mom_counts.almost,file.name.key)

mom.counts.lme <- lme(usv.count ~ strain * label*recording,
                      random = ~1|rat.id/label,
                      data = mom_counts)
anova.lme(mom.counts.lme)
mom.counts.sum <- summary(lsmeans(mom.counts.lme, pairwise~strain|label, adjust="Tukey"))[["contrasts"]]
View(mom.counts.sum[mom.counts.sum$p.value<0.05,])

##mom counts by category



#start times!!
mom_start_1 <- mom_counts_1

