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
msx3_beh_in <- read.csv("data/Exp2_MB_Scores.csv",stringsAsFactors = F)
#remove outliers
msx3_beh_scores <- msx3_beh_in %>% filter(!(file.name=="T0000399"|
                                              file.name=="T0000410"|
                                              file.name=="T0000291"))

#start with WKY
#active scores
msx3_beh_wk <- msx3_beh_scores %>% filter(strain=="WK")

wilcox.test(msx3_beh_wk$retrieval ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$corporal ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$anogenital ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$nest.building ~ msx3_beh_wk$recording)

#latencies
wilcox.test(msx3_beh_wk$lat.retrieve ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$lat.group ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$lat.hover ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$lat.nurse ~ msx3_beh_wk$recording)

#durations
wilcox.test(msx3_beh_wk$dur.hover ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$dur.nurse ~ msx3_beh_wk$recording)

#non maternal
wilcox.test(msx3_beh_wk$crossing ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$rearing ~ msx3_beh_wk$recording)

wilcox.test(msx3_beh_wk$self.groom ~ msx3_beh_wk$recording)


#try SD
msx3_beh_sd <- msx3_beh_scores %>% filter(strain=="SD")

wilcox.test(msx3_beh_sd$dur.nurse ~ msx3_beh_sd$recording)

#compare behaviors between strains within treatments
#vehicle
msx3_beh_veh <- msx3_beh_scores %>% filter(recording=="VEH")

wilcox.test(msx3_beh_veh$retrieval ~ msx3_beh_veh$strain)

wilcox.test(msx3_beh_veh$mouthing ~ msx3_beh_veh$strain)

wilcox.test(msx3_beh_veh$corporal ~ msx3_beh_veh$strain)

wilcox.test(msx3_beh_veh$anogenital ~ msx3_beh_veh$strain)

wilcox.test(msx3_beh_veh$nest.building ~ msx3_beh_veh$strain)

#msx3
msx3_beh_msx3 <- msx3_beh_scores %>% filter(recording=="MSX3")

wilcox.test(msx3_beh_msx3$corporal ~ msx3_beh_msx3$strain)

wilcox.test(msx3_beh_msx3$anogenital ~ msx3_beh_msx3$strain)


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
