##individual anovas

flat.3 <- count_frame_3 %>% filter(categories.allowed=="flat")

flat.lme <- lme(total.counts~strain*recording,
                random=~1|rat.id,
                data=flat.3)
anova.lme(flat.lme)
flat.sum <- summary(lsmeans(flat.lme, pairwise~strain*recording, adjust="Dunnet"))[["contrasts"]]
View(flat.sum)


short.3 <- count_frame_3 %>% filter(categories.allowed=="short")

short.lme <- lme(total.counts~strain*recording,
                random=~1|rat.id,
                data=short.3)
anova.lme(short.lme)
short.sum <- summary(lsmeans(short.lme, pairwise~strain, adjust="Dunnet"))[["contrasts"]]
View(short.sum[short.sum$p.value<0.05,])



trill.3 <- count_frame_3 %>% filter(categories.allowed=="trill")

trill.lme <- lme(total.counts~strain*recording,
                random=~1|rat.id,
                data=trill.3)
anova.lme(trill.lme)
trill.sum <- summary(lsmeans(trill.lme, pairwise~strain*recording, adjust="Dunnet"))[["contrasts"]]
View(trill.sum[trill.sum$p.value<0.05,])

trillc.3 <- count_frame_3 %>% filter(categories.allowed=="trill-c")

trillc.lme <- lme(total.counts~strain*recording,
                 random=~1|rat.id,
                 data=trillc.3)
anova.lme(trillc.lme)
trillc.sum <- summary(lsmeans(trillc.lme, pairwise~strain*recording, adjust="Dunnet"))[["contrasts"]]
View(trillc.sum[trillc.sum$p.value<0.05,])


flatz.3 <- count_frame_3 %>% filter(categories.allowed=="flat-z")

flatz.lme <- lme(total.counts~strain*recording,
                  random=~1|rat.id,
                  data=flatz.3)
anova.lme(flatz.lme)



###mom usvs
mom.flat <- mom_counts %>% filter(label=="flat")
m.flat.lme <- lme(usv.count~strain*recording,
                  random=~1|rat.id,
                  data=mom.flat)
anova.lme(m.flat.lme)

mom.trill <- mom_counts %>% filter(label=="trill")
m.trill.lme <- lme(usv.count~strain*recording,
                  random=~1|rat.id,
                  data=mom.trill)
anova.lme(m.trill.lme)


mom.trillc <- mom_counts %>% filter(label=="trill-c")
m.trillc.lme <- lme(usv.count~strain*recording,
                  random=~1|rat.id,
                  data=mom.trillc)
anova.lme(m.trillc.lme)








