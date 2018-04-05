#Analysis by individual ANOVAs

#####mom alone and pups separated
maternal_counts <- count_frame_2 %>% filter(recording == "MomAlone" | recording == "PupsSep")

maternal_flat <- maternal_counts %>% filter(categories.allowed=="flat")
lme_maternal_flats <- lme(total.counts ~ strain * recording, 
                    random = ~1| rat.id/recording,
                    data = maternal_flat)
anova.lme(lme_maternal_flats)

maternal_trillc <- maternal_counts %>% filter(categories.allowed=="trill-c")
lme_maternal_trillc <- lme(total.counts ~ strain * recording, 
                          random = ~1| rat.id/recording,
                          data = maternal_trillc)
anova.lme(lme_maternal_trillc)

maternal_short <- maternal_counts %>% filter(categories.allowed=="short")
lme_maternal_short <- lme(total.counts ~ strain * recording, 
                           random = ~1| rat.id/recording,
                           data = maternal_short)
anova.lme(lme_maternal_short)

