##prepare data
data_counts <- read.csv(file.choose(),stringsAsFactors = F)
source("C:/Users/ituncali/Documents/Master's Thesis/USV-Data-Management/src/count_total.R")
library(stringr)
library(dplyr)
allowed.categories <- c("flat", "flat-z", "flat-mz", "short", "short-su", "short-sd",
                        "short-ur", "short-dr", "short-c", "complex", "upward ramp",
                        "downward ramp", "step up", "step down", "multi-step", "multi-step-s",
                        "trill", "trill-c", "trill-f", "inverted-U", "unclear") 
count_list <- by(data = data_counts, INDICES = data_counts$file.name, FUN = function(x) count_total(x, allowed.categories))
count_frame <- do.call(rbind, count_list)
count_frame$file.name <- row.names(count_frame)
row.names(count_frame) <- NULL
count_frame$file.name <- str_replace(count_frame$file.name, pattern = ".[0-9]+$", replacement =  "" )
count_frame <- count_frame %>% group_by(file.name) %>%
  mutate(total.filecounts=sum(total.counts), rel.filecount=total.counts/total.filecounts)


library(reshape2)
file.name.key.unmelted <- read.csv(file.choose(),stringsAsFactors = F)

file.name.key <- melt(data=file.name.key.unmelted,id = c("strain","rat.id"),
                      variable.name="recording",value.name = "file.name")

library(dplyr)
count_frame_2 <- left_join(count_frame,file.name.key)
pup_counts <- count_frame_2 %>% filter(recording == "Fpupiso" | recording == "Mpupiso")
#on 3/9 discovered that WK95 F is missing from this data set b/c no usvs!!
#need to add her
wk95.f <- data.frame(categories.allowed = allowed.categories, 
                     total.counts = c(rep(0,21)),
                     file.name = c(rep("T0000305",21)),
                     total.filecounts = c(rep(0,21)),
                     rel.filecount = c(rep(0,21)),
                     strain = c(rep("WK",21)),
                     rat.id = c(rep("WK95", 21)),
                     recording = c(rep("Fpupiso",21)))

pup_counts <- rbind.data.frame(pup_counts, wk95.f)

pup_counts$rat.id.fix <- as.factor(ifelse(pup_counts$recording=="Fpupiso",paste(pup_counts$rat.id,"F"),paste(pup_counts$rat.id,"M")))
pup_counts$strain <- as.factor(pup_counts$strain)

pup_counts$t.counts <- mapply(function (x) ifelse(x > 0, log(x),NA), pup_counts$total.counts)

shapiro.test(pup_counts[pup_counts$categories.allowed=="flat" & pup_counts$strain=="SD",]$t.counts)
shapiro.test(pup_counts[pup_counts$categories.allowed=="flat" & pup_counts$strain=="WK",]$t.counts)
shapiro.test(pup_counts[pup_counts$categories.allowed=="flat-z" & pup_counts$strain=="SD",]$t.counts)
shapiro.test(pup_counts[pup_counts$categories.allowed=="flat-z" & pup_counts$strain=="WK",]$t.counts)



#kruskal-wallis
kruskal.test(total.counts ~ strain,
             data = pup_counts)

kruskal.test(total.counts ~ categories.allowed,
             data = pup_counts)

pairwise.wilcox.test(pup_counts$total.counts, pup_counts$categories.allowed,
                     p.adjust.method = "BH")


##see which distribution best fits data
library(car)
library(MASS)
pup_counts$total.counts.t <- pup_counts$total.counts + 1

qqp(pup_counts$total.counts.t, "norm")

qqp(pup_counts$total.counts.t, "lnorm")

nbinom <- fitdistr(pup_counts$total.counts.t, "Negative Binomial")
qqp(pup_counts$total.counts.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(pup_counts$total.counts.t, "Poisson")
qqp(pup_counts$total.counts.t, "pois", poisson$estimate)

gamma <- fitdistr(pup_counts$total.counts.t, "gamma")
qqp(pup_counts$total.counts.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#none of these fit. let's try pql

PQL <- glmmPQL(total.counts ~ strain * recording * categories.allowed, ~1 | rat.id.fix, family = gaussian(link = "identity"),
               data = pup_counts, verbose = FALSE)

#anova for pup call counts with label as repeated-measures

pup.aov <- aov(total.counts ~ (strain*recording*categories.allowed) + Error(rat.id.fix/categories.allowed), data = pup_counts)
#must differentiate between M and F rat.ids!!! ^^
#in order to do post-hoc comparisons of repeated measures, must use lme
library(nlme)
#use compound symmetric correlation structure
lme_pups <- lme(total.counts ~ strain * recording * categories.allowed, 
                random = ~1| rat.id.fix,
                correlation = corAR1(form = ~1|rat.id.fix),
                data = pup_counts)
library(lsmeans)
lsmeans(lme_pups, pairwise ~ categories.allowed, adjust = "Tukey")



library(multcomp)
anova(lme_pups)
K <- diag(length(coef(lme_pups)))[-1,]
rownames(K) <- names(coef(lme_pups))[-1]
summary(glht(lme_pups,linfct=K))

#anova for pup call frequencies JUST FLAT, FLAT-Z AND SHORT with label as repeated-measures
data_freqs <- read.csv(file.choose(), stringsAsFactors = F)
m.freq <- apply(data_freqs[,c(1:50)],1,function(y) mean(y))
data_freqs <- mutate(data_freqs, m.freq = m.freq)
pup_freqs <- data_freqs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  group_by(strain, label, recording, file.name, rat.id) %>% 
  summarise(mean.freq = mean(m.freq), sem = sd(m.freq)/sqrt(length(m.freq)),
            count = length(m.freq)) %>%
  filter(label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short")

pup_freqs$rat.id.fix <- ifelse(pup_freqs$recording=="Fpupiso",paste(pup_freqs$rat.id,"F"),paste(pup_freqs$rat.id,"M"))
pup_freqs$strain <- as.factor(pup_freqs$strain)

##check normality
pup_freqs_sh <- data_freqs %>% filter(recording == "Mpupiso" | recording == "Fpupiso")
shapiro.test(pup_freqs_sh[pup_freqs_sh$strain == "WK",]$m.freq)

#not normal....

#lme

lme_pup_freq <- lme(mean.freq ~ strain * recording * label, 
                    random = ~1|rat.id.fix, 
                    correlation = corCompSymm(form = ~1|rat.id.fix),
                    data = pup_freqs)
anova(lme_pup_freq)

K <- diag(length(coef(lme_pup_freq)))[-1,]
rownames(K) <- names(coef(lme_pup_freq))[-1]
summary(glht(lme_pup_freq,linfct=mcp(Material = "Tukey")))



  
#anova for pup call durations with label as repeated-measures
data_durs <- read.csv(file.choose(),stringsAsFactors = F)


pup_durs <- data_durs %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% group_by(strain, label,recording,file.name,rat.id) %>%
summarise(mean.dur = mean(duration), SEM = sd(duration)/sqrt(length(duration)),
count = length(duration)) %>%
filter(label == "flat" | label == "flat-z" | label == "flat-mz" | label == "short")
pup_durs$rat.id.fix <- ifelse(pup_durs$recording=="Fpupiso",paste(pup_durs$rat.id,"F"),paste(pup_durs$rat.id,"M"))


lme_pup_dur <- lme(mean.dur ~ strain * recording * label, 
                   random = ~label|rat.id.fix, 
                   data = pup_durs,
                   correlation = corCompSymm(form = ~label|rat.id.fix)
                   )

#start times

pup_starts <- data_counts %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  group_by(strain, recording, file.name, rat.id) %>% 
  summarise(lat.to.start = min(start.time), lat.to.stop = max(start.time))
pup_starts$rat.id.fix <- ifelse(pup_starts$recording=="Fpupiso",paste(pup_starts$rat.id,"F"),paste(pup_starts$rat.id,"M"))
to.join <- pup_counts %>% group_by(file.name) %>% summarise(total.count = unique(total.filecounts))
pup_starts <- left_join(pup_starts,to.join)
pup_starts <- mutate(pup_starts, rate = total.count/(lat.to.stop-lat.to.start))

lme_pup_starts <- lme(lat.to.start ~ strain * recording, 
                      random = ~1|rat.id.fix,
                      data = pup_starts)
anova.lme(lme_pup_starts)

lme_pup_stops <- lme(lat.to.stop ~ strain * recording, 
                      random = ~1|rat.id.fix,
                      data = pup_starts)
anova.lme(lme_pup_stops)

lme_pup_rates <- pup_starts %>% filter(rate != Inf) %>%
  with(lme(rate ~ strain * recording, 
                     random = ~1|rat.id.fix))
anova.lme(lme_pup_rates)

#by 1 min bins
bins.pups <- cut(start_data[start_data$recording == "Mpupiso" | start_data$recording == "Fpupiso",]$start.time,5,include.lowest=T, labels = c("1","2","3","4","5"))

bin.data.pup <- start_data %>% filter(recording == "Mpupiso" | recording == "Fpupiso") %>% 
  mutate(bin = bins.pups) %>% 
  group_by(bin,strain,file.name, rat.id, recording) %>% 
  summarise(count = sum(no.counts))
bin.data.pup$rat.id.fix <- ifelse(bin.data.pup$recording=="Fpupiso",paste(bin.data.pup$rat.id,"F"),paste(bin.data.pup$rat.id,"M"))

bin.pup.lme <- lme(count ~ strain * bin, random = ~1|file.name, data = bin.data.pup)
anova.lme(bin.pup.lme)

###pup weights and temps
#see if weights differ by strain
pup.w.lme.data <- pup_durs_w %>% group_by(strain, rat.id.fix, recording) %>% summarise(BW = unique(BW))
wk95.f <- data.frame(strain = "WK", rat.id.fix = "WK95 F", recording = "Fpupiso", BW = 8.32)
pup.w.lme.data <- rbind.data.frame(pup.w.lme.data, wk95.f)

pup.w.lme <- lme(BW ~ strain * recording, random = ~1|rat.id.fix, data = pup.w.lme.data)
anova.lme(pup.w.lme)

pup.w <- read.csv(file.choose(),stringsAsFactors = F)
pup_counts_w <- left_join(pup_counts, pup.w)
pup_counts_w <- pup_counts_w %>% filter(!is.na(BW))

pup_counts_w  %>% group_by(strain, rat.id.fix) %>%
  summarise(total.an = sum(total.counts), BW = unique(BW)) %>%
  ggplot(aes(x = BW, y = total.an)) + 
  geom_point(aes(x = BW, y = total.an, colour = strain))

lme.w.count.data <- pup_counts_w %>% group_by(strain, rat.id.fix) %>%
  summarise(total.an = sum(total.counts), BW = unique(BW))

lme.pup.count.w <- lme(total.an ~ strain * BW, random = ~1|rat.id.fix, data = lme.w.count.data)
anova.lme(lme.pup.count.w)

pup_freqs_w <- left_join(pup_freqs, pup.w)
pup_freqs_w <- pup_freqs_w %>% filter(!is.na(BW))

pup_freqs_w  %>% group_by(strain, rat.id.fix) %>%
  summarise(mean.freq.an = mean(mean.freq), sd.freq.an = sd(mean.freq)/sqrt(length(mean.freq)), BW = unique(BW)) %>%
  ggplot(aes(x = BW, y = mean.freq.an)) + 
  geom_point(aes(x = BW, y = mean.freq.an, colour = strain))

lme.w.freq.data <- data_freqs %>% filter(recording == "Fpupiso" | recording == "Mpupiso")
lme.w.freq.data$rat.id.fix <- as.factor(ifelse(lme.w.freq.data$recording=="Fpupiso",paste(lme.w.freq.data$rat.id,"F"),paste(lme.w.freq.data$rat.id,"M")))
lme.w.freq.data <- left_join(lme.w.freq.data, pup.w)
lme.w.freq.data <- lme.w.freq.data %>% filter(!is.na(BW))

lme.pup.freq.w <- lme(m.freq ~ strain * BW, random = ~1|rat.id.fix, data = lme.w.freq.data)
anova.lme(lme.pup.freq.w)

library(ggpubr)
ggscatter(lme.w.freq.data, x = "BW", y = "m.freq", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "pup weight (g)", ylab = "mean frequency (Hz)")
cor.test(lme.w.freq.data$BW, lme.w.freq.data$m.freq, method = "pearson")

pup_durs_w <- left_join(pup_durs, pup.w)
pup_durs_w <- pup_durs_w %>% filter(!is.na(BW))

pup_durs_w  %>% group_by(strain, rat.id.fix, label) %>%
  summarise(mean.dur.an = mean(mean.dur), sd.freq.an = sd(mean.dur)/sqrt(length(mean.dur)), BW = unique(BW)) %>%
  ggplot(aes(x = BW, y = mean.dur.an)) + 
  geom_point(aes(x = BW, y = mean.dur.an, colour = strain))

lme.w.dur.data <- data_durs %>% filter(recording == "Fpupiso" | recording == "Mpupiso")
lme.w.dur.data$rat.id.fix <- as.factor(ifelse(lme.w.dur.data$recording=="Fpupiso",paste(lme.w.dur.data$rat.id,"F"),paste(lme.w.dur.data$rat.id,"M")))
lme.w.dur.data <- left_join(lme.w.dur.data, pup.w)
lme.w.dur.data <- lme.w.dur.data %>% filter(!is.na(BW))

lme.pup.dur.w <- lme(duration ~ strain * BW, random = ~1|rat.id.fix, data = lme.w.dur.data)
anova.lme(lme.pup.dur.w)

library(ggpubr)
ggscatter(lme.w.dur.data, x = "BW", y = "duration", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "pup weight (g)", ylab = "duration (s)")

#pup litter weights
litter.w <- read.csv(file.choose(),stringsAsFactors = F)

litter.w.lme <- lme(bw.pup ~ strain, random = ~1|rat.id, data = litter.w)
anova.lme(litter.w.lme)


##pup total calling times
data_counts <- left_join(data_counts, file.name.key)
pup_start_counts <- filter(data_counts, recording == "Fpupiso" | recording == "Mpupiso")
pup_start_counts$rat.id.fix <- as.factor(ifelse(pup_start_counts$recording=="Fpupiso",paste(pup_start_counts$rat.id,"F"),paste(pup_start_counts$rat.id,"M")))
pup_start_counts <- pup_start_counts %>% group_by(strain, rat.id.fix, recording) %>%
  summarise(tot.call.time = max(start.time) - min(start.time), lat.to.start.call = min(start.time))

lme.tot.call <- lme(tot.call.time ~ strain * recording,
                    random = ~1|rat.id.fix,
                    data = pup_start_counts)
anova.lme(lme.tot.call)

