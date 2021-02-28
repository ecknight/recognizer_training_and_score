#title: Analysis and figures for Knight and Bayne 2018. Classification threshold and training data affect the quality and utility of focal species data processed with automated audio-recognition software. Bioacoustics 10.1080/09524622.2018.1503971.
#author: Elly C. Knight
#date: Oct. 1, 2018

library(stringi)
library(tidyverse)
library(readxl)
library(aspace)
library(AICcmodavg)
library(ggforce)
library(colorspace)
library(mrds)
library(Distance)
library(lme4)
library(grDevices)
library(gridExtra)

#Proprietary library of functions from Zuur & Ieno
source("/Users/ellyknight/Documents/R/Scripts/HighstatLibV10.R")

options(scipen=999)

my.theme <- theme_classic() +
  theme(text=element_text(size=14, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

mround <- function(x,base){ 
  base*ceiling(x/base) 
} 

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#A. WRANGLING-----

#1. Read in validated data----

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Analysis/Score")

r1 <- read.table("Clips7_CONIPeent0m_0_20_validated.txt", header=FALSE, sep="\t")
r2 <- read.table("Clips7_CONIPeent120m_0_20_validated.txt", header=FALSE, sep="\t")
r3 <- read.table("Clips7_CONIPeent240m_0_20_validated.txt", header=FALSE, sep="\t")
r4 <- read.table("Clips7_CONIPeentMixed_0_20_validated.txt", header=FALSE, sep="\t")

r5 <- rbind(r1, r2, r3, r4)
colnames(r5) <- c("filepath", "start", "duration", "level", "quality", "score", "recognizer", "comments")

r5$file <- stri_sub(gsub("\\", "/", r5$filepath, fixed=TRUE), from=28, length = 60)
levels(r5$recognizer) <- c("near training data", "midrange training data", "far training data", "mixed training data")

r6 <- r5 %>% 
  separate(file, c("clip", "s", "aruID", "date", "time"), sep="_", remove=FALSE) %>% 
  separate(aruID, c("project", "transect", "aru"), remove=FALSE) %>% 
  mutate(transect=as.numeric(transect),
         aru=as.numeric(aru),
         clip=as.numeric(clip))

#2. Filter out stuff----
rout <- r6 %>% 
  filter(comments != "y") %>% 
  dplyr::select(clip, aru) %>% 
  unique()

r7 <- r6 %>% 
  anti_join(rout) %>% 
  filter(aru !=0,
         aru <=500)

#3. Add the distances----
bearing <- read.csv("EDRTransectBearings.csv")
dir.correct <- read.csv("DirectionConversions.csv")
det <- read_excel("qry_observations.xlsx")

colnames(det) <- c("clip", "horizontal", "cardinal", "vertical", "boom", "transect", "zero", "UTME", "UTMN", "use", "DateTime")

r8 <- r7 %>% 
  left_join(det, by=c("clip", "transect")) %>% 
  left_join(bearing, by="transect") %>% 
  full_join(dir.correct, by="cardinal") %>% 
  mutate(distance = sqrt(aru^2 + horizontal^2 - 2*aru*horizontal*cos_d(abs(bearing-direction))),
         distance=ifelse(horizontal==0, aru, distance))

#4. Remove training clips----

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Processing/Score/Annotations/2.0/SongScopeNotes0m/Use")
files0 <- as.data.frame(list.files(pattern=".ssn"))
files0$recognizer <- "near training data"

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Processing/Score/Annotations/2.0/SongScopeNotes120m/Use")
files120 <- as.data.frame(list.files(pattern=".ssn"))
files120$recognizer <- "midrange training data"

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Processing/Score/Annotations/2.0/SongScopeNotes240m/Use")
files240 <- as.data.frame(list.files(pattern=".ssn"))
files240$recognizer <- "far training data"

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Processing/Score/Annotations/2.0/SongScopeNotesMixed/Use")
filesmixed <- as.data.frame(list.files(pattern=".ssn"))
filesmixed$recognizer <- "mixed training data"

files <- rbind(files0, files120, files240, filesmixed)
colnames(files) <- c("file.join", "recognizer")
files$file.join <- as.character(stri_sub(files$file.join, -42, -5))

r9 <- r8 %>% 
  mutate(file.join = stri_sub(file, 1, 42),
         file.join = stri_sub(file.join, -42, -5))

r10 <- r9[!r9$file.join %in% files$file.join,]


#5. Remove clips from cleaning----

r <- r10 %>% 
  filter(clip != 3821,
         clip != 3657, 
         clip != 3626,
         clip != 3532,
         clip != 3480,
         clip != 3297 | aru != 500,
         clip != 3324 | aru != 400, 
         clip != 3814 | aru != 500)

#6. Final list of clips----

clips <- r %>% 
  dplyr::select(clip) %>% 
  unique()
colnames(clips) <- "clip"
nrow(clips)
nrow(clips)*11

#B. DISTANCE MODELLING----

#1. Zuur data vis & cleaning----

MyVar <- c("score", "level", "distance")

#look for outliers
Mydotplot(r[,MyVar])

#collinearity
Mypairs(r[, MyVar])
par(mfrow = c(1, 2))

#try log transforming distance
r$distance.log <- log(r$distance+1)
r$score.log <- log(r$score)
r$level.log <- log(r$level)

#2. Model----

#split into recognizers
r.0 <- subset(r, recognizer=="near training data")
r.120 <- subset(r, recognizer=="midrange training data")
r.240 <- subset(r, recognizer=="far training data")
r.mix <- subset(r, recognizer=="mixed training data")

#Visualize
ggplot(r, aes(x=distance, y=score, colour=as.factor(transect))) +
  geom_point() +
  geom_smooth(colour="grey60", method="gam") +
  facet_wrap(~recognizer)

#2a. 0 m----

lmsd.0 <- list()
lmsd.0[[1]] <- lm(score ~ 1, data=r.0)
lmsd.0[[2]] <- lm(score ~ distance, data=r.0)
lmsd.0[[3]] <- lm(score ~ distance + I(distance^2), data=r.0)
lmsd.0[[4]] <- lm(score ~ distance + I(distance^2) + I(distance^3), data=r.0)
aictab(lmsd.0)
#choose lmsd.0[[3]]

#homogeneity
lm.E <- resid(lmsd.0[[3]])
lm.F <- fitted(lmsd.0[[3]])
par(mfrow=c(1,1))
plot(x = lm.F,
     y = lm.E,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

#influence
par(mfrow = c(1, 1))
plot(cooks.distance(lmsd.0[[3]]), type = "h", ylim = c(0, 1))
abline(h = 1)

#Normality
hist(lm.E, breaks = 15)

#independence
plot(x = r.0$distance,
     y = lm.E)
abline(h = 0, lty = 2)

#2b. 120 m----
lmsd.120 <- list()
lmsd.120[[1]] <- lm(score ~ 1, data=r.120)
lmsd.120[[2]] <- lm(score ~ distance, data=r.120)
lmsd.120[[3]] <- lm(score ~ distance + I(distance^2), data=r.120)
lmsd.120[[4]] <- lm(score ~ distance + I(distance^2) + I(distance^3), data=r.120)
aictab(lmsd.120)
#choose lmsd.0[[4]]

#homogeneity
lm.E <- resid(lmsd.0[[4]])
lm.F <- fitted(lmsd.0[[4]])
par(mfrow=c(1,1))
plot(x = lm.F,
     y = lm.E,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

#influence
par(mfrow = c(1, 1))
plot(cooks.distance(lmsd.0[[4]]), type = "h", ylim = c(0, 1))
abline(h = 1)

#Normality
hist(lm.E, breaks = 15)

#independence
plot(x = r.120$distance,
     y = lm.E)
abline(h = 0, lty = 2)

#2c. 240 m----

lmsd.240 <- list()
lmsd.240[[1]] <- lm(score ~ 1, data=r.240)
lmsd.240[[2]] <- lm(score ~ distance, data=r.240)
lmsd.240[[3]] <- lm(score ~ distance + I(distance^2), data=r.240)
lmsd.240[[4]] <- lm(score ~ distance + I(distance^2) + I(distance^3), data=r.240)
aictab(lmsd.240)
#choose lmsd.240[[4]]

#homogeneity
lm.E <- resid(lmsd.240[[4]])
lm.F <- fitted(lmsd.240[[4]])
par(mfrow=c(1,1))
plot(x = lm.F,
     y = lm.E,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

#influence
par(mfrow = c(1, 1))
plot(cooks.distance(lmsd.240[[4]]), type = "h", ylim = c(0, 1))
abline(h = 1)

#Normality
hist(lm.E, breaks = 15)

#independence
plot(x = r.240$distance,
     y = lm.E)
abline(h = 0, lty = 2)

#2d. Mixed----

lmsd.mix <- list()
lmsd.mix[[1]] <- lm(score ~ 1, data=r.mix)
lmsd.mix[[2]] <- lm(score ~ distance, data=r.mix)
lmsd.mix[[3]] <- lm(score ~ distance + I(distance^2), data=r.mix)
lmsd.mix[[4]] <- lm(score ~ distance + I(distance^2) + I(distance^3), data=r.mix)
aictab(lmsd.mix)
#choose lmsd.mix[[3]]

#homogeneity
lm.E <- resid(lmsd.mix[[3]])
lm.F <- fitted(lmsd.mix[[3]])
par(mfrow=c(1,1))
plot(x = lm.F,
     y = lm.E,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, v = 0, lty = 2)

#influence
par(mfrow = c(1, 1))
plot(cooks.distance(lmsd.mix[[2]]), type = "h", ylim = c(0, 1))
abline(h = 1)

#Normality
hist(lm.E, breaks = 15)

#independence
plot(x = r.mix$distance,
     y = lm.E)
abline(h = 0, lty = 2)

#3. Plot----

lmsd.0.use <- lm(score ~ distance + I(distance^2), data=r.0)
summary(lmsd.0.use)
lmsd.120.use <- lm(score ~ distance + I(distance^2) + I(distance^3), data=r.120)
summary(lmsd.120.use)
lmsd.240.use <- lm(score ~ distance + I(distance^2) + I(distance^3), data=r.240)
summary(lmsd.240.use)
lmsd.mix.use <- lm(score ~ distance + I(distance^2), data=r.mix)
summary(lmsd.mix.use)

MyData1 <- expand.grid(distance = seq(from = min(r.0$distance),
                                     to = max(r.0$distance),
                                     length = 1000))
P.0 <- predict(lmsd.0.use, newdata = MyData1, se = TRUE)
MyData1$mu    <- P.0$fit                #Fitted values
MyData1$selow <- P.0$fit - 2 * P.0$se.fit #lower bound
MyData1$seup  <- P.0$fit + 2 * P.0$se.fit #lower bound
MyData1$recognizer <- "near training data"

MyData2 <- expand.grid(distance = seq(from = min(r.120$distance),
                                   to = max(r.120$distance),
                                   length = 1000))
P.120 <- predict(lmsd.120.use, newdata = MyData2, se = TRUE)
MyData2$mu    <- P.120$fit                #Fitted values
MyData2$selow <- P.120$fit - 2 * P.120$se.fit #lower bound
MyData2$seup  <- P.120$fit + 2 * P.120$se.fit #lower bound
MyData2$recognizer <- "midrange training data"

MyData3 <- expand.grid(distance = seq(from = min(r.240$distance),
                                   to = max(r.240$distance),
                                   length = 1000))
P.240 <- predict(lmsd.240.use, newdata = MyData3, se = TRUE)
MyData3$mu    <- P.240$fit                #Fitted values
MyData3$selow <- P.240$fit - 2 * P.240$se.fit #lower bound
MyData3$seup  <- P.240$fit + 2 * P.240$se.fit #lower bound
MyData3$recognizer <- "far training data"

MyData4 <- expand.grid(distance = seq(from = min(r.mix$distance),
                                   to = max(r.mix$distance),
                                   length = 1000))
P.mix <- predict(lmsd.mix.use, newdata = MyData4, se = TRUE)
MyData4$mu    <- P.mix$fit                #Fitted values
MyData4$selow <- P.mix$fit - 2 * P.mix$se.fit #lower bound
MyData4$seup  <- P.mix$fit + 2 * P.mix$se.fit #lower bound
MyData4$recognizer <- "mixed training data"

MyDataAllD <- rbind(MyData1, MyData2, MyData3, MyData4)
MyDataAllD$recognizer <- factor(MyDataAllD$recognizer, levels=c("near training data", "midrange training data", "far training data", "mixed training data"))
summary(MyDataAllD$recognizer)

fig3 <- ggplot() +
  geom_point(aes(y=score, x=distance),
             data=r, colour="grey40", alpha=0.3) +
  geom_ribbon(aes(x=distance, ymax=seup, ymin=selow),
              data=MyDataAllD,
              alpha=0.3, fill="grey40") +
  geom_line(aes(y=mu, x=distance),
            data=MyDataAllD, size=1, colour="grey20") +
  facet_wrap(~recognizer) +
  ylab("Recognizer score") +
  xlab("Distance (m)") +
  my.theme +
  theme(legend.position="none")
fig3

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Writing/Bioacoustics_Submission")
ggsave("Fig3.png", fig3, width=6, height=6, dpi=600)

#C. RELATIONSHIP RELIABILITY----

#1. Wrangling----

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Analysis/Score")

covs <- read.csv("qry_covs.csv") %>% 
  dplyr::select(pkey_ObservationID, Temp, Wind, Humidity) %>% 
  rename(clip=pkey_ObservationID) %>% 
  unique()

r.0.cov <- r.0 %>% 
  left_join(covs, by="clip") %>% 
  mutate(transect=as.factor(transect))

#2. VIF----

library(usdm)

covs.weather <- r.0.cov %>% 
  dplyr::select(Temp, Wind, Humidity) #extract just the habitat covariates
usdm::vifstep(covs.weather, th=10) #run stepwise VIF removal

#3. Model----

#library(nlme)
#lme.0 <- lme(score~distance + I(distance^2) + Temp + Wind + Humidity, data=r.0.cov, random= ~1|transect/clip)

library(nlme)
lme.0.cov <- lme(score ~ distance + I(distance^2) + Temp + Wind + Humidity, data=r.0.cov, random= ~1|transect)
summary(lme.0.cov)
lm.0.cov <- lm(score ~ distance + I(distance^2) + Temp + Wind + Humidity, data=r.0.cov)
summary(lm.0.cov)

library(lme4)
lme.0 <- list()
lme.0[[1]] <- lmer(score ~ distance + I(distance^2) + Temp + Wind + Humidity + (1|transect/clip), data=r.0.cov, REML=FALSE)
lme.0[[2]] <- lmer(score ~ distance + I(distance^2) + Temp + Wind + (1|transect/clip), data=r.0.cov, REML=FALSE)
lme.0[[3]] <- lmer(score ~ distance + I(distance^2) + Temp + Humidity + (1|transect/clip), data=r.0.cov, REML=FALSE)
lme.0[[4]] <- lmer(score ~ distance + I(distance^2) + Wind + Humidity + (1|transect/clip), data=r.0.cov, REML=FALSE)
lme.0[[5]] <- lmer(score ~ distance + I(distance^2) + Temp + (1|transect/clip), data=r.0.cov, REML=FALSE)
lme.0[[6]] <- lmer(score ~ distance + I(distance^2) + Wind + (1|transect/clip), data=r.0.cov, REML=FALSE)
lme.0[[7]] <- lmer(score ~ distance + I(distance^2) + Humidity + (1|transect/clip), data=r.0.cov, REML=FALSE)
lme.0[[8]] <- lmer(score ~ distance + I(distance^2) + (1|transect/clip), data=r.0.cov, REML=FALSE)
aictab(lme.0)
#choose null model

#4. ICC----
lme.0.use <- lmer(score ~ distance + I(distance^2) + (1|transect/clip), data=r.0.cov, REML=TRUE)
summary(lme.0.use)
theta.clip <- attributes(VarCorr(lme.0.use)$clip)$stddev
theta.transect <- attributes(VarCorr(lme.0.use)$transect)$stddev
sigma <- lme.0.use@devcomp$cmp[10]
icc.clip <- (theta.clip^2 + theta.transect^2)/(theta.clip^2 + theta.transect^2 + sigma^2)
icc.clip
icc.transect <- (theta.transect^2)/(theta.clip^2 + theta.transect^2 + sigma^2)
icc.transect

#5. Plot----
#https://stats.stackexchange.com/questions/98958/plots-to-illustrate-results-of-linear-mixed-effect-model

lme.0.use <- lmer(score ~ distance + I(distance^2) + (1|transect/clip), data=r.0.cov, REML=TRUE)

fit <- predict(lme.0.use, r.0.cov)
MyData5 <- cbind(r.0.cov, fit)
MyData5$transect <- factor(MyData5$transect, labels=c("Individual 1", "Individual 2", "Individual 3", "Individual 4", "Individual 5"))

colours <- grey.colors(5, start = 0.1, end = 0.8, gamma = 2.2, alpha = NULL)

fig4 <- ggplot(MyData5) +
  geom_point(aes(x=distance, y=score), alpha = 0.3, colour="grey40") +
  geom_line(aes(x=distance, y=fit, group=clip), alpha = 1.0, colour="grey20") +
  facet_wrap(~transect, ncol=5) +
  ylab("Recognizer score") +
  xlab("Distance (m)") +
  my.theme
fig4


setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/Score/Writing/Bioacoustics_Revision")
ggsave("Fig4.png", fig4, width=15, height=3, dpi=600)

#D. PROBABILITY OF DETECTION----

#1. Wrangle----

#Read in processed data
setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Analysis/Score")
r1 <- read.table("Clips7_CONIPeent0m_0_20_validated.txt", header=FALSE, sep="\t")
r2 <- read.table("Clips7_CONIPeent120m_0_20_validated.txt", header=FALSE, sep="\t")
r3 <- read.table("Clips7_CONIPeent240m_0_20_validated.txt", header=FALSE, sep="\t")
r4 <- read.table("Clips7_CONIPeentMixed_0_20_validated.txt", header=FALSE, sep="\t")

r5 <- rbind(r1, r2, r3, r4)
colnames(r5) <- c("filepath", "start", "duration", "level", "quality", "score", "recognizer", "comments")

r5$file <- stri_sub(gsub("\\", "/", r5$filepath, fixed=TRUE), from=28, length = 60)
levels(r5$recognizer) <- c("near training data", "midrange training data", "far training data", "mixed training data")

r6 <- r5 %>% 
  separate(file, c("clip", "s", "aruID", "date", "time"), sep="_", remove=FALSE) %>% 
  separate(aruID, c("project", "transect", "aru"), remove=FALSE) %>% 
  filter(clip != 3297 | aru != 500,
         clip != 3324 | aru != 400, 
         clip != 3814 | aru != 500)

rout <- r6 %>% 
  filter(comments != "y") %>% 
  dplyr::select(clip, aru) %>% 
  unique()

r7 <- r6 %>% 
  anti_join(rout) %>% 
  filter(aru !=0,
         aru <=500)

#Add in undetected clips
files <- read.csv("ClipList.csv")
files.exp <- expand.grid(file=files$file, recognizer=c("near training data", "midrange training data", "far training data", "mixed training data")
)
r8 <- left_join(files.exp, r7, by=c("file", "recognizer")) %>% 
  dplyr::select(file, level, score, comments, recognizer) %>% 
  separate(file, c("clip", "s", "aruID", "date", "time"), sep="_", remove=FALSE) %>% 
  separate(aruID, c("project", "transect", "aru"), remove=FALSE) %>% 
  mutate(transect=as.numeric(transect),
         aru=as.numeric(aru),
         clip=as.numeric(clip))

#Add in distances
bearing <- read.csv("EDRTransectBearings.csv")
dir.correct <- read.csv("DirectionConversions.csv")
det <- read_excel("qry_observations.xlsx")

colnames(det) <- c("clip", "horizontal", "cardinal", "vertical", "boom", "transect", "zero", "UTME", "UTMN", "use", "DateTime")

r9 <- r8 %>% 
  left_join(det, by=c("clip", "transect")) %>% 
  left_join(bearing, by="transect") %>% 
  full_join(dir.correct, by="cardinal") %>% 
  mutate(distance = sqrt(aru^2 + horizontal^2 - 2*aru*horizontal*cos_d(abs(bearing-direction))),
         distance=ifelse(horizontal==0, aru, distance)) %>% 
  filter(aru !=0) %>% #take out the 0 m clips
  filter(aru <= 500) %>% #take out the clips about 500 m
  mutate(hit=ifelse(is.na(score), 0, 1),
         p=(distance^2)/((max(distance))^2),
         distancebin=mround(distance, 50))

#Filter by final clip list
r9clips <- r9 %>% 
  dplyr::select(clip) %>% 
  unique

r10 <- r9 %>% 
  inner_join(clips, by="clip") %>% 
  mutate(distance2=-distance^2)

table(r9$clip, r9$aru)

write.csv(r10, "FinalClipListAndRecognizerResults.csv", row.names = FALSE)

#2. Descriptive stats----

r.10.0 <- subset(r10, recognizer=="near training data")
summary(as.factor(r.10.0$hit))
r.10.120 <- subset(r10, recognizer=="midrange training data")
summary(as.factor(r.10.120$hit))
r.10.240 <- subset(r10, recognizer=="far training data")
summary(as.factor(r.10.240$hit))
r.10.mix <- subset(r10, recognizer=="mixed training data")
summary(as.factor(r.10.mix$hit))

#3. Model----

#3a. No score threshold----
p.0 <- glm(hit~distance2, family=binomial("cloglog"), data=r.10.0)
#logit2prob(confint(p.0))
edr.0 <- sqrt(1/coef(p.0)[2])
edr.0
p.120 <- glm(hit~distance2, family=binomial("cloglog"), data=r.10.120)
#logit2prob(confint(p.120))
edr.120 <- sqrt(1/coef(p.120)[2])
edr.120
p.240 <- glm(hit~distance2, family=binomial("cloglog"), data=r.10.240)
#logit2prob(confint(p.240))
edr.240 <- sqrt(1/coef(p.240)[2])
edr.240
p.mix <- glm(hit~distance2, family=binomial("cloglog"), data=r.10.mix)
#logit2prob(confint(p.mix))
edr.mix <- sqrt(1/coef(p.mix)[2])
edr.mix

MyData6 <- expand.grid(distance2 = seq(from = min(r.10.0$distance2),
                                      to = max(r.10.0$distance2),
                                      length = 1000))
p.p.0 <- predict(p.0, newdata = MyData6, se = TRUE, type="response")
MyData6$mu    <- p.p.0$fit                   #Fitted values
MyData6$selow <- p.p.0$fit - 2 * p.p.0$se.fit #lower bound
MyData6$seup  <- p.p.0$fit + 2 * p.p.0$se.fit #lower bound
MyData6$recognizer <- "near training data"
MyData6$distance <- sqrt(-MyData6$distance2)

MyData7 <- expand.grid(distance2 = seq(from = min(r.10.120$distance2),
                                      to = max(r.10.120$distance2),
                                      length = 1000))
p.p.120 <- predict(p.120, newdata = MyData7, se = TRUE, type="response")
MyData7$mu    <- p.p.120$fit                #Fitted values
MyData7$selow <- p.p.120$fit - 2 * p.p.120$se.fit #lower bound
MyData7$seup  <- p.p.120$fit + 2 * p.p.120$se.fit #lower bound
MyData7$recognizer <- "midrange training data"
MyData7$distance <- sqrt(-MyData7$distance2)

MyData8 <- expand.grid(distance2 = seq(from = min(r.10.240$distance2),
                                      to = max(r.10.240$distance2),
                                      length = 1000))
p.p.240 <- predict(p.240, newdata = MyData8, se = TRUE, type="response")
MyData8$mu    <- p.p.240$fit                #Fitted values
MyData8$selow <- p.p.240$fit - 2 * p.p.240$se.fit #lower bound
MyData8$seup  <- p.p.240$fit + 2 * p.p.240$se.fit #lower bound
MyData8$recognizer <- "far training data"
MyData8$distance <- sqrt(-MyData8$distance2)

MyData9 <- expand.grid(distance2 = seq(from = min(r.10.mix$distance2),
                                      to = max(r.10.mix$distance2),
                                      length = 1000))
p.p.mix <- predict(p.mix, newdata = MyData9, se = TRUE, type="response")
MyData9$mu    <- p.p.mix$fit                #Fitted values
MyData9$selow <- p.p.mix$fit - 2 * p.p.mix$se.fit #lower bound
MyData9$seup  <- p.p.mix$fit + 2 * p.p.mix$se.fit #lower bound
MyData9$recognizer <- "mixed training data"
MyData9$distance <- sqrt(-MyData9$distance2)

MyData <- rbind(MyData6, MyData7, MyData8, MyData9)
MyData$recognizer <- factor(MyData$recognizer, levels=c("near training data", "midrange training data", "far training data", "mixed training data"))
summary(-MyData$recognizer)

colours <- heat_hcl(4, h=c(50, 320), l=c(80,50), c=c(100,50), power=1)

r.10.pres <- subset(r.10.0, hit==1)

ggplot() +
  geom_ribbon(aes(x=distance, ymax=seup, ymin=selow, group=factor(recognizer)),
              data=MyData,
              alpha=0.5, fill="grey40", colour=NA) +
  geom_line(aes(y=mu, x=distance, colour=factor(recognizer)),
            data=MyData, size=1) +
  ylab("Probability of detection") +
  xlab("Distance (m)") +
  scale_y_continuous(limits=c(0,1))+
  my.theme +
  scale_colour_manual(name="Recognizer", values=colours)

ggplot(data=r.10.pres) +
  geom_histogram(aes(distance), binwidth=30) +
  xlim(c(30,500))

#3b. Score threshold----

score <- c(40, 45, 50, 55, 60, 65, 70)

edr.all <- data.frame()
det.all <- data.frame()
for (i in 1:length(score)) {
  
  score.i <- score[i]
  
  r <- r10 %>% 
    mutate(hit.score = ifelse(score < score.i, 0, hit),
           hit.score = ifelse(is.na(hit.score), 0, hit.score))
  
  #split into recognizers
  r.0 <- subset(r, recognizer=="near training data")
  r.120 <- subset(r, recognizer=="midrange training data")
  r.240 <- subset(r, recognizer=="far training data")
  r.mix <- subset(r, recognizer=="mixed training data")
  
  p.0 <- glm(hit.score~distance2, binomial("cloglog"), data=r.0)
  edr.0 <- sqrt(1/coef(p.0)[2])
  p.120 <- glm(hit.score~distance2, binomial("cloglog"), data=r.120)
  edr.120 <- sqrt(1/coef(p.120)[2])
  p.240 <- glm(hit.score~distance2, binomial("cloglog"), data=r.240)
  edr.240 <- sqrt(1/coef(p.240)[2])
  p.mix <- glm(hit.score~distance2, family="binomial", data=r.mix)
  edr.mix <- sqrt(1/coef(p.mix)[2])
  
  MyData6 <- expand.grid(distance2 = seq(from = min(r.0$distance2),
                                        to = max(r.0$distance2),
                                        length = 1000))
  p.p.0 <- predict(p.0, newdata = MyData6, se = TRUE, type="response")
  MyData6$mu    <- p.p.0$fit                #Fitted values
  MyData6$selow <- p.p.0$fit - 2 * p.p.0$se.fit #lower bound
  MyData6$seup  <- p.p.0$fit + 2 * p.p.0$se.fit #lower bound
  MyData6$recognizer <- "near training data"
  MyData6$distance <- sqrt(-MyData6$distance2)
  
  MyData7 <- expand.grid(distance2 = seq(from = min(r.120$distance2),
                                        to = max(r.120$distance2),
                                        length = 1000))
  p.p.120 <- predict(p.120, newdata = MyData7, se = TRUE, type="response")
  MyData7$mu    <- p.p.120$fit                #Fitted values
  MyData7$selow <- p.p.120$fit - 2 * p.p.120$se.fit #lower bound
  MyData7$seup  <- p.p.120$fit + 2 * p.p.120$se.fit #lower bound
  MyData7$recognizer <- "midrange training data"
  MyData7$distance <- sqrt(-MyData7$distance2)
  
  MyData8 <- expand.grid(distance2 = seq(from = min(r.240$distance2),
                                        to = max(r.240$distance2),
                                        length = 1000))
  p.p.240 <- predict(p.240, newdata = MyData8, se = TRUE, type="response")
  MyData8$mu    <- p.p.240$fit                #Fitted values
  MyData8$selow <- p.p.240$fit - 2 * p.p.240$se.fit #lower bound
  MyData8$seup  <- p.p.240$fit + 2 * p.p.240$se.fit #lower bound
  MyData8$recognizer <- "far training data"
  MyData8$distance <- sqrt(-MyData8$distance2)
  
  MyData9 <- expand.grid(distance2 = seq(from = min(r.mix$distance2),
                                        to = max(r.mix$distance2),
                                        length = 1000))
  p.p.mix <- predict(p.mix, newdata = MyData9, se = TRUE, type="response")
  MyData9$mu    <- p.p.mix$fit                #Fitted values
  MyData9$selow <- p.p.mix$fit - 2 * p.p.mix$se.fit #lower bound
  MyData9$seup  <- p.p.mix$fit + 2 * p.p.mix$se.fit #lower bound
  MyData9$recognizer <- "mixed training data"
  MyData9$distance <- sqrt(-MyData9$distance2)
  
  MyData <- rbind(MyData6, MyData7, MyData8, MyData9)
  MyData$score <- score.i
  
  edr <- data.frame(rbind(edr.0, edr.120, edr.240, edr.mix)) %>% 
    dplyr::rename(edr=distance2)
  edr$score <- score.i
  
  det.all <- rbind(det.all, MyData)
  edr.all <- rbind(edr, edr.all)
  
}

#4. Plot----

det.all$recognizer <- factor(det.all$recognizer, levels=c("near training data", "midrange training data", "far training data", "mixed training data"))
summary(det.all$recognizer)

colours <- grey.colors(7, start = 0, end = 0.9, gamma = 1, alpha = NULL)
lines <- c(1:7)

fig6 <- ggplot() +
  geom_ribbon(aes(x=distance, ymax=seup, ymin=selow, group=factor(score)),
              data=det.all,
              alpha=0.3, fill="grey40") +
  geom_line(aes(y=mu, x=distance, colour=factor(score), linetype=factor(score)),
            data=det.all, size=1) +
  facet_wrap(~recognizer) +
  ylab("Probability of detection") +
  xlab("Distance (m)") +
  my.theme +
  scale_colour_manual(name="Score\nthreshold", values=colours) +
  scale_linetype_manual(name="Score\nthreshold", values=lines)
fig6

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Writing/Bioacoustics_Submission")
ggsave("Fig6.png", fig6, width=8, height=6, dpi=600)

#E. RECALL----

score <- c(40, 45, 50, 55, 60, 65, 70)

metrics <- data.frame()
for (i in 1:length(score)) {
  
  score.i <- score[i]
  
  r <- r10 %>% 
    mutate(hit.score = ifelse(score < score.i, 0, hit),
           hit.score = ifelse(is.na(hit.score), 0, hit.score))
  
  metrics.i <- r %>% 
    group_by(recognizer) %>% 
    dplyr::summarize(hits=sum(hit.score),
              total=n()) %>% 
    mutate(recall=hits/total,
           score=score.i)
  
  metrics <- rbind(metrics, metrics.i)
}

metrics$recognizer <- factor(metrics$recognizer, levels=c("near training data", "midrange training data", "far training data", "mixed training data"))
summary(metrics$recognizer)

colours <- grey.colors(4, start = 0.1, end = 0.8, gamma = 2.2, alpha = NULL)

fig5 <- ggplot() +
  geom_line(aes(y=recall, x=score, colour=factor(recognizer)),
            data=metrics, size=1) +
  geom_point(aes(y=recall, x=score, colour=factor(recognizer)),
             data=metrics, size=2) +
  ylab("Recall") +
  xlab("Score threshold") +
  ylim(c(0,0.8)) +
  my.theme +
  scale_colour_manual(name="Recognizer", values=colours)
fig5

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Writing/Bioacoustics_Submission")
ggsave("Fig5.png", fig5, width=5.5, height=3, dpi=600)

#F. TRAINING DATA DESCRIPTION----

files1 <- files %>% 
  separate(file.join, c("clip", "s", "aruID", "date", "time"), sep="_", remove=FALSE) %>% 
  separate(aruID, c("project", "transect", "aru"), remove=FALSE) %>% 
  mutate(transect=as.numeric(transect),
         aru=as.numeric(aru),
         clip=as.numeric(clip))

files2 <- files1 %>% 
  left_join(det, by=c("clip", "transect")) %>% 
  left_join(bearing, by="transect") %>% 
  full_join(dir.correct, by="cardinal") %>% 
  mutate(distance = sqrt(aru^2 + horizontal^2 - 2*aru*horizontal*cos_d(abs(bearing-direction))),
         distance=ifelse(horizontal==0, aru, distance))

files2$recognizer <- factor(files2$recognizer, levels=c("near training data", "midrange training data", "far training data", "mixed training data"))
summary(files2$recognizer)

hist <- ggplot(files2) +
  geom_histogram(aes(distance), bins=50, fill="grey40") +
  my.theme +
  facet_wrap(~recognizer) +
  ylab("Number of clips") +
  xlab("Distance (m)") +
  scale_fill_manual(guide=FALSE) +
  my.theme
hist

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/StatisticalValidation/Writing/Bioacoustics_Submission")
ggsave("Fig2.png", hist, width=6, height=6, dpi=600)

files.0 <- subset(files2, recognizer=="near training data")
summary(files.0$distance)
sd(files.0$distance)

#G. SPECTROGRAM----

library(tuneR)
library(signal)
library(seewave)
library(sound)

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/Score/Analysis")

near <- readWave("3272_s_EDR-04-0000_20160713_220000.wav")
mid <- readWave("3272_s_EDR-04-0120_20160713_214400.wav")
far <- readWave("3272_s_EDR-04-0240_20160713_214400.wav")

#rec.down <- downsample(rec, 11025)
#spectro(rec.down)
#spectro(rec, wl=256, wn="hanning", flim=c(2.5,7.5), ovlp = 50, palette = reverse.gray.colors.2)

n <- ggspectro(near, flim=c(2.5, 7.5), wl=256) +
  geom_tile(aes(fill=amplitude)) +
  scale_fill_continuous(name="Amplitude\n(dB)\n", na.value="transparent", low="white", high="black", limits=c(-45,0)) +
  my.theme

m <- ggspectro(mid, flim=c(2.5, 7.5), wl=256) +
  geom_tile(aes(fill=amplitude)) +
  scale_fill_continuous(name="Amplitude\n(dB)\n", na.value="transparent", low="white", high="black", limits=c(-45,0)) +
  my.theme

f <- ggspectro(far, flim=c(2.5, 7.5), wl=256) +
  geom_tile(aes(fill=amplitude)) +
  ggtitle("Far") +
  scale_fill_continuous(name="Amplitude\n(dB)\n", na.value="transparent", low="white", high="black", limits=c(-45,0)) +
  my.theme

ggarrange(n,m,f, ncol=3, common.legend=TRUE, legend="right")

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/Score/Analysis")
all <- appendSample("3272_s_EDR-04-0000_20160713_220000.wav", "3272_s_EDR-04-0120_20160713_214400.wav", "3272_s_EDR-04-0240_20160713_214400.wav")
saveSample(all, "3272_concatenated.wav")
all.2 <- readWave("3272_concatenated.wav")

all.3 <- appendSample("3751_s_EDR-07-0000_20160716_214000.wav", "3751_s_EDR-07-0120_20160716_214000.wav", "3751_s_EDR-07-0240_20160716_214000.wav")
saveSample(all.3, "3751_concatenated.wav")
all.4 <- readWave("3751_concatenated.wav")

fig1 <- ggspectro(all.2, flim=c(2.5, 7.5), wl=256) +
  geom_tile(aes(fill=amplitude)) +
  scale_fill_continuous(name="Amplitude\n(dB)\n", na.value="transparent", low="white", high="black", limits=c(-45,0)) +
  annotate("text", x=0, y=8, label="Near (5 m)", hjust=0, size=4) +
  annotate("text", x=0.5, y=8, label="Midrange (117 m)", hjust=0, size=4) +
  annotate("text", x=1.0, y=8, label="Far (127 m)", hjust=0, size=4) +
  my.theme
fig1

setwd("/Users/ellyknight/Documents/UoA/Projects/Projects/Score/Writing/Bioacoustics_Revision")
ggsave("Fig1.png", fig1, width=6, height=3, dpi=600)

websitefig <- ggspectro(all.2, flim=c(2.5, 7.5), wl=256) +
  geom_tile(aes(fill=amplitude)) +
  scale_fill_continuous(name="Amplitude\n(dB)\n", na.value="transparent", low="white", high="black", limits=c(-45,0)) +
  my.theme +
  theme(legend.position = "none")

setwd("/Users/ellyknight/Documents/UoA/Presentations/Materials/Figures")
ggsave("Distance.png", websitefig, width=6, height=6, dpi=600)
