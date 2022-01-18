##

#Analysis for lithium chloride manuscript

#Contact: cprouty@ufl.edu    https://scholar.google.com/citations?user=PpeDx78AAAAJ&hl=en

##

#load packages
library(lme4)
library(lsmeans)
library(multcomp)
library(afex)
library(lme4)
##

#Import files
OAMite <- read.csv("OAMitewashes.csv")
OASticky <- read.csv("OAStickyBoards.csv")
##

#Data organization
OAMite$Time <- as.factor(OAMite$Time)
OAMite$Hive <- as.factor(OAMite$Hive)
OASticky$Time <- as.factor(OASticky$Time)
OASticky$Hive <- as.factor(OASticky$Hive)
OAMiteT1 <- subset(OAMite, Time == 1)
OAMiteT2 <- subset(OAMite, Time == 2)
OASticky1 <- subset(OASticky, Time == 1)
OASticky2 <- subset(OASticky, Time == 2)
OASticky3 <- subset(OASticky, Time == 3)
OASticky4 <- subset(OASticky, Time == 4)
###

#Effect of treatment on mite infestation - alcohol wash
Mite <- lmer(Per100 ~ Treatment * Time + (1|Hive), data= OAMite)
anova(Mite)

lsm<-lsmeans (Mite, list( ~ Treatment + Time))
cld(lsm)
##

Mite1 <- lm(Per100 ~ Treatment, data= OAMiteT1)
anova(Mite1)

lsm<-lsmeans (Mite1, list( ~ Treatment))
cld(lsm)
##

Mite2 <- lm(Per100 ~ Treatment, data= OAMiteT2)
anova(Mite2)

lsm<-lsmeans (Mite2, list( ~ Treatment))
cld(lsm)
###

#Effect of treatment on mite drop - sticky board data
Sticky <- lmer(Mites ~ Treatment * Time + (1|Hive), data= OASticky)
anova(Sticky)

lsm<-lsmeans (Sticky, list( ~ Treatment + Time))
cld(lsm)
##

Sticky1 <- lmer(Mites ~ Treatment, data= OASticky1)
anova(Sticky1)

lsm<-lsmeans (Sticky1, list( ~ Treatment))
cld(lsm)
##

Sticky2 <- lmer(Mites ~ Treatment, data= OASticky2)
anova(Sticky2)

lsm<-lsmeans (Sticky2, list( ~ Treatment))
cld(lsm)
##

Sticky3 <- lmer(Mites ~ Treatment, data= OASticky3)
anova(Sticky3)

lsm<-lsmeans (Sticky3, list( ~ Treatment))
cld(lsm)
##

Sticky4 <- lmer(Mites ~ Treatment, data= OASticky4)
anova(Sticky4)

lsm<-lsmeans (Sticky4, list( ~ Treatment))
cld(lsm)
###