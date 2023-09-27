OAManuscript
================

This markdown file contain analyses from both experiments in the
manuscript entitled: “Oxalic acid application method and treatment
intervals for reduction of Varroa destructor populations in western
honey bee (Apis mellifera) colonies”. Data collection for both
experiments were conducted by Branden Stanford and Hossam Abou-Sharra,
and advised by Jamie Ellis and Cameron Jack. Cody Prouty led data
analysis and writing.

\#Load packages

``` r
library(lme4)
library(lmerTest)
library(DHARMa)
library(multcomp)
library(lsmeans)
```

\#Function for data summary (mean, sd, se, ci)

``` r
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )
    datac <- rename(datac, c("mean" = measurevar))
    datac$se <- datac$sd / sqrt(datac$N)
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    return(datac)
}
#Credit to Ania Majewska
```

Experiment 1: Optimal application method

In this experiment, we tested whether the method (OA dribble, fogger, or
vaporizer with a control) significantly affected the amount of bees and
brood in colonies, the number of dead bees in front of hives, and the
infestation of Varroa destructor (measured both by mite washes and mite
fall using sticky boards).

\#File organization and data manipulation

``` r
Bees <- read.csv("Bees.csv")
Brood <- read.csv("Brood.csv")
DeadBees <- read.csv("DeadBees.csv")
MiteWash <- read.csv("MiteWash.csv")
StickyBoard <- read.csv("StickyBoard.csv")
MiteWashW <- read.csv("MiteWashW.csv")

Bees$Time <- as.factor(Bees$Time)
Bees$Colony <- as.factor(Bees$Colony)
Brood$Time <- as.factor(Brood$Time)
Brood$Colony <- as.factor(Brood$Colony)
DeadBees$Time <- as.factor(DeadBees$Time)
DeadBees$Colony <- as.factor(DeadBees$Colony)
MiteWash$Colony <- as.factor(MiteWash$Colony)
MiteWash$Time <- as.factor(MiteWash$Time)
StickyBoard$Colony <- as.factor(StickyBoard$Colony)
StickyBoard$Time <- as.factor(StickyBoard$Time)
```

Tests for significance and plots for the manuscript are below.

\#Effect of estimated number of bees

``` r
Bee <- lmer(Bees ~ Treatment*Time + (1|Colony), data=Bees)
anova(Bee)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                  Sum Sq  Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Treatment       7254386  2418129     3 35.842  2.0754   0.12076    
    ## Time           24163425 24163425     1 33.436 20.7384 6.658e-05 ***
    ## Treatment:Time 14001880  4667293     3 33.393  4.0057   0.01537 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(Bee)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Bees ~ Treatment * Time + (1 | Colony)
    ##    Data: Bees
    ## 
    ## REML criterion at convergence: 1205
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.2392 -0.5412 -0.1046  0.2860  2.9329 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Colony   (Intercept) 1508823  1228    
    ##  Residual             1165152  1079    
    ## Number of obs: 76, groups:  Colony, 40
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)             4743.78     517.10   52.56   9.174  1.7e-12 ***
    ## TreatmentDribble        -782.48     731.30   52.56  -1.070  0.28952    
    ## TreatmentFogger         -747.64     731.30   52.56  -1.022  0.31130    
    ## TreatmentVape           -199.04     731.30   52.56  -0.272  0.78655    
    ## Time2                    293.58     482.73   32.20   0.608  0.54734    
    ## TreatmentDribble:Time2   851.06     715.28   33.50   1.190  0.24248    
    ## TreatmentFogger:Time2    315.12     697.36   32.80   0.452  0.65433    
    ## TreatmentVape:Time2     2254.66     697.36   32.80   3.233  0.00279 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) TrtmnD TrtmnF TrtmnV Time2  TrD:T2 TrF:T2
    ## TrtmntDrbbl -0.707                                          
    ## TretmntFggr -0.707  0.500                                   
    ## TreatmentVp -0.707  0.500  0.500                            
    ## Time2       -0.467  0.330  0.330  0.330                     
    ## TrtmntDr:T2  0.315 -0.445 -0.223 -0.223 -0.675              
    ## TrtmntFg:T2  0.323 -0.228 -0.457 -0.228 -0.692  0.467       
    ## TrtmntVp:T2  0.323 -0.228 -0.228 -0.457 -0.692  0.467  0.479

``` r
#Significant interaction, going to split by treatment

sim_bmond <- simulateResiduals(fittedModel = Bee, n = 250)
hist(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plot(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
#Residuals look good

Bees1 <- subset(Bees, Treatment=="Control")
Bees2 <- subset(Bees, Treatment=="Vape")
Bees3 <- subset(Bees, Treatment=="Dribble")
Bees4 <- subset(Bees, Treatment=="Fogger")

Bee1 <- lm(Bees~ Time, data=Bees1) #change number here to flip through treatments
anova(Bee1)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Bees
    ##           Df   Sum Sq Mean Sq F value Pr(>F)
    ## Time       1   430958  430958  0.1845 0.6726
    ## Residuals 18 42038042 2335447

``` r
#Vaporizer is the only significant treatment
```

We opted for a table in the manuscript, these plots are a look behind
the curtain.

\#Number of estimated bees plot

``` r
Beez <- na.omit(Bees[,1:5])
tgcBees <- summarySE(Beez, measurevar="Bees", groupvars=c("Treatment", "Time"))
Plot <- data.frame(Time=c("1", "2"), Control = tgcBees[c(1:2), "Bees"], Dribble = tgcBees[c(3:4), "Bees"], 
                   Fogger = tgcBees[c(5:6), "Bees"], Vape = tgcBees[c(7:8), "Bees"])

#png("BeesOAFV.png", width=6, height=6, units="in", res=400)
#^Code for writing png for figures
plot(.8:1.8, Plot[,2], xaxt="n", ylim=c(3200, 7500), xlim=c(.6,2.4), type="b", pch=16, col="blue", xlab="Time", ylab="Number of bees", cex=1.5, lty=4, cex.lab=1.5, lwd=2)
lines(.9:1.9, Plot[,3], pch=17, col="darkgreen", type="b", cex=1.5,lty=5, lwd=2)
lines(1.0:2.0, Plot[,4], pch=15, col="black", type="b", cex=1.5,lty=6, lwd=2)
lines(1.1:2.1, Plot[,5], pch=18, col="red", type="b", cex=1.5,lty=7, lwd=2)
axis(1, 1:2, labels=c("Pre-Treatment", "Post-Treatment"))
arrows(.8:1.8, Plot[,2]-tgcBees[c(1:2),"se"], .8:1.8, Plot[,2]+tgcBees[c(1:2), "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
arrows(.9:1.9, Plot[,3]-tgcBees[c(3:4),"se"], .9:1.9, Plot[,3]+tgcBees[c(3:4), "se"], code=3, angle=90, length=0.05, lwd=2, col="darkgreen")
arrows(1.0:2.0, Plot[,4]-tgcBees[c(5:6),"se"], 1.0:2.0, Plot[,4]+tgcBees[c(5:6), "se"], code=3, angle=90, length=0.05, lwd=2, col="black")
arrows(1.1:2.1, Plot[,5]-tgcBees[c(7:8),"se"], 1.1:2.1, Plot[,5]+tgcBees[c(7:8), "se"], code=3, angle=90, length=0.05, lwd=2, col="red")
legend("topleft", c("Control", "Dribble", "Fogger", "Vaporizer"), col=c("blue", "darkgreen", "black", "red"), lty=c(4,5,6,7), title="Treatment", pch=c(16,17,15,18, cex=5))
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
#dev.off()
```

Brood estimations were similar to bee estimations.

\#Tests for significance

``` r
Broo <- lmer(Brood ~ Treatment*Time + (1|Colony), data=Brood)
anova(Broo)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Treatment       887736  295912     3 34.711  1.2728  0.298932    
    ## Time           7442837 7442837     1 31.627 32.0128 3.045e-06 ***
    ## Treatment:Time 4737007 1579002     3 31.598  6.7915  0.001156 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant interaction

sim_bmond <- simulateResiduals(fittedModel = Broo, n = 250)
hist(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
plot(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
Brood1 <- subset(Brood, Treatment=="Control")
Brood2 <- subset(Brood, Treatment=="Vape")
Brood3 <- subset(Brood, Treatment=="Dribble")
Brood4 <- subset(Brood, Treatment=="Fogger")

Broo1 <- lm(Brood~ Time, data=Brood1)
anova(Broo1)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Brood
    ##           Df   Sum Sq Mean Sq F value Pr(>F)
    ## Time       1   662316  662316  0.8502 0.3687
    ## Residuals 18 14021414  778967

``` r
#Vaporizer and dribble were both significant, control and fogger ns.
```

\#Plot for brood estimation

``` r
Brooz <- na.omit(Brood)
tgcBrood <- summarySE(Brooz, measurevar="Brood", groupvars=c("Treatment", "Time"))


Plot <- data.frame(Time=c("1", "2"), Control = tgcBrood[c(1:2), "Brood"], Dribble = tgcBrood[c(3:4), "Brood"], 
                   Fogger = tgcBrood[c(5:6), "Brood"], Vape = tgcBrood[c(7:8), "Brood"])

#png("BroodOAFV.png", width=6, height=6, units="in", res=400)
plot(.8:1.8, Plot[,2], xaxt="n", ylim=c(2000, 4500), xlim=c(.6,2.4), type="b", pch=16, col="blue", xlab="Time", ylab="Area of brood (cm^2)", cex=1.5, lty=4, cex.lab=1.5, lwd=2)
lines(.9:1.9, Plot[,3], pch=17, col="darkgreen", type="b", cex=1.5,lty=5, lwd=2)
lines(1.0:2.0, Plot[,4], pch=15, col="black", type="b", cex=1.5,lty=6, lwd=2)
lines(1.1:2.1, Plot[,5], pch=18, col="red", type="b", cex=1.5,lty=7, lwd=2)
axis(1, 1:2, labels=c("Pre-Treatment", "Post-Treatment"))
arrows(.8:1.8, Plot[,2]-tgcBrood[c(1:2),"se"], .8:1.8, Plot[,2]+tgcBrood[c(1:2), "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
arrows(.9:1.9, Plot[,3]-tgcBrood[c(3:4),"se"], .9:1.9, Plot[,3]+tgcBrood[c(3:4), "se"], code=3, angle=90, length=0.05, lwd=2, col="darkgreen")
arrows(1.0:2.0, Plot[,4]-tgcBrood[c(5:6),"se"], 1.0:2.0, Plot[,4]+tgcBrood[c(5:6), "se"], code=3, angle=90, length=0.05, lwd=2, col="black")
arrows(1.1:2.1, Plot[,5]-tgcBrood[c(7:8),"se"], 1.1:2.1, Plot[,5]+tgcBrood[c(7:8), "se"], code=3, angle=90, length=0.05, lwd=2, col="red")
legend("topleft", c("Control", "Dribble", "Fogger", "Vaporizer"), col=c("blue", "darkgreen", "black", "red"), lty=c(4,5,6,7), title="Treatment", pch=c(16,17,15,18, cex=5))
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
#dev.off()
```

Dead bee counts in front of hives

\#Tests for significance

``` r
DB <- lmer(Bees ~ Treatment*Time + (1|Colony), data = DeadBees)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
anova(DB)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
    ## Treatment      132.12  44.039     3   104  0.3086 0.8191
    ## Time           617.11 308.555     2   104  2.1624 0.1202
    ## Treatment:Time 108.63  18.104     6   104  0.1269 0.9928

``` r
#Not significant

sim_bmond <- simulateResiduals(fittedModel = DB, n = 250)
hist(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
plot(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
#It is still good to check residuals :)
```

\#Dead bee plot

``` r
DBs <- na.omit(DeadBees)
tgcDB <- summarySE(DBs, measurevar="Bees", groupvars=c("Treatment", "Time"))

#png("DBOAFV.png", width=6, height=6, units="in", res=400)
Plot <- data.frame(Time=c("1", "2", "1"),Control = tgcDB[c(1:1), "Bees"], 
                   Dribble= tgcDB[c(4:6), "Bees"], Fogger = tgcDB[c(7:9), "Bees"], Vape = tgcDB[c(10:12), "Bees"])

plot(1:3, Plot[,2], xaxt="n", ylim=c(10,30),xlim=c(.8,3.4), type="b", pch=16, col="blue", xlab="Application", ylab="Number of dead bees", cex=1.5, lty=2, cex.lab=1.5, lwd=2)
lines(1:3, Plot[,3], pch=17, col="black", type="b", cex=1.5,lty=2, lwd=2)
lines(1.1:3.1, Plot[,4], pch=15, col="darkgreen", type="b", cex=1.5,lty=2, lwd=2)
lines(1.1:3.1, Plot[,5], pch=18, col="red", type="b", cex=1.5,lty=2, lwd=2)
axis(1, 1:3, labels=c("1", "2", "3"))
arrows(1:3, Plot[,2]-tgcDB[c(1:3),"se"], 1:3, Plot[,2]+tgcDB[c(1:3), "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
arrows(1:3, Plot[,3]-tgcDB[c(4:6),"se"], 1:3, Plot[,3]+tgcDB[c(4:6), "se"], code=3, angle=90, length=0.05, lwd=2, col="black")
arrows(1.1:3.1, Plot[,4]-tgcDB[c(7:9),"se"], 1.1:3.1, Plot[,4]+tgcDB[c(7:9), "se"], code=3, angle=90, length=0.05, lwd=2, col="darkgreen")
arrows(1.1:3.1, Plot[,5]-tgcDB[c(10:12),"se"], 1.1:3.1, Plot[,5]+tgcDB[c(10:12), "se"], code=3, angle=90, length=0.05, lwd=2, col="red")
legend("topright", c("Control", "Dribble", "Fogger", "Vaporizer"), col=c("blue", "black", "darkgreen", "red"), lty=0, title="Time", pch=c(16,17,15,18))
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
#dev.off()
```

Mite infestation - Number of mites per 100 bees

\#Tests for significance -

``` r
MW <- lmer(Per100 ~ Treatment*Time + (1|Colony), data = MiteWash)
anova(MW)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## Treatment       10.292   3.431     3 34.580  0.3099 0.818026   
    ## Time            98.308  98.308     1 33.747  8.8822 0.005309 **
    ## Treatment:Time 146.880  48.960     3 33.727  4.4235 0.009977 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant interaction

sim_bmond <- simulateResiduals(fittedModel = MW, n = 250)
hist(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
plot(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
MiteWash1 <- subset(MiteWash, Treatment=="Control")
MiteWash2 <- subset(MiteWash, Treatment=="Vape")
MiteWash3 <- subset(MiteWash, Treatment=="Dribble")
MiteWash4 <- subset(MiteWash, Treatment=="Fogger")

MW1 <- lm(Per100~ Time, data=MiteWash1)
anova(MW1)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Per100
    ##           Df  Sum Sq Mean Sq F value Pr(>F)
    ## Time       1   0.782  0.7822  0.0641 0.8029
    ## Residuals 18 219.548 12.1971

``` r
#Significant effects of vaporizer and dribble
```

\#Plot for mite wash

``` r
MWs <- na.omit(MiteWash[,1:6])
tgcMW <- summarySE(MWs, measurevar="Per100", groupvars=c("Treatment", "Time"))

Plot <- data.frame(Time=c("1", "2"), Control = tgcMW[c(1:2), "Per100"], Dribble = tgcMW[c(3:4), "Per100"], 
                   Fogger = tgcMW[c(5:6), "Per100"], Vape = tgcMW[c(7:8), "Per100"])

#png("MWOAFV.png", width=6, height=6, units="in", res=400)
plot(.8:1.8, Plot[,2], xaxt="n", ylim=c(2, 12), xlim=c(.6,2.4), type="b", pch=16, col="blue", xlab="Time", ylab="Number of bees", cex=1.5, lty=4, cex.lab=1.5, lwd=2)
lines(.9:1.9, Plot[,3], pch=17, col="darkgreen", type="b", cex=1.5,lty=5, lwd=2)
lines(1.0:2.0, Plot[,4], pch=15, col="black", type="b", cex=1.5,lty=6, lwd=2)
lines(1.1:2.1, Plot[,5], pch=18, col="red", type="b", cex=1.5,lty=7, lwd=2)
axis(1, 1:2, labels=c("Pre-Treatment", "Post-Treatment"))
arrows(.8:1.8, Plot[,2]-tgcMW[c(1:2),"se"], .8:1.8, Plot[,2]+tgcMW[c(1:2), "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
arrows(.9:1.9, Plot[,3]-tgcMW[c(3:4),"se"], .9:1.9, Plot[,3]+tgcMW[c(3:4), "se"], code=3, angle=90, length=0.05, lwd=2, col="darkgreen")
arrows(1.0:2.0, Plot[,4]-tgcMW[c(5:6),"se"], 1.0:2.0, Plot[,4]+tgcMW[c(5:6), "se"], code=3, angle=90, length=0.05, lwd=2, col="black")
arrows(1.1:2.1, Plot[,5]-tgcMW[c(7:8),"se"], 1.1:2.1, Plot[,5]+tgcMW[c(7:8), "se"], code=3, angle=90, length=0.05, lwd=2, col="red")
legend("topright", c("Control", "Dribble", "Fogger", "Vaporizer"), col=c("blue", "darkgreen", "black", "red"), lty=c(4,5,6,7), title="Treatment", pch=c(16,17,15,18, cex=5))
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
#dev.off()
```

We also tested the difference between pre- and post-treatment.

\#Difference in MW significance

``` r
MiteWashW$Change <- (MiteWashW$Per1002 - MiteWashW$Per1001)
MWW <- lm(Change ~ Treatment, data = MiteWashW)
anova(MWW)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Change
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## Treatment  3 295.84  98.613  4.5492 0.008938 **
    ## Residuals 33 715.34  21.677                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant effect

sim_bmond <- simulateResiduals(fittedModel = MWW, n = 250)
hist(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
plot(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

\#Plot for the difference in mite wash

``` r
OAMiteW1 <- na.omit(MiteWashW)
tgcMiteW1 <- summarySE(OAMiteW1, measurevar = "Change", groupvars = c("Treatment"))

#png("MitesChangeT.png", width=6, height=6, units="in", res=400)
plot(1:4, tgcMiteW1[,3], xaxt="n", ylim=c(-10, 10), xlim=c(1,4), type="b", pch=16, col="blue", xlab="Treatment", ylab="Difference in mite infestation", cex=1.5, lty=0, cex.lab=1.5, lwd=2)
lines(1:4, c(0,0,0,0), lty=2)
axis(1, 1:4, labels=c("Control", "Dribble", "Fogger", "Vaporizer"))
arrows(1:4, tgcMiteW1[,3]-tgcMiteW1[,"se"], 1:4, tgcMiteW1[,3]+tgcMiteW1[, "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
points(OAMiteW1$Change ~ OAMiteW1$Treatment, pch=19, col="grey", cex=.5)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
#dev.off()
```

Sticky boards are used to catch mite fall. The boards are placed during
treatment and counted 72-hours after each treatment is applied.

\#Tests for significance

``` r
SB <- lmer(Mites ~ Treatment*Time + (1|Colony), data = StickyBoard)
anova(SB)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Treatment       69788   23263     3 30.889  24.070 3.170e-08 ***
    ## Time            39384   19692     2 64.730  20.375 1.369e-07 ***
    ## Treatment:Time  39942    6657     6 64.705   6.888 1.098e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant interaction, since there are multiple time points we will be separating by time instead of treatment.

sim_bmond <- simulateResiduals(fittedModel = SB, n = 250)
hist(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
plot(sim_bmond)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
StickyBoard1 <- subset(StickyBoard, Time==1)
StickyBoard2 <- subset(StickyBoard, Time==2)
StickyBoard3 <- subset(StickyBoard, Time==3)

SB1 <- lm(Mites~ Treatment, data=StickyBoard1)
anova(SB1)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Mites
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## Treatment  3  61437 20479.0   11.33 2.221e-05 ***
    ## Residuals 36  65068  1807.4                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant difference between treatments

SB2 <- lm(Mites~ Treatment, data=StickyBoard2)
anova(SB2)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Mites
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## Treatment  3  85653   28551  21.148 5.425e-08 ***
    ## Residuals 35  47252    1350                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant difference between treatments

SB3 <- lm(Mites~ Treatment, data=StickyBoard3)
anova(SB3)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Mites
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## Treatment  3 4689.7 1563.22  8.9614 0.0001743 ***
    ## Residuals 33 5756.5  174.44                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant difference between treatments
```

\#Sticky board plot

``` r
SBs <- na.omit(StickyBoard[,1:4])
tgcSB <- summarySE(SBs, measurevar="Mites", groupvars=c("Treatment", "Time"))

#png("SBOAFV.png", width=6, height=6, units="in", res=400)
Plot <- data.frame(Time=c("1", "2", "3"),Control = tgcSB[c(1:3), "Mites"], 
                   Dribble= tgcSB[c(4:6), "Mites"], Fogger = tgcSB[c(7:9), "Mites"], Vape = tgcSB[c(10:12), "Mites"])

plot(1:3, Plot[,2], xaxt="n", ylim=c(0,150),xlim=c(.8,3.4), type="b", pch=16, col="blue", xlab="Application", ylab="72-hr mite fall", cex=1.5, lty=2, cex.lab=1.5, lwd=2)
lines(1:3, Plot[,3], pch=17, col="black", type="b", cex=1.5,lty=2, lwd=2)
lines(1.1:3.1, Plot[,4], pch=15, col="darkgreen", type="b", cex=1.5,lty=2, lwd=2)
lines(1.1:3.1, Plot[,5], pch=18, col="red", type="b", cex=1.5,lty=2, lwd=2)
axis(1, 1:3, labels=c("1", "2", "3"))
arrows(1:3, Plot[,2]-tgcSB[c(1:3),"se"], 1:3, Plot[,2]+tgcSB[c(1:3), "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
arrows(1:3, Plot[,3]-tgcSB[c(4:6),"se"], 1:3, Plot[,3]+tgcSB[c(4:6), "se"], code=3, angle=90, length=0.05, lwd=2, col="black")
arrows(1.1:3.1, Plot[,4]-tgcSB[c(7:9),"se"], 1.1:3.1, Plot[,4]+tgcSB[c(7:9), "se"], code=3, angle=90, length=0.05, lwd=2, col="darkgreen")
arrows(1.1:3.1, Plot[,5]-tgcSB[c(10:12),"se"], 1.1:3.1, Plot[,5]+tgcSB[c(10:12), "se"], code=3, angle=90, length=0.05, lwd=2, col="red")
legend("topright", c("Control", "Dribble", "Fogger", "Vaporizer"), col=c("blue", "black", "darkgreen", "red"), lty=0, title="Time", pch=c(16,17,15,18))
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
#dev.off()
```

Experiment 2: Vaporization interval

Unfortunately, there were not as many response variables collected for
this experiment. Chronologically, this experiment was conducted prior to
Experiment 1 and we realized that some information would have been
useful to collect (such as colony estimations and dead bee counts) so
those data were collected in Experiment 1.

\#Data organization

``` r
OAMite <- read.csv("OAMitewashes.csv")
OASticky <- read.csv("OAStickyBoards.csv")
OAMiteW <- read.csv("MiteWashesW.csv")

OAMite$Time <- as.factor(OAMite$Time)
OAMite$Hive <- as.factor(OAMite$Hive)
OASticky$Time <- as.factor(OASticky$Time)
OASticky$Hive <- as.factor(OASticky$Hive)
```

Mite wash - Interaction design

\#Mite wash

``` r
Mite <- lmer(Per100 ~ Treatment * Time + (1|Hive), data= OAMite)
anova(Mite)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Treatment      111.26   37.09     3 35.073  1.6749 0.1901355    
    ## Time           320.81  320.81     1 33.970 14.4889 0.0005624 ***
    ## Treatment:Time 385.46  128.49     3 33.915  5.8029 0.0025860 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant interaction - split by treatment

OAMiteT1 <- subset(OAMite, Treatment == "control")
OAMiteT2 <- subset(OAMite, Treatment == "day3")
OAMiteT3 <- subset(OAMite, Treatment == "day5")
OAMiteT4 <- subset(OAMite, Treatment == "day7")

Mite <- lm(Per100 ~ Time, data= OAMiteT1)
anova(Mite)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Per100
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## Time       1  13.60  13.602  0.2586 0.6176
    ## Residuals 17 894.23  52.602

``` r
#5 and 7 day intervals are significant
```

\#Mite wash plot

``` r
OAMite2 <- na.omit(OAMite[,1:8])
tgcMiteW <- summarySE(OAMite2, measurevar = "Per100", groupvars = c("Treatment","Time"))

Plot <- data.frame(Time=c("1", "2"), Control = tgcMiteW[c(1:2), "Per100"], Day1 = tgcMiteW[c(3:4), "Per100"], 
                   Day5 = tgcMiteW[c(5:6), "Per100"], Day7 = tgcMiteW[c(7:8), "Per100"])

#png("MitesFieldOA100.png", width=6, height=6, units="in", res=400)
plot(.8:1.8, Plot[,2], xaxt="n", ylim=c(0, 14), xlim=c(.6,2.4), type="b", pch=16, col="blue", xlab="Time", ylab="Mites per 100 Bees", cex=1.5, lty=4, cex.lab=1.5, lwd=2)
lines(.9:1.9, Plot[,3], pch=17, col="darkgreen", type="b", cex=1.5,lty=5, lwd=2)
lines(1.0:2.0, Plot[,4], pch=15, col="black", type="b", cex=1.5,lty=6, lwd=2)
lines(1.1:2.1, Plot[,5], pch=18, col="red", type="b", cex=1.5,lty=7, lwd=2)
axis(1, 1:2, labels=c("Pre-Treatment", "Post-Treatment"))
arrows(.8:1.8, Plot[,2]-tgcMiteW[c(1:2),"se"], .8:1.8, Plot[,2]+tgcMiteW[c(1:2), "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
arrows(.9:1.9, Plot[,3]-tgcMiteW[c(1:4),"se"], .9:1.9, Plot[,3]+tgcMiteW[c(1:4), "se"], code=3, angle=90, length=0.05, lwd=2, col="darkgreen")
arrows(1.0:2.0, Plot[,4]-tgcMiteW[c(5:6),"se"], 1.0:2.0, Plot[,4]+tgcMiteW[c(5:6), "se"], code=3, angle=90, length=0.05, lwd=2, col="black")
arrows(1.1:2.1, Plot[,5]-tgcMiteW[c(7:8),"se"], 1.1:2.1, Plot[,5]+tgcMiteW[c(7:8), "se"], code=3, angle=90, length=0.05, lwd=2, col="red")
legend("topright", c("Control", "1 Day", "5 Day", "7 Day"), col=c("blue", "darkgreen", "black", "red"), lty=c(4,5,6,7), title="Treatment", pch=c(16,17,15,18, cex=5))
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
#dev.off()
```

Again in this experiment we also tested the difference in mite
infestations for pre- and post-application.

\#Significance

``` r
OAMiteW$Change <- (OAMiteW$Per1002 - OAMiteW$Per1001) 

MiteW <- lm(Change ~ Treatment, data = OAMiteW)
anova(MiteW)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Change
    ##           Df  Sum Sq Mean Sq F value   Pr(>F)   
    ## Treatment  3  821.96  273.99  6.2326 0.001866 **
    ## Residuals 32 1406.73   43.96                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant effect
```

\##Difference plot

``` r
OAMiteW1 <- na.omit(OAMiteW)

tgcMiteW <- summarySE(OAMiteW1, measurevar = "Change", groupvars = c("Treatment"))

#png("MitesChange.png", width=6, height=6, units="in", res=400)
plot(1:4, tgcMiteW[,3], xaxt="n", ylim=c(-15, 10), xlim=c(1,4), type="b", pch=16, col="blue", xlab="Treatment", ylab="Difference in mite infestation", cex=1.5, lty=0, cex.lab=1.5, lwd=2)
lines(1:4, c(0,0,0,0), lty=2)
axis(1, 1:4, labels=c("Control", "1-Day", "5-Day", "7-Day"))
arrows(1:4, tgcMiteW[,3]-tgcMiteW[,"se"], 1:4, tgcMiteW[,3]+tgcMiteW[, "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
points(OAMiteW1$Change ~ OAMiteW1$Treatment, pch=19, col="grey", cex=.5)
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
#dev.off()
```

Sticky board data for Experiment 2 were analyzed the same as Experiment
1.

\#Significance

``` r
Sticky <- lmer(Mites ~ Treatment * Time + (1|Hive), data= OASticky)
anova(Sticky)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## Treatment      121492   40497     3  35.434  5.1251  0.004761 ** 
    ## Time           745265  248422     3 105.339 31.4386 1.355e-14 ***
    ## Treatment:Time 394527   43836     9 105.288  5.5476 3.062e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant intercation

OASticky1 <- subset(OASticky, Time == 1)
OASticky2 <- subset(OASticky, Time == 2)
OASticky3 <- subset(OASticky, Time == 3)
OASticky4 <- subset(OASticky, Time == 4)

Sticky1 <- lm(Mites ~ Treatment, data= OASticky1)
anova(Sticky1)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Mites
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## Treatment  3 479588  159863  6.1089 0.001808 **
    ## Residuals 36 942085   26169                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant difference between treatments


Sticky2 <- lm(Mites ~ Treatment, data= OASticky2)
anova(Sticky2)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Mites
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## Treatment  3 110650   36883  4.6543 0.007532 **
    ## Residuals 36 285285    7925                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant difference between treatments

Sticky3 <- lm(Mites ~ Treatment, data= OASticky3)
anova(Sticky3)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Mites
    ##           Df Sum Sq Mean Sq F value  Pr(>F)  
    ## Treatment  3  21037  7012.5  3.0929 0.03901 *
    ## Residuals 36  81623  2267.3                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Significant difference between treatments

Sticky4 <- lm(Mites ~ Treatment, data= OASticky4)
anova(Sticky4)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Mites
    ##           Df  Sum Sq Mean Sq F value Pr(>F)
    ## Treatment  3    97.8   32.62  0.0678 0.9766
    ## Residuals 33 15869.7  480.90

``` r
#No significant difference between treatments
```

\##SB Graph

``` r
OASB2 <- na.omit(OASticky[,1:6])

tgcSB <- summarySE(OASB2, measurevar = "Mites", groupvars = c("Treatment", "Time"))

Plot <- data.frame(Time=c("1", "2", "3", "4"), Control = tgcSB[c(1:4), "Mites"], Day3 = tgcSB[c(5:8), "Mites"], 
                   Day5 = tgcSB[c(9:12), "Mites"], Day7 = tgcSB[c(13:16), "Mites"])


#png("MitesStickyboardOA.png", width=6, height=6, units="in", res=400)
plot(.8:3.8, Plot[,2], xaxt="n", ylim=c(20,450), xlim=c(.6,4.4), type="b", pch=16, col="blue", xlab="OA Treatment Number", ylab="Mites dropped", cex=1.5, lty=4, cex.lab=1.5, lwd=2)
lines(.9:3.9, Plot[,3], pch=17, col="darkgreen", type="b", cex=1.5,lty=5, lwd=2)
lines(1.0:4.0, Plot[,4], pch=15, col="black", type="b", cex=1.5,lty=6, lwd=2)
lines(1.1:4.1, Plot[,5], pch=18, col="red", type="b", cex=1.5,lty=7, lwd=2)
axis(1, 1:4, labels=c("1", "2", "3", "4"))
arrows(.8:3.8, Plot[,2]-tgcSB[c(1:4),"se"], .8:3.8, Plot[,2]+tgcSB[c(1:4), "se"], code=3, angle=90, length=0.05, lwd=2, col="blue")
arrows(.9:3.9, Plot[,3]-tgcSB[c(5:8),"se"], .9:3.9, Plot[,3]+tgcSB[c(5:8), "se"], code=3, angle=90, length=0.05, lwd=2, col="darkgreen")
arrows(1.0:4.0, Plot[,4]-tgcSB[c(9:12),"se"], 1.0:4.0, Plot[,4]+tgcSB[c(9:12), "se"], code=3, angle=90, length=0.05, lwd=2, col="black")
arrows(1.1:4.1, Plot[,5]-tgcSB[c(13:16),"se"], 1.1:4.1, Plot[,5]+tgcSB[c(13:16), "se"], code=3, angle=90, length=0.05, lwd=2, col="red")
legend("topright", c("Control", "3 Day", "5 Day", "7 Day"), col=c("blue", "darkgreen", "black", "red"), lty=c(4,5,6,7), title="Treatment", pch=c(16,17,15,18, cex=5))
```

![](OAManuscript_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
#dev.off()
```
