library(MASS)
data(oats)
names(oats) = c('block', 'variety', 'nitrogen', 'yield')
oats$mainplot = oats$variety
oats$subplot = oats$nitrogen

summary(oats)
str(oats)

library(nlme)
m1.nlme = lme(yield ~ variety*nitrogen,
              random = ~ 1|block/mainplot,
              data = oats)

summary(m1.nlme)

anova(m1.nlme)

## linear model
library(nlme)
m1.nlme = lme(Time1stShock  ~ APA * Year ,
              random = ~ 1|ID,
              data = maddyWTtrained)
summary(m1.nlme)
anova(m1.nlme)

m1.nlme = lme(Time1stShock  ~ APA * Genotype ,
              random = ~ 1|ID/Year,
              data = maddy)
summary(m1.nlme)
anova(m1.nlme)


