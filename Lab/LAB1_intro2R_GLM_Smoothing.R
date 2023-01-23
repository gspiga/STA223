rm(list = ls())

#set your own path
setwd("~/Documents/223TA/lab")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set the working path to the path of the current script
#### Intro to R ####

# order of loading
# If two packages use the same function name, then the package loaded last will hide the function from earlier packages
library(plyr)
library(Hmisc)

?dplyr::summarize
?Zelig::summarize


longley
write.table(longley, row.names = FALSE, file='longley.txt')

d = read.table('longley.txt', header=TRUE)

d # View the data
view(d)
names(d) # Check the variable names in the data
summary(d) # summary of mean, extremes, and quartiles for all variables
pairs(d) # Obtain the pairwise scatter plot

varnames = names(d)
names(d) = c(paste0('x', 1:6), 'y')

fit1 = lm(y ~ . , data = d)
fit1 = lm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data = d)

fit2 = lm(y ~ 0 + ., data = d)
fit2 = lm(y ~ . - 1, data = d)

fit3 = lm(y ~ x1 + x3, data = d)

fit4 = lm(y ~ x1 * x2 * x3, data = d)
fit4 = lm(y ~ x1 + x2 + x3 + x1:x2 + x1:x3 + x2:x3 + x1:x2:x3, data = d)

fit5 = lm(y ~ x1 + x2 + x5 + x1:x2, data = d)

fit1
names(fit1)
fit1$fitted.values
fit1$residuals

summary(fit1)
summary(aov(fit1))

library(MASS)
stepAIC(fit1)

par(mfrow= c(1,2))
plot(fit1$fitted.values, fit1$residuals, xlab = "Fitted Values",
     ylab = "Residuals", main = "Residuals Vs Fitted Values")
abline(0,0, col = "red")
qqnorm(fit1$residuals)
qqline(fit1$residuals, col = "red")

plot(fit1,1)
plot(fit1,2)

par(mfrow=c(2,2))
plot(fit1)

save(fit1, file = "longleyfit1.Rdata") # in current working directory

save(list=ls(), file = "longleyfit.Rdata")

load("longleyfit.Rdata")

objects()

rm(fit2)

x = c(2,3,4,NA)
mode(x)
y = c(2,3,4,'2')
mode(y)

is.na(NA)
is.na(NaN)
is.nan(NA)

#### GLM Example 1: Poisson Regression ####

## Dobson (1990) Page 93: Randomized Controlled Trial :

counts = c(18,17,15,20,10,20,25,13,12)
outcome = gl(3,1,9) #generate a factor of 3 levels
treatment = gl(3,3) #generate a factor of 3 levels
d = print(data.frame(treatment, outcome, counts))

glm.D93 = glm(counts ~ outcome + treatment, family=poisson())

#model selection through deviance table
anova(glm.D93, test = "Chisq") # test can take "F" or "Cp" as well

library(MASS)
stepAIC(glm.D93) # model selection through AIC criterion
summary(glm.D93) # parameter estimates

residuals(glm.D93)
print(glm.D93)
coef(glm.D93)
deviance(glm.D93)
AIC(glm.D93)
logLik(glm.D93)
fitted(glm.D93)

?summary.glm
?residuals.glm

?family

vcov(glm.D93)

#### GLM Example 2: Binomial Regression ####

?esoph
summary(esoph)

# effects of alcohol, tobacco and interaction, age adjusted
model1 = glm(cbind(ncases, ncontrols) ~ agegp +  tobgp * alcgp,
             data = esoph, family = binomial())
anova(model1)

#### Smoothing ####

## Cubic Smoothing Splines:
?smooth.spline

set.seed(1)
x = rnorm(100); y = x^2+rnorm(100) # simulate a quadratic relationship
fit.cv = smooth.spline(x,y, cv = TRUE) # CV fit
par(mfrow=c(1,1))
plot(x,y, pch=16, cex=0.6); lines(smooth.spline(x, y, cv = TRUE), col = "red")
fit.cv$spar # the smoothing parameter corresponding to the CV choice
newspar = fit.cv$spar + 0.5
lines(smooth.spline(x,y, spar = newspar), col = "blue") # newspar fit
#degrees of freedom
# plot(x,y, pch=16, cex=0.6); lines(smooth.spline(x, y, df = 8, lambda = 0.001), col = "red")
# fit.cv2 = smooth.spline(x, y, df = 8,lambda = 0.001)
# fit.cv2
## Nadaraya-Watson Kernel Estimators:
?ksmooth

attach(cars); plot(speed, dist, main="Kernel Estimator Example", pch=18)
lines(ksmooth(speed, dist, "normal", bandwidth=2), col="red", lwd=2)
lines(ksmooth(speed, dist, "normal", bandwidth=5), col="green", lwd=2)
lines(ksmooth(speed, dist, "box", bandwidth=10), col="blue", lwd=2)

## B-Spline
library(splines)
?bs

attach(women) # Women height-weight data
# Visualization of the spline basis functions
plot(58:72,bs(women$height, df = 5)[,1], type='l',ylim=c(0,1),
     xlab='x', ylab='spline basis functions')
for(i in 2:5){lines(58:72,bs(women$height, df = 5)[,i])}
# Fit a linear model with B-splines
summary(fm1 <- lm(weight ~ bs(height, df = 5), data = women))
# check the prediction
plot(women, xlab = "Height (in)", ylab = "Weight (lb)",
     main="B-Spline Example", pch=18)
ht <- seq(57, 73, len = 200)
lines(ht, predict(fm1, data.frame(height=ht)), lwd=2, col="red")

## Local Polynomial Smoothing
library(locfit)
?locfit
?lp

data(ethanol, package="locfit") # ethanol gas emission data
head(ethanol)
fit1 <- locfit(NOx ~ E, data=ethanol)
summary(fit1)
plot(fit1, get.data=TRUE, main="local poly example")
# a bivariate local regression with smaller smoothing parameter
fit2 <- locfit(NOx~lp(E,C,nn=0.4,scale=0), data=ethanol)
plot(fit2, main="2D local poly fit")


library(KernSmooth)
?locpoly
mean = locpoly(ethanol$E, ethanol$NOx, bandwidth = 0.1, drv = 0)
deriv = locpoly(ethanol$E, ethanol$NOx, bandwidth = 0.1, drv = 1)
with(ethanol, {plot(E, NOx, ylim=range(deriv), pch=16, cex=0.6); lines(mean$x, mean$y, col='green'); lines(deriv$x, deriv$y, col='red')})

## Local Polynomial Ridge Regression
library(lpridge)
?lpridge

attach(cars); plot(speed, dist, main = "lpRIDGE Regression")
myfit <- lpridge(speed,dist,bandw = 5, ridge=0) # no ridging
lines(myfit$x.out,myfit$est,col=2)
myridge <- lpridge(speed,dist,bandw = 5) # slight ridging
lines(myridge$x.out,myridge$est,col=3); mtext("bandw = 5")
legend(5,120, c("ridge = 0", "default ridging"), col = 2:3, lty = 1)

