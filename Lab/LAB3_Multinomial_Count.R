#set your own path
setwd("~")
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Multinomial Regression ---------------


# * Proportional Odds Model ----------------

# performs logistic regression for the first m categories combined versus the rest, m = 1, ..., M-1.
# restricts the slope parameters to be the same for all submodels;
# only the intercept coefficient increases in submodel/category index m.

library(MASS)
?polr
# The response must be a factor!

## two formats of multinomial data
data(housing)
head(housing)
Sat <- matrix(housing$Freq, byrow=TRUE, ncol=3)
colnames(Sat) <- housing$Sat[1:3]
housing2 <- data.frame(Sat, housing[seq(1,nrow(housing),3),2:4])
head(housing2)
## Note that: if the observations for some categories are missing, you can impute it by 0 

# if the dataset is in the same format as housing2, you need to turn it into the format as housing to be fitted by polr().

house.plr <- polr(formula = Sat ~ Infl + Type + Cont, data = housing, weights = Freq)
summary(house.plr)
# n by M matrix of predicted prob:
prd_prob_po = predict(house.plr, type='prob')
# equivalently
prd_prob_po = fitted(house.plr)
# vector of predicted labels:
prd_labl_po = predict(house.plr)

# ** Pearson residuals ----------------

if(0){
  # wrong code:
  obslabel <- matrix(0, nrow = nrow(housing), ncol = length(levels(housing$Sat)))
  for (i in seq_len(nrow(housing))) obslabel[i,which(colnames(prd_prob_po)==housing$Sat[i])] = 1
  head(data.frame(obslabel,housing$Sat))
  resP.plr <- (obslabel - prd_prob_po) / sqrt(prd_prob_po * (1 - prd_prob_po))
}

## key: base on your proportional odds model to construct pearson residuals
obslabel <- t(apply(housing2[,1:3], 1, function(x) {
  res <- numeric(3)
  res[which.max(x)] <- 1
  res
}))

resP.plr <- sapply(1:(ncol(obslabel)-1), function(m) {
  obs_m <- rowSums(as.matrix(obslabel[,1:m]))
  fit_m <- rowSums(as.matrix(prd_prob_po[seq_len(nrow(housing2))*3,1:m]))
  (obs_m - fit_m) / sqrt(fit_m * (1 - fit_m))
})

## diagnostics based on Pearson residual ##
fitted_m = sapply(1:(ncol(obslabel)-1), function(m){
  rowSums(as.matrix(prd_prob_po[seq_len(nrow(housing2))*3,1:m]))
})

plot(fitted_m[,1], resP.plr[,1], pch=16, cex=0.6, ylab='Pearson Residuals', xlab='Fitted Values')
abline(h=0, lty=2, col='grey')

plot(fitted_m[,2], resP.plr[,2], pch=16, cex=0.6, ylab='Pearson Residuals', xlab='Fitted Values')
abline(h=0, lty=2, col='grey')


# * Baseline Odds Model ----------------

library(nnet)
?nnet::multinom
# The response can be:
# either a factor vector
# or (grouped data format) a matrix of numbers of subjects falling in certain categories (aggregated format according to shared predictor levels).

housing.bo <- multinom(Sat ~ Infl + Type + Cont, data=housing, weights = Freq)

housing.bo2 <- multinom(cbind(Low, Medium, High) ~ ., data = housing2)

## compare with fitting 2 logistic regression models ##
margin_1 = glm(cbind(Medium, Low) ~ Infl + Type + Cont,
    data = housing2, family = binomial())
margin_2 = glm(cbind(High, Low) ~ Infl + Type + Cont,
    data = housing2, family = binomial())

summary(margin_1)
summary(margin_2)

# The summary table only provides coefficient estimates and corresponding standard errors.
# You can utilize the Gaussian approximation to perform simple hypothesis testing and obtain approximate p-values.
summary(housing.bo, digit=3)
summary(housing.bo2, digit=3)
# z values
zval.bo <- coef(housing.bo) / summary(housing.bo)$standard.errors
# two-sided p-values
pval.bo <- 2 * pnorm(abs(zval.bo), lower.tail=FALSE)

prd_prob_bo = predict(housing.bo, type = 'prob')
prd_prob_bo = fitted(housing.bo)
head(prd_prob_bo)
prd_labl_bo = predict(housing.bo)
head(prd_labl_bo)

prd_prob_bo2 <- fitted(housing.bo2)

# ** Pearson residuals ----------------

# wrong code:
#resP.bo = (obslabel - prd_prob_bo) / sqrt(prd_prob_bo * (1 - prd_prob_bo))

# a list of (M-1) elements, each element contains the Pearson residuals for one submodel
resP.bo <- sapply(2:ncol(obslabel), function(m) {
  # baseline is column 1 here 
  # otherwise you should replace "1" with the corresponding index and adjust the range of "m" accordingly
  obs_m <- obslabel[rowSums(obslabel[,c(1,m)]) > 0, m]
  fit_m <- prd_prob_bo2[rowSums(obslabel[,c(1,m)]) > 0,c(1,m)]
  fit_m <- fit_m[,2] / rowSums(fit_m)
  (obs_m - fit_m) / sqrt(fit_m * (1 - fit_m))
})

# Don't use housing.bo$residuals!
head(resP.bo)
head(housing.bo2$residuals)
head(obslabel - prd_prob_bo2)

## model diagnostics ##
# follow the model diagnostics for logistic regression and Cook's distance in Lab2.

# * Domain shift issue in Multinomial regression model -------------
n = 100
x = c(rnorm(n, -3, 1), rnorm(n, 0, 1), rnorm(n, 3, 1))
y = rep(c(1,2,3), each = n)

library(ggplot2)
library(dplyr)
plot_df = data.frame(x = x,y=as.factor(y))
ggplot(plot_df) + geom_histogram(aes(x = x,  fill = y),color="#e9ecef", binwidth = 0.2, alpha=0.6, position='identity') 

## use the figure to understand the distribution of your predictor associated with each category
# if the domain does not overlap, there might be some convergence issue for glm (check the code below)
glm(y ~ x, family = binomial(), data = plot_df %>% dplyr::filter(y %in% c(1,2)))
glm(y ~ x, family = binomial(), data = plot_df %>% dplyr::filter(y %in% c(1,3))) # convergence issue 

glm(y ~ x, family = binomial(), data = plot_df %>% dplyr::filter(y %in% c(1,2)))
glm(y ~ x, family = binomial(), data = plot_df %>% dplyr::filter(y %in% c(2,3))) # change baseline will work
# several way to resolve the issue
# option1: change the baseline to category 2
# option2: use a threshold to filter on the dataset before the modeling 
# option3: use proportional odds model


# Count data ----------------

# * Overdispersion in Poisson Model ----------------

# Estimate of the overdispersion parameter = Deviance / (n-p)  #method of moment

# ** Negative Binomial Model ----------------
library(MASS)
?glm.nb()

# glm(formula, family = negative.binomial(theta))
# must specify the overdispersion parameter theta

library(pscl) # containing the bioChemists data
data("bioChemists")
n <- nrow(bioChemists)
# compare Poisson and Negative Binomial Model
poisson <- glm(art ~ ., data = bioChemists, family = poisson())
sigma2 <- poisson$deviance/(n-length(coef(poisson)))
sigma2 # 1.798
# fitting nb without specifying theta
nb1 <- glm.nb(art ~ ., data=bioChemists)
summary(nb1)
nb1$deviance / (n-length(coef(nb1))) # 1.105
#estimate theta based on ML estimator
nb1$theta

# fitting nb with a specified theta
nb2 <- glm(art ~ ., data = bioChemists, family=negative.binomial(nb1$theta))
summary(nb2)
nb2$deviance / (n-length(coef(nb2))) # 1.105

# *** Why are the AICs of nb1 and nb2 different? ----------------
# nb1 needs to estimate theta while theta is specified in nb2

# ** Quasi-Poisson Model ----------------

quasip <- glm(art ~ ., data = bioChemists, family = quasipoisson())
summary(quasip)
quasip$deviance / (n-length(coef(quasip))) # 1.798, same as Poisson model
summary(quasip)$dispersion # 1.829, not 1 as fixed in Poisson model

# * Zero-Inflation ----------------

# ** Zero-Inflated Poisson and Negative Binomial Model ----------------

library(pscl)
?zeroinfl

hist(bioChemists$art)
boxplot(bioChemists$art)
quantile(bioChemists$art)

# Zero-Inflated Poisson Model
zip <- zeroinfl(art ~ . | ., data = bioChemists, dist = 'poisson')
summary(zip)
# Zero-Inflated Negative Binomial Model
zinb <- zeroinfl(art ~ . | ., data = bioChemists, dist = 'negbin')
summary(zinb)

residuals(zip, type = "response")
residuals(zip, type = "pearson")

plot(zip$fitted.values, residuals(zip, type = "pearson"))
plot(zinb$fitted.values, residuals(zinb, type = "pearson"))

plot(zip$fitted.values, residuals(zip, type = "response"))
plot(zip$fitted.values, zip$residuals)

# ** Test for the Existence of Zero-Inflation in R ----------------

# likelihood ratio test
# 2(l(full) - l(null)) converges in distribution to \chi^2_{df(full) - df(null)} under the null.

poi <- glm(art ~., data = bioChemists, family = poisson())

1 - pchisq(as.numeric(2*(logLik(zip) - logLik(poi))),
           df = length(coef(zip)) - length(coef(poi)))
