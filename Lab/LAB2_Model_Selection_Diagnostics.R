#set your own path
setwd("~")


# Model Selection ----------------

# * AIC Type Criteria ----------------

# local optimum
# May need to try different initial models
# Where do AIC and BIC comes from? KL divergence/ Bayesian theorem
# When do you want to use AIC? What about BIC?

library(MASS)
example(birthwt) # birth weight example
birthwt.glm <- glm(low~age+lwt+race+smoke+ptd,
                   family = binomial(), data = bwt)
birthwt.step <- stepAIC(birthwt.glm, trace = FALSE)  # default direction is backward when scope is missing
birthwt.step$anova
stepAIC(birthwt.glm, trace = TRUE) # trace = TRUE to print out fitting process

# define a wider scope with more predictors in the upper model
Scope = list(upper = ~age+lwt+race+smoke+ptd+ht+ui+ftv, lower = ~1)
birthwt.step2 <- stepAIC(birthwt.glm, trace = TRUE, scope = Scope) # default direction is both when scope is specified
birthwt.step2$anova

#by choosing k = log(n), we can use BIC for the model selection


# * Likelihood Ratio Test ----------------

# D(smaller) - D(larger) ~ \chi^2_{df(larger) - df(smaller)}

bwtfit <- glm(formula = low ~ lwt + race + smoke + ptd, family = binomial(),
              data = bwt)
h0.fit = glm(low~lwt + race + smoke, family=binomial(),
             data = bwt)
anova(h0.fit, bwtfit, test="Chi")


# * Deviance Table ----------------

anova(glm(low~smoke+ptd+lwt, family = binomial(), data = bwt), test="Chi")

anova(glm(low~lwt+ptd+smoke, family = binomial(), data = bwt), test="Chi")

# Inconsistent result!
# Using deviance tables for model selection only when p is small;
# Needs to compare all the deviance tables corresponding to different orders of predictors

# Model Diagnostics ----------------

bwtfit <- glm(formula = low ~ lwt + race + smoke + ptd, family = binomial(),
              data = bwt)

# * Deviance Residuals and Pearson Residuals ----------------

# We expect these two types of residuals have similar distributions.
# no lack-of-fit => similar boxplots
# similar boxplots -> next step: Residual Plots

res.P = residuals(bwtfit, type="pearson")
res.D = residuals(bwtfit, type="deviance") #or residuals(fit), by default
boxplot(cbind(res.P, res.D), names = c("Pearson", "Deviance"))

# * Residual Plots ----------------

# no lack-of-fit => no systematic pattern
# next step: Runs Test

par(mfrow=c(1,2))
plot(bwtfit$fitted.values, res.P, pch=16, cex=0.6, ylab='Pearson Residuals', xlab='Fitted Values')
lines(smooth.spline(bwtfit$fitted.values, res.P, spar=0.9), col=2)
abline(h=0, lty=2, col='grey')
plot(bwtfit$fitted.values, res.D, pch=16, cex=0.6, ylab='Deviance Residuals', xlab='Fitted Values')
lines(smooth.spline(bwtfit$fitted.values, res.D, spar=0.9), col=2)
abline(h=0, lty=2, col='grey')

# * Runs Test ----------------

# Null hypothesis: no systematic pattern
# Reject H0 => lack-of-fit
# In the plot: consecutive positives/negatives => systematic pattern

library(lawstat)
# please pay attention to the library
# there are different runs.test() functions in different packages
# we are specifically using this one!

runs.test(y = res.P, plot.it = TRUE)
title(main='Pearson Residual Runs Test')
runs.test(y = res.D, plot.it = TRUE)
title(main='Deviance Residual Runs Test')

# * Check Influential Points & Outliers ----------------

# ** Leverage Points ----------------

# leverage points => influential points

leverage = hatvalues(bwtfit)

W = diag(bwtfit$weights)
X = cbind(rep(1,nrow(bwt)), bwt[['lwt']], bwt[['race']]=='black',
          bwt[['race']]=='other', bwt[['smoke']], bwt[['ptd']])
Hat = sqrt(W) %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% sqrt(W)
all(abs(leverage - diag(Hat)) < 1e-15)

plot(names(leverage), leverage, xlab="Index", type="h")
points(names(leverage), leverage, pch=16, cex=0.6)
p <- length(coef(bwtfit))
n <- nrow(bwt)
abline(h=2*p/n,col=2,lwd=2,lty=2)
infPts <- which(leverage>2*p/n)

# ** Cook's Distance ----------------

# high Cook's distance => influential points/outliers
# leverage points with high Cook's distance => suspicious influential points & outliers
#                    may need to be deleted -> check scatterplots

cooks = cooks.distance(bwtfit)

plot(cooks, ylab="Cook's Distance", pch=16, cex=0.6)
points(infPts, cooks[infPts], pch=17, cex=0.8, col=2)
susPts <- as.numeric(names(sort(cooks[infPts], decreasing=TRUE)[1:3]))
text(susPts, cooks[susPts], susPts, adj=c(-0.1,-0.1), cex=0.7, col=4)

dispersion <- 1
all(abs(cooks - (res.P/(1 - leverage))^2 * leverage/(dispersion * p) < 1e-15))

# Alternative Formulation of Logistic Regression in R ----------------

# In the previous example, the binary response is of individual (0-1) form
bwt$low

# Alternative grouped form: (# of successes, # of failures)
if(0){
  example <- data.frame(lapply(esoph,as.numeric))
  example <- example[,c(4:5,1:3)]
  write.table(example, "example.dat", row.names = FALSE, col.names = TRUE)
}
# Read the data
d = read.table("example.dat", header=TRUE)
names(d) = c("Yes", "No", "A", "B", "C") # set var names
for(i in 3:5){ d[,i] = as.factor(d[,i]) } # convert to factor

# Alternatively, this can be done by
d = read.table("example.dat", header=TRUE, 
               colClasses=c("integer","integer", "factor","factor","factor"))
names(d) = c("Yes", "No", "A", "B", "C")

# Fit logistic regression via glm
# main effect only
fit1 = glm(cbind(Yes, No) ~ A + B + C, data = d, family = binomial())
# main effect & up to 2-way interactions
fit2 = glm(cbind(Yes, No) ~ (A + B + C) * (A + B + C),
           data = d, family = binomial())
# Recall that A * B <=> A + B + A:B, main effect plus interaction
# A : B : C means three-way interaction term only
# A * B * C means main effect, all 2-way plus 3-way interaction

residuals(fit1, type="pearson")

