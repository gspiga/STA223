
# Classification ----------------

# Both LDA and QDA requires Gaussian assumption.

# decision boundary for LDA / QDA / logistic regression is linear or not? 

# * Linear Discriminant Analysis ----------------
# requires that the two groups share the same covariance matrix
# Sometimes we can only use LDA but not QDA due to degrees of freedom issue.

library(MASS)
?lda

head(iris)
table(iris$Species)

# n=150, p=4, N_s = N_c = N_v = 50
train <- sample(1:150, 75)
z <- lda(Species ~ ., iris, prior = c(1,1,1)/3, subset = train)
pred <- predict(z, iris[-train, ])$class
cl.test <- iris$Species[-train]
table(pred, cl.test)
# To get the misclassification rate
misrate <- length(which(cl.test != pred))/length(cl.test)

zcv <- lda(Species ~ ., iris, prior = c(1,1,1)/3, CV = TRUE)
# it contains CV-predicted class label and fitted probabilities
pred_cv <- zcv$class
table(pred_cv, iris$Species)

# * Quadratic Discriminant Analysis ----------------
# when the two groups have different covariance matrices

tr <- sample(1:50, 25)
train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
test <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
cl <- factor(c(rep("Setosa",25), rep("Versicolor",25), rep("Virginica",25)))
# class labels are the same for training and testing
zq <- qda(train, cl)
predq <- predict(zq,test)$class
table(cl, predq)
# To get misclassification rate
misrateq <- length(which(cl != predq))/nrow(test)

zqcv <- qda(train, cl, CV = TRUE)
pred_qcv <- zqcv$class
table(cl, pred_qcv)

# Cross Validation for Model Selection and Model Evaluation ----------------

# * Elastic Net GLM for Model Selection ----------------

library(glmnet)
?glmnet # fit a model with a set of lambda
?cv.glmnet # K-fold cross validation for a saturated model to choose lambda

# generate data
set.seed(1010)
n = 1000; p = 100
nzc = trunc(p/10) # number of nonzero coefficient columns
x = matrix(rnorm(n*p),n,p) # design matrix without intercept
beta = rnorm(nzc) # nonzero coefficients
fx= x[,seq(nzc)] %*% beta # linear predictor
 # Only the first nzc coefficients are nonzero.
px = exp(fx)
px = px/(1+px) # true probability
ly = rbinom(n=length(px),prob=px,size=1) # binary response vector

set.seed(1011)
cvob_bin=cv.glmnet(x,ly,family="binomial", alpha=1) # default: alpha=1 (lasso)
plot(cvob_bin)
title("Binomial Family",line=2.5) # See Figure 9.
# top row of numbers: number of predictors selected

# select the optimal regularization parameter lambda
opt_lambda_idx = which(cvob_bin$lambda == cvob_bin$lambda.min)
opt_lambda = cvob_bin$lambda[opt_lambda_idx]

# check first 15 of the corresponding coefficients
head(cvob_bin$glmnet.fit$beta[,opt_lambda_idx], n = 15)
# We see that some coefficients are shrunk to zero.

# check the corresponding number of nonzero coefficients
which(abs(cvob_bin$glmnet.fit$beta[,opt_lambda_idx] )> 0)

# Fit with optimal lambda:
# DO NOT do this in practice, as it is double-dipping.
# Fit the model on TEST SET instead.
doubledippingfit = glmnet(x,ly,family="binomial",lambda=opt_lambda)
coef(doubledippingfit)
as.matrix(coef(doubledippingfit))


# Multinomial regression with regularization 
data(MultinomialExample)
x <- MultinomialExample$x
y <- MultinomialExample$y
fit <- glmnet(x, y, family = "multinomial", type.multinomial = "grouped")
coef_fit = coef(fit)
coef_fit1 = as.matrix(coef_fit$`1`)
coef_fit2 = as.matrix(coef_fit$`2`)
coef_fit3 = as.matrix(coef_fit$`3`)

plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")
apply(as.matrix(coef_fit$`1`), 2, sum)


# * Embedded Cross-Validation in R ----------------

#library(cvTools)
library(boot)
?cv.glm

data(mammals, package="MASS")
mammals.glm <- glm(log(brain) ~ log(body), data = mammals)
# actually is a linear model regressing log(brain) on log(body) here
cv.err <- cv.glm(mammals, mammals.glm)$delta # leave-one-out CV
cv.err

cv.err.6 <- cv.glm(mammals, mammals.glm, K = 6)$delta # 6-fold CV
cv.err.6


