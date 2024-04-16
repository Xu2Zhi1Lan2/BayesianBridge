setwd("D:/桌面/Research Project/Bayesian bridge")
library(BayesBridge)
library(horseshoe)
library(glmnet)
library(MASS)
library(boot)
library(rbridge)
source("BBLR.R")
set.seed(20240313)


# convergence
## example 2 (example 1 in test_example.R)
### p > n
beta <- c(rep(5,10),rep(0,10))
sigma <- c(3,1)
n <- c(10,20)
SigA <- diag(rep(2.05,n[1]))+matrix(rep(0.95,n[1]^2),n[1],n[1])
A <- mvrnorm(20,rep(0,n[1]),SigA)
SigB <- diag(rep(0.05,n[1]))+matrix(rep(0.95,n[1]^2),n[1],n[1])
B <- mvrnorm(20,rep(0,n[1]),SigB)
yA <- as.vector(t(A)%*%beta + rnorm(n[1],0,sigma[1]))
yB <- as.vector(t(B)%*%beta + rnorm(n[1],0,sigma[2]))

BBLR.A <- BBLR(yA, t(A), tuning=1000,Tb=5,burn=30000,nmc=30000,thin=15,method.alpha="beta")
BBLR.B <- BBLR(yB, t(B), tuning=1000,Tb=5,burn=30000,nmc=30000,thin=15,method.alpha="beta")
## converge

# RMSE
## example 3: low dimension
set.seed(20240313)
library(mlbench)
data(Ozone)
# missing data processing by KNN, K=3
## missing data detection
library(mice)
library(VIM)
aggr(Ozone, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(Ozone), cex.axis=.7, gap=3, ylab=c("Missing data","Pattern"))
## imputation by KNN, K=3
Ozone.imp <- kNN(Ozone, k=3)
## data: response: V13; predictors: V1-V12
y <- Ozone.imp$V13
X <- Ozone.imp[,1:12]
X$V1 <- as.numeric(X$V1)
X$V2 <- as.numeric(X$V2)
X$V3 <- as.numeric(X$V3)
X <- as.matrix(X)
train <- sample(1:nrow(Ozone.imp),round(nrow(Ozone.imp)*0.5))
y.train <- y[train]
X.train <- X[train,]
y.test <- y[-train]
X.test <- X[-train,]

## lasso with bootstrap
B <- 2000
rmse1.boot <- rep(0,B)
for (i in 1:B){
  boot.index <- sample(1:nrow(X.train),nrow(X.train),replace=TRUE)
  X.boot <- X.train[boot.index,]
  y.boot <- y.train[boot.index]
  cv.ex1.lasso.boot <- cv.glmnet(X.boot, y.boot, alpha=1)
  lambda1.boot <- cv.ex1.lasso.boot$lambda.min
  ex1.lasso.boot <- glmnet(X.boot, y.boot, alpha=1, lambda=lambda1.boot)
  yhat1.boot <- predict(ex1.lasso.boot, s=lambda1.boot, newx=X.test)
  rmse1.boot[i] <- sqrt(mean((yhat1.boot-y.test)^2))
}
rmse1.mean <- mean(rmse1.boot)
rmse1.sd <- sd(rmse1.boot)
cat("RMSE of lasso with bootstrap: mean=",rmse1.mean,"; sd=",rmse1.sd,"\n")
boxplot(rmse1.boot, main="RMSE of lasso with bootstrap", ylab="RMSE")

## elastic net with bootstrap
rmse2.boot <- rep(0,B)
for (i in 1:B){
  boot.index <- sample(1:nrow(X.train),nrow(X.train),replace=TRUE)
  X.boot <- X.train[boot.index,]
  y.boot <- y.train[boot.index]
  cv.ex1.elastic.boot <- cv.glmnet(X.boot, y.boot, alpha=0.5)
  lambda2.boot <- cv.ex1.elastic.boot$lambda.min
  ex1.elastic.boot <- glmnet(X.boot, y.boot, alpha=0.5, lambda=lambda2.boot)
  yhat2.boot <- predict(ex1.elastic.boot, s=lambda2.boot, newx=X.test)
  rmse2.boot[i] <- sqrt(mean((yhat2.boot-y.test)^2))
}
rmse2.mean <- mean(rmse2.boot)
rmse2.sd <- sd(rmse2.boot)
cat("RMSE of elastic net with bootstrap: mean=",rmse2.mean,"; sd=",rmse2.sd,"\n")
boxplot(rmse2.boot, main="RMSE of elastic net with bootstrap", ylab="RMSE")

## horseshoe
burn_in <- 30000
nmc <- 30000
thin <- 15
ex3.horseshoe <- horseshoe(y.train,
                           X.train,
                           method.tau = "halfCauchy",
                           method.sigma = "Jeffreys",
                           Sigma2 = 1,
                           burn = burn_in,
                           nmc = nmc,
                           thin = thin)
rmse3 <- rep(0,nmc/thin)
for (i in 1:nmc/thin){
  beta.horseshoe <- ex3.horseshoe$BetaSamples[,i]
  yhat3 <- X.test%*%beta.horseshoe
  rmse3[i] <- sqrt(mean((yhat3-y.test)^2))
}
rmse3.mean <- mean(rmse3)
rmse3.sd <- sd(rmse3)
cat("RMSE of horseshoe: mean=",rmse3.mean,"; sd=",rmse3.sd,"\n")
boxplot(rmse3, main="RMSE of horseshoe", ylab="RMSE")

## BayesBridge
ex3.BBR <- bridge.reg.tri(y.train, X.train, nsamp=2000, alpha=0.5,
                          sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0, burn = burn_in)
for (i in 1:2000){
  beta.BBR <- ex3.BBR$beta[i,]
  yhat4 <- X.test%*%beta.BBR
  rmse4[i] <- sqrt(mean((yhat4-y.test)^2))
}
rmse4.mean <- mean(rmse4)
rmse4.sd <- sd(rmse4)
cat("RMSE of BayesBridge: mean=",rmse4.mean,"; sd=",rmse4.sd,"\n")
boxplot(rmse4, main="RMSE of BayesBridge", ylab="RMSE")

## BBLR
ex3.BBLR <- BBLR(y.train, X.train, tuning=20,Tb=2,burn=30000,nmc=30000,thin=15,method.alpha="fixed")
rmse5 <- rep(0,nmc/thin)
for (i in 1:nmc/thin){
  beta.BBLR <- ex3.BBLR$BetaSamples[,i]
  yhat5 <- X.test%*%beta.BBLR
  rmse5[i] <- sqrt(mean((yhat5-y.test)^2))
}
rmse5.mean <- mean(rmse5)
rmse5.sd <- sd(rmse5)
cat("RMSE of BBLR: mean=",rmse5.mean,"; sd=",rmse5.sd,"\n")
boxplot(rmse5, main="RMSE of BBLR", ylab="RMSE")


rmse <- cbind(rmse1.boot,rmse2.boot,rmse3,rmse4,rmse5)
boxplot(rmse, main="RMSE of different methods", ylab="RMSE", names=c("lasso","elastic net","horseshoe","BayesB","BBLR"))


## example 4: high dimension
set.seed(20240313)
library(chemometrics)
data(NIR)
y <- as.vector(NIR$yGlcEtOH$Glucose)
X <- as.matrix(NIR$xNIR)

train <- sample(1:nrow(X),round(nrow(X)*0.5))
y.train <- y[train]
X.train <- X[train,]
y.test <- y[-train]
X.test <- X[-train,]

## lasso with bootstrap
B <- 2000
rmse1.boot <- rep(0,B)
for (i in 1:B){
  print(i)
  boot.index <- sample(1:nrow(X.train),nrow(X.train),replace=TRUE)
  X.boot <- X.train[boot.index,]
  y.boot <- y.train[boot.index]
  cv.ex1.lasso.boot <- cv.glmnet(X.boot, y.boot, alpha=1)
  lambda1.boot <- cv.ex1.lasso.boot$lambda.min
  ex1.lasso.boot <- glmnet(X.boot, y.boot, alpha=1, lambda=lambda1.boot)
  yhat1.boot <- predict(ex1.lasso.boot, s=lambda1.boot, newx=X.test)
  rmse1.boot[i] <- sqrt(mean((yhat1.boot-y.test)^2))
}
rmse1.mean <- mean(rmse1.boot)
rmse1.sd <- sd(rmse1.boot)
cat("RMSE of lasso with bootstrap: mean=",rmse1.mean,"; sd=",rmse1.sd,"\n")
boxplot(rmse1.boot, main="RMSE of lasso with bootstrap", ylab="RMSE")

## elastic net with bootstrap
rmse2.boot <- rep(0,B)
for (i in 1:B){
  print(i)
  boot.index <- sample(1:nrow(X.train),nrow(X.train),replace=TRUE)
  X.boot <- X.train[boot.index,]
  y_boot <- y.train[boot.index]
  cv.ex1.elastic.boot <- cv.glmnet(X.boot, y_boot, alpha=0.5)
  lambda2.boot <- cv.ex1.elastic.boot$lambda.min
  ex1.elastic.boot <- glmnet(X.boot, y_boot, alpha=0.5, lambda=lambda2.boot)
  yhat2.boot <- predict(ex1.elastic.boot, s=lambda2.boot, newx=X.test)
  rmse2.boot[i] <- sqrt(mean((yhat2.boot-y.test)^2))
}
rmse2.mean <- mean(rmse2.boot)
rmse2.sd <- sd(rmse2.boot)
cat("RMSE of elastic net with bootstrap: mean=",rmse2.mean,"; sd=",rmse2.sd,"\n")
boxplot(rmse2.boot, main="RMSE of elastic net with bootstrap", ylab="RMSE")

## horseshoe
burn_in <- 30000
nmc <- 30000
thin <- 15
ex3.horseshoe <- horseshoe(y.train,
                           X.train,
                           method.tau = "halfCauchy",
                           method.sigma = "Jeffreys",
                           Sigma2 = 1,
                           burn = burn_in,
                           nmc = nmc,
                           thin = thin)
rmse3 <- rep(0,nmc/thin)
for (i in 1:nmc/thin){
  beta.horseshoe <- ex3.horseshoe$BetaSamples[,i]
  yhat3 <- X.test%*%beta.horseshoe
  rmse3[i] <- sqrt(mean((yhat3-y.test)^2))
}
rmse3.mean <- mean(rmse3)
rmse3.sd <- sd(rmse3)
cat("RMSE of horseshoe: mean=",rmse3.mean,"; sd=",rmse3.sd,"\n")
boxplot(rmse3, main="RMSE of horseshoe", ylab="RMSE")

## BayesB  (not converge)
ex3.BayesB <- bridge.reg.tri(y.train, X.train, nsamp=2000, alpha=0.5,
                             sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0, burn = burn_in)
rmse4 <- rep(0,nmc/thin)
for (i in 1:nmc/thin){
  beta.BayesB <- ex3.BayesB$beta[i,]
  yhat4 <- X.test%*%beta.BayesB
  rmse4[i] <- sqrt(mean((yhat4-y.test)^2))
}
rmse4.mean <- mean(rmse4)
rmse4.sd <- sd(rmse4)
cat("RMSE of BayesB: mean=",rmse4.mean,"; sd=",rmse4.sd,"\n")
boxplot(rmse4, main="RMSE of BayesB", ylab="RMSE")


## BBLR
ex3.BBLR <- BBLR(y.train, X.train, tuning=20,Tb=2,burn=30000,nmc=30000,thin=15,method.alpha="fixed")
rmse5 <- rep(0,nmc/thin)
for (i in 1:nmc/thin){
  beta.BBLR <- ex3.BBLR$BetaSamples[,i]
  yhat5 <- X.test%*%beta.BBLR
  rmse5[i] <- sqrt(mean((yhat5-y.test)^2))
}
rmse5.mean <- mean(rmse5)
rmse5.sd <- sd(rmse5)
cat("RMSE of BBLR: mean=",rmse5.mean,"; sd=",rmse5.sd,"\n")
boxplot(rmse5, main="RMSE of BBLR", ylab="RMSE")


rmse <- cbind(rmse1.boot,rmse2.boot,rmse3,rmse5)
boxplot(rmse, main="RMSE of different methods", ylab="RMSE", names=c("lasso","elastic net","horseshoe","BBLR"))


