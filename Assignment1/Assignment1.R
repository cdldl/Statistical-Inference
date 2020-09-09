install.packages("rmarkdown")
install.packages('HistData')
library(rmarkdown)
library(data.table)
library(HistData)
rm(list=ls())


# Baseline model
data(GaltonFamilies)
nrep <- 50
n <- dim(GaltonFamilies)[1]
ntest <- 200
MSE <- data.frame(matrix(0,nrep,1))
names(MSE) <- c("lm")
# partion training data and testing data 
for (i in 1:nrep){
  #i=1
  train = sample(1:n, n-ntest, replace = FALSE)
  traindata <- GaltonFamilies[train,]
  testdata <- GaltonFamilies[-train,]
  # Fit a baseline model
  fit_baseline <- lm(I(childHeight - mean(childHeight))~gender + I(midparentHeight -mean(midparentHeight))-1, traindata)
  # make prediction based on the fitted model
  height_baseline_pred <- predict(fit_baseline,testdata) + mean(traindata$childHeight)
  # Evaluate model performance
  MSE[i,1] <- mean((testdata$childHeight - height_baseline_pred)^2)
  # Warning: You should not use childNum as a predictor in your model
  # because children within a family are listed in decreasing order of height
  # for boys followed by girls
}
boxplot(MSE)
mean(MSE[,1])


Second approach with James Stein estimator

js = function(draw) {
  l2 <- sum(draw^2)
  return((1 - (length(draw) - 2) / l2) *draw)
}

MSE <- data.frame(matrix(0,nrep,1))
names(MSE) <- c("lm")
for (i in 1:nrep){
  #i=1
  # partion training data and testing data
  train = sample(1:n, n-ntest, replace = FALSE)
  traindata <- GaltonFamilies[train,]
  testdata <- GaltonFamilies[-train,]
  
  # Fit a JS model
  fit_js =  lm(js(childHeight)~gender + js(midparentHeight), traindata)
  #fit_baseline <- lm(childHeight~gender + I(midparentHeight -js(midparentHeight))-1, traindata)
  # make prediction based on the fitted model
  height_baseline_pred <- predict(fit_js,testdata) #+ (1 - 1/(sd(traindata$childHeight)^2+1))*mean(traindata$childHeight)
  # Evaluate model performance
  MSE[i,1] <- mean((testdata$childHeight - height_baseline_pred)^2)
  # Warning: You should not use childNum as a predictor in your model
  # because children within a family are listed in decreasing order of height
  # for boys followed by girls
}
boxplot(MSE)
mean(MSE[,1])

# Third approach:  Exponential GLM following Geometric mean regression
MSE <- data.frame(matrix(0,nrep,1))
names(MSE) <- c("lm")
for (i in 1:nrep){
  train = sample(1:n, n-ntest, replace = FALSE)
  # partion training data and testing data
  traindata <- GaltonFamilies[train,]
  testdata <- GaltonFamilies[-train,]
  # Fit a model
  fit_gme <- glm(js(childHeight)~gender + js(midparentHeight),
                      family = Gamma(link="log"),data=traindata)
  # make prediction based on the fitted model
  height_baseline_pred <- predict(fit_gme,testdata,type='response')
  # Evaluate model performance
  MSE[i,1] <- mean((testdata$childHeight - height_baseline_pred)^2)
  # Warning: You should not use childNum as a predictor in your model
  # because children within a family are listed in decreasing order of height
  # for boys followed by girls
}
boxplot(MSE)
mean(MSE[,1])


# Fourth approach: Tweedie
library(cplm)
library(tweedie)
library(statmod)
xi.vec <- seq(1, 3, by=0.2)
out <- tweedie.profile(GaltonFamilies$childHeight ~1, xi.vec=xi.vec, do.plot=TRUE, verbose=TRUE)
MSE <- data.frame(matrix(0,nrep,1))
names(MSE) <- c("lm")
for (i in 1:nrep){
  #i=1
  # partion training data and testing data
  train = sample(1:n, n-ntest, replace = FALSE)
  traindata <- GaltonFamilies[train,]
  testdata <- GaltonFamilies[-train,]
  # Fit a model
  fit_baseline <- glm(js(childHeight)~gender + js(midparentHeight),
                      data =traindata, 
                      family=tweedie(var.power=out$xi.max,link.power=0))
  # make prediction based on the fitted model
  height_baseline_pred <- predict(fit_baseline,testdata,type='response')
  # Evaluate model performance
  MSE[i,1] <- mean((testdata$childHeight - height_baseline_pred)^2)
}
boxplot(MSE)
mean(MSE[,1])


# Fifth Approach: Stan
library(rstanarm)
glm_mcmc = stan_glm(js(childHeight) ~gender + js(midparentHeight) -1,data=GaltonFamilies, family=Gamma(link="log"),chains=2)
summary(glm_mcmc)
nrep <- 50
n <- dim(GaltonFamilies)[1]
ntest <- 200
MSE <- data.frame(matrix(0,nrep,1))
names(MSE) <- c("lm")
for (i in 1:nrep){
  #i=1
  # partion training data and testing data
  train = sample(1:n, n-ntest, replace = FALSE)
  traindata <- GaltonFamilies[train,]
  testdata <- GaltonFamilies[-train,]
  for(j in 1:nrow(testdata)) testdata$midparentHeight[j] = js(c(traindata$midparentHeight,testdata$midparentHeight[j]))[length(traindata$midparentHeight)+1]
  
  # make prediction based on the fitted model
  height_baseline_pred <- posterior_predict(glm_mcmc,testdata)
  # Evaluate model performance
  MSE[i,1] <- mean((testdata$childHeight - apply(height_baseline_pred,2,mean))^2)
}
par(mfrow=c(2,2))
b = boxplot(MSE,plot=F)
mean(MSE[,1])

