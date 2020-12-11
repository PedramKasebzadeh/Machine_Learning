# Task 1 - Implementing GP Regression

library("mvtnorm")
library("kernlab")
library("ggplot2")
# Squared exponential kernel function given as a hint
# Varying the Hyperparameters
SquaredExp <- function(x1,x2,sigmaF,l){
  n1 <- length(x1)
  n2 <- length(x2)
  K <- matrix(NA,n1,n2)
  for (i in 1:n2){
    K[,i] <- sigmaF^2*exp(((-0.5)/(l^2))* (x1-x2[i])^2 )
  }
  return(K)
}

# Gaussian process regression  of Rasmussen's book
PosteriorGP <- function(X,y,Xstar,hyperParam,sigmaNoise){
  # K(x,x)
  k <- SquaredExp(X,X,hyperParam[1],hyperParam[2])
  L <- t(chol(k + diag(sigmaNoise^2,length(X)))) #the R function returns an upper triangular
  # matrix. So, we need to transpose the output
  
  #predictive mean eq
  alphatemp <- solve(L,y)
  alpha <- solve(t(L),alphatemp)
    # K(x,xstar)
    #predictive variance
    kstar <- SquaredExp(X,Xstar,hyperParam[1],hyperParam[2])
    fstar<- t(kstar)%*%alpha
    v <- solve(L,kstar)
    # K(xstar,xstar)
    kstarstar <- SquaredExp(Xstar,Xstar,hyperParam[1],hyperParam[2])
    var <- kstarstar - t(v)%*%v
    out <- list("mean" = fstar, "var" = var)
    return(out)
}


# Task 2 - update the prior 

#simulating from the posterior distribution
# given values 
sigmaf <- 1
l <- 0.3
hyperParam<-c(sigmaf,l)
sigmaNoise <- 0.1
X <- 0.4
y <- 0.719
Xstar <- seq(-1,1,length=20)
res = PosteriorGP(X,y,Xstar,hyperParam,sigmaNoise)
#95% probability band
lbound = c()
ubound = c()
for (i in 1:length(res$mean)){
    lbound[i] <- res$mean[i] - 1.96*sqrt(res$var[i,i])
    ubound[i] <- res$mean[i] + 1.96*sqrt(res$var[i,i])
  }
# ploting 
plot(x=Xstar,y=res$mean,type="l",ylim=c(-3,3))
# Thanks to https://stat.ethz.ch/pipermail/r-help/2008-July/168230.html
polygon(c(Xstar,rev(Xstar)),c(ubound,rev(lbound)), col = "lightgray" ,border = NA)
points(X,y,col="Red")
lines(x=Xstar,res$mean,col="green")
legend("bottom",c("posterior mean", "Measurement"),col=c("green","red"),lty=c(1,NA),pch=c(NA,1))


# Task 3 - Update the posterior


#simulating from the posterior distribution (Two observations)
# given values 
sigmaf <- 1
l <- 0.3
hyperParam<-c(sigmaf,l)
sigmaNoise <- 0.1
X <- c(-0.6,0.4)
y <- c(-0.044,0.719)
Xstar <- seq(-1,1,length=20)
res = PosteriorGP(X,y,Xstar,hyperParam,sigmaNoise)
#95% probability band
lbound = c()
ubound = c()
for (i in 1:length(res$mean)){
    lbound[i] <- res$mean[i] - 1.96*sqrt(res$var[i,i])
    ubound[i] <- res$mean[i] + 1.96*sqrt(res$var[i,i])
  }
# ploting 
plot(x=Xstar,y=res$mean,type="l",ylim=c(-3,3))
# Thanks to https://stat.ethz.ch/pipermail/r-help/2008-July/168230.html
polygon(c(Xstar,rev(Xstar)),c(ubound,rev(lbound)), col = "lightgray" ,border = NA)
points(X,y,col="Red")
lines(x=Xstar,res$mean,col="green")
legend("bottom",c("posterior mean", "Measurement"),col=c("green","red"),lty=c(1,NA),pch=c(NA,1))


# Task 4 - Copmute the posterior distribution 

#simulating from the posterior distribution (five observations)
# given values 
sigmaf <- 1
l <- 0.3
hyperParam<-c(sigmaf,l)
sigmaNoise <- 0.1
X <- c(-1,-0.6,-0.2,0.4,0.8)
y <-c(0.768,-0.044,-0.940,0.719,-0.664)
Xstar <- seq(-1,1,length=20)
res = PosteriorGP(X,y,Xstar,hyperParam,sigmaNoise)
#95% probability band
lbound = c()
ubound = c()
for (i in 1:length(res$mean)){
    lbound[i] <- res$mean[i] - 1.96*sqrt(res$var[i,i])
    ubound[i] <- res$mean[i] + 1.96*sqrt(res$var[i,i])
  }
# ploting 
plot(x=Xstar,y=res$mean,type="l",ylim=c(-1.5,1.5))
# Thanks to https://stat.ethz.ch/pipermail/r-help/2008-July/168230.html
polygon(c(Xstar,rev(Xstar)),c(ubound,rev(lbound)), col = "lightgray" ,border = NA)
points(X,y,col="Red")
lines(x=Xstar,res$mean,col="green")
legend("topright",c("posterior mean", "Measurement"),col=c("green","red"),lty=c(1,NA),pch=c(NA,1))


# Task 5 - Task 4 with new parameters 


# Five observations and new hyper parameters
sigmaf <- 1
l <- 1
hyperParam<-c(sigmaf,l)
sigmaNoise <- 0.1
X <- c(-1,-0.6,-0.2,0.4,0.8)
y <-c(0.768,-0.044,-0.940,0.719,-0.664)
Xstar <- seq(-1,1,length=20)
res = PosteriorGP(X,y,Xstar,hyperParam,sigmaNoise)
#95% probability band
lbound = c()
ubound = c()
for (i in 1:length(res$mean)){
    lbound[i] <- res$mean[i] - 1.96*sqrt(res$var[i,i])
    ubound[i] <- res$mean[i] + 1.96*sqrt(res$var[i,i])
  }
# ploting 
plot(x=Xstar,y=res$mean,type="l",ylim=c(-1.5,1.5))
# Thanks to https://stat.ethz.ch/pipermail/r-help/2008-July/168230.html
polygon(c(Xstar,rev(Xstar)),c(ubound,rev(lbound)), col = "lightgray" ,border = NA)
points(X,y,col="Red")
lines(x=Xstar,res$mean,col="green")
legend("topright",c("posterior mean", "Measurement"),col=c("green","red"),lty=c(1,NA),pch=c(NA,1))

# Task 6 - Stockholm weather 

# preparing data 
data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/TempTullinge.csv", header=TRUE, sep=";")
date <- data$date
temp <- data$temp
time <- seq(1,2190,5)
day <- rep(seq(1,365,5),6)

# Task 6.1
# The squared exponential kernel function
SqExKernel <- function(sigmaf, ell)
{
  rval <- function(x, y = NULL) {
    r = sqrt(crossprod(x-y));
    return(sigmaf^2*exp(-0.5*( (r)/ell)^2 ))
  }
  class(rval) <- "kernel"
  return(rval)
}
SEKernelFunc <- SqExKernel(sigmaf = 1, ell = 1)
# SEKernelFunc(1,2) # evaluate kernel at (1,2)
K <- kernelMatrix(kernel = SEKernelFunc, x = matrix(c(1,3,4), 3, 1), y = matrix(c(2,3,4), 3, 1))
K


# Task 6.2
# Estimating the noise variance from a quadratic regression
tempe = temp[time]
fit <- lm(tempe ~  time + I(time^2) )
sigmaNoise = sd(fit$residuals)
# Fit the GP with the defined kernel function
sigmaf <- 20
ell <- 0.2
GPfit <- gausspr(time, tempe, kernel = SqExKernel(sigmaf,ell), var = sigmaNoise^2)
meanPred <- predict(GPfit, time)
# plotting 
plot(x = time, y = tempe)
lines(x = time, y = meanPred, col = "red")
legend("bottomright",c("Temperatures","Prediction"),pch = c(1,NA),lty = c(NA,1),
       col=c("black","red"))
Temperature has a seasonality so it makes sense to have a good prediction.

# Task 6.3
tempScaled <- scale(tempe)
timeScaled <- scale(time)
xGrid <- seq(1,2190,by=5)
xGridScaled <- scale(xGrid)
hyperParam<-c(sigmaf,ell)
est <- PosteriorGP(timeScaled,tempScaled,xGridScaled,hyperParam,sigmaNoise)
# We need to undo the effect of scale to get the correct result
posteriorMean_ <- est$mean
posteriorMean <- posteriorMean_ * attr(tempScaled, 'scaled:scale') + attr(tempScaled, 'scaled:center')
lbound <- matrix(NA,length(posteriorMean),1)
ubound <- matrix(nrow = length(posteriorMean), ncol = 1 ,0)
for (i in 1:length(est$mean)){
  var <- est$var[i,i]
  var <- var * attr(tempScaled, 'scaled:scale') + attr(tempScaled, 'scaled:center')
  lbound[i] <- posteriorMean[i] - 1.96*sqrt(var)
  ubound[i] <- posteriorMean[i] + 1.96*sqrt(var)
}
plot(time, tempe,ylim=c(-30,28),main="95% probability")
polygon(c(time, rev(time)),c(ubound, rev(lbound)),col = "gray")
points(x = time, y = tempe,pch=16,col="black")
lines(x = time, y = meanPred, col = "red")
legend("bottomright",c("Temp","Prediction"),col=c("black","red"),
       lty=c(NA,1),pch=c(16,NA))


# Task 6.4
# Estimating the noise variance from a quadratic regression
fit2 <- lm(tempe ~  day + I(day^2) )
sigmaNoise = sd(fit2$residuals)
GDfit <- gausspr(day, tempe, kernel = SqExKernel(sigmaf,ell), var = sigmaNoise^2)
daypredict <- predict(GDfit, day)
lineWidth = 0.6
dfplot<-as.data.frame(xGrid)
colnames(dfplot)<-"x"
dfplot$lbound<-lbound
dfplot$ubound<-ubound
plot(time, meanPred,ylim=c(-30,28),main="95% probability")
polygon(c(dfplot$x, rev(dfplot$x)),c(ubound, rev(lbound)),col = "gray")
points(x = time, y = tempe,pch=16,col="black")
lines(x = time, y = daypredict, col = "red")
lines(x = time, y = meanPred, col = "green")
legend("bottomright",c("Temp","Day","time"),col=c("black","red","green"),
       lty=c(NA,1,1),pch=c(16,NA,NA))

# Task 6.5
sigmaf <- 20
ell <- c(1,10)
quadFit <- lm(tempe ~  time + I(time^2) )
sigmaNoise = sd(quadFit$residuals)
PeriodicKernel <- function(sigmaf, ell1,ell2)
{
  rval <- function(x, y = NULL) {
    r = sqrt(crossprod(x-y));
    d = 365/sd(time)
    return( sigmaf^2 * exp(-0.5*( (r)/ell[2])^2 ) * exp(-2*( (sin(pi*r/d))/ell[1])^2) )
  }
  class(rval) <- "kernel"
  return(rval)
}
GPfitPeriodic <- gausspr(time, tempe, kernel = PeriodicKernel(sigmaf,ell[1],ell[2]), var = sigmaNoise^2)
meanPredPeriodic  <- predict(GPfitPeriodic, time)
plot(x=time,y=meanPredPeriodic,type="l",col="red",ylim=c(-20,25))
points(x=time,y=tempe)
# Training errors
errSETime = GPfit@error
errSEDay = GDfit@error
errPrTime = GPfitPeriodic@error


# Task 7 - 



data <- read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")
names(data) <- c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] <- as.factor(data[,5])
set.seed(111);
SelectTraining <- sample(1:dim(data)[1], size = 1000,replace = FALSE)
train <- data[SelectTraining,]
test <- data[-SelectTraining,]
# Task 3.1
GPfitbanknote <- gausspr(fraud ~  varWave + skewWave, data=train)
GPfitbanknote
library(AtmRay)
x1 <- seq(min(train[,1])-1,max(train[,1])+1,length=200)
x2 <- seq(min(train[,2])-1,max(train[,2])+1,length=200)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))
gridPoints <- data.frame(gridPoints)
names(gridPoints) <- names(train)[1:2]
probPreds <- predict(GPfitbanknote, gridPoints, type="probabilities")
# Confusion matrix and accuracy on the training set
predVarSkewTrain <- predict(GPfitbanknote,train[,1:2])
confMatVarSkewTrain <- table(predVarSkewTrain, train[,5]) # confusion matrix
#confMatVarSkewTrain 
accVarSkewTrain <- sum(diag(confMatVarSkewTrain))/sum(confMatVarSkewTrain) 
# Misclassification rate
accVarSkewTrain
# Plotting for 
contour(x1,x2,matrix(probPreds[,1],200,byrow = TRUE), 20, xlab = "varWave", ylab = "skewWave")
points(train[train[,5]==1,1],train[train[,5]==1,2],col="blue")
points(train[train[,5]==0,1],train[train[,5]==0,2],col="red")


# Task 7.2
# Confusion matrix and accuracy on the test set
predVarSkewTest <- predict(GPfitbanknote,test[,1:2])
confMatVarSkewTest <- table(predVarSkewTest, test[,5]) # confusion matrix
confMatVarSkewTest
accVarSkewTest <- sum(diag(confMatVarSkewTest))/sum(confMatVarSkewTest) # Misclassification rate
accVarSkewTest

# Task 7.3
# Classification using all four covariates
GPfitbanknote <- gausspr(fraud ~  varWave + skewWave + kurtWave + entropyWave, data=train)
GPfitbanknote
# Confusion matrix and accuracy on the training set
predAllTrain <- predict(GPfitbanknote,train[,1:4])
confMatAllTrain <- table(predAllTrain, train[,5]) # confusion matrix
confMatAllTrain 
accAllTrain <- sum(diag(confMatAllTrain))/sum(confMatAllTrain) # Misclassification rate
accAllTrain
# Confusion matrix and accuracy on the test set
predAllTest <- predict(GPfitbanknote,test[,1:4])
confMatAllTest <- table(predAllTest, test[,5]) # confusion matrix
confMatAllTest
accAllTest <- sum(diag(confMatAllTest))/sum(confMatAllTest) # Misclassification rate
accAllTest
