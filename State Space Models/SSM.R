

# Task 1

set.seed(12345)
library(ggplot2)
ts = 100 #time steps 
sigma= 1
states = matrix(0,nrow=ts)
# Transition function
transition = function(mu,sigma=sigma){
  temp = sample(c(1,2,3),1)
  newmu= switch(temp,  # Different possible steps 
        z=mu,
        z=mu+1,
        z = mu+2)
  z = rnorm(1,newmu,sigma) 
  return(z)
}
states[1,1]=runif(1,min = 0, max = 100) # given uniform dist for the first step
# state simulation 
for ( i in 2:ts){
states[i,] = transition(mu=states[i-1,],sigma = 1)
}
plot(states,type = "l")

# Emission 
Emission = function(mu,sigma=sigma){
  temp = sample(c(1,2,3),1) # devided by 3 so equall probability 
                            # for any of them to happen 
  newmu= switch(temp,  # Different possible steps 
        z=mu,       
        z=mu-1,
        z = mu+1)
  x = rnorm(1,newmu,sigma) 
  return(x)
}
observations = apply(states,1,Emission,sigma=1)
plot(states,type = "l",col="red")
lines(observations,col="blue")
legend("topleft",c("red","blue"),c("states","observations"),lty = c(1,1),col = c("red","blue"))


# weight update function 
Wup = function(X,z,sd=1){
    W <- (1/3)*(dnorm(X, mean = z, sd = sd)+
    dnorm(X, mean = z-1, sd = sd)+
    dnorm(X, mean = z+1, sd = sd))
    return(W)
}
# Particle filter 
Particel = function(observations,ts=100,partc=100,sigma=1){
  pointest= rep(0,100)
  W = matrix(0,nrow=100,ncol=100)
x = matrix(0,nrow=100,ncol=100)
    x0 = runif(partc,0,100)
    w = rep(1/partc,partc)
    for(t in 1:ts){
    for(p in 1:partc){
    wtemp = sample(1:partc,1, replace=TRUE,prob = w)  # random choice based on weights
      temp = sample(c(1,2,3),1) 
      newmu= switch(temp,  # Different possible steps 
        z=x0[wtemp],       
        z=x0[wtemp]+1,
        z = x0[wtemp]-1)
  x[p,t] = rnorm(1,newmu,sigma)               # position 
  W[p,t] = Wup(X=x[p,t],z=observations[t],sigma)  #updating weight matrix
                      }
  w= W[,t]/sum(W[,t]) # normalizing the weights(should always sum up to 1)
  pointest[t] = sum(w*x[,t]) # point estimator for this step
  x0 = x[,t]
 }
 results = data.frame("point estimators"=pointest,"positions"=x,"weights"=W)
 results
}
Party =Particel(observations)
point_estimator = Party$point.estimators
Particels = Party[2:101]

plot(states,type="l",col="blue",ylim = c(-2,190))
lines(observations,col="red")
lines(x=1:100,point_estimator,col="green")
points(x=rep(1,100),y=Particels[,1])
points(x=rep(33,100),y=Particels[,33])
points(x=rep(66,100),y=Particels[,66])
points(x=rep(100,100),y=Particels[,100])
legend("topleft",c("blue","red","green","black"),c("States","Observations","point estimates","Particels"),col=c("blue","red","green","black"),lty=c(1,1,1,NA),bg="white",pch=c(NA,NA,NA,1))


# Task 2 


# standard devision for emission = 5
observations5 = apply(states,1,Emission,sigma=5)
Party5 = Particel(observations,sigma = 5)
#plot
point_estimator5 = Party5$point.estimators
Particels5 = Party5[2:101]
plot(states,type="l",col="blue")
lines(observations5,col="red")
lines(x=1:100,point_estimator5,col="green")
points(x=rep(1,100),y=Particels5[,1])
points(x=rep(33,100),y=Particels5[,33])
points(x=rep(66,100),y=Particels5[,66])
points(x=rep(100,100),y=Particels5[,100])
legend("topleft",c("blue","red","green","black"),c("States","Observations","point estimates","Particels"),col=c("blue","red","green","black"),lty=c(1,1,1,NA),bg="white",pch=c(NA,NA,NA,1))

# standard devision for emission = 5
observations50 = apply(states,1,Emission,sigma=50)
Party50 = Particel(observations,sigma = 50)
#plot
point_estimator50 = Party50$point.estimators
Particels50 = Party50[2:101]
plot(states,type="l",col="blue",ylim = c(-100,400))
lines(observations50,col="red")
lines(x=1:100,point_estimator50,col="green")
points(x=rep(1,100),y=Particels50[,1])
points(x=rep(33,100),y=Particels50[,33])
points(x=rep(66,100),y=Particels50[,66])
points(x=rep(100,100),y=Particels50[,100])
legend("topleft",c("blue","red","green","black"),c("States","Observations","point estimates","Particels"),col=c("blue","red","green","black"),lty=c(1,1,1,NA),bg="white",pch=c(NA,NA,NA,1))
```
# As we can see increasing the standard deviation of the emission will result in a greater range of error,
# however the mean of the observations (sensor readings) stays the same and close to actual state of the robot.

# task 3

set.seed(12345)
Particel2 = function(observations,ts=100,partc=100,sigma=1){
  pointest= rep(0,100)
  W = matrix(0,nrow=100,ncol=100)
x = matrix(0,nrow=100,ncol=100)
    x0 = runif(partc,0,100)
    w = rep(1,partc)
    for(t in 1:ts){
    for(p in 1:partc){
    wtemp = sample(1:partc,1, replace=TRUE,prob = w)  # random choice based on weights
      temp = sample(c(1,2,3),1) 
      newmu= switch(temp,  # Different possible steps 
        z=x0[wtemp],       
        z=x0[wtemp]+1,
        z = x0[wtemp]+2)
  x[p,t] = rnorm(1,newmu,sigma)               # position 
  W[p,t] = 1  #updating weight matrix
                      }
  w= W[,t]/sum(W[,t]) # normalizing the weights(should always sum up to 1)
  pointest[t] = sum(w*x[,t]) # point estimator for this step
  x0 = x[,t]
 }
 results = data.frame("point estimators"=pointest,"positions"=x,"weights"=W)
 results
}
weightlessParty =Particel2(observations)
point_estimator = weightlessParty$point.estimators
Particelsw = weightlessParty[,2:101]
# Plotting 
plot(states,type="l",col="blue",ylim = c(0,250))
lines(observations,col="red")
lines(x=1:100,point_estimator,col="green")
points(x=rep(1,100),y=Particelsw[,1])
points(x=rep(33,100),y=Particelsw[,33])
points(x=rep(66,100),y=Particelsw[,66])
points(x=rep(100,100),y=Particelsw[,100])
legend("topleft",c("States","Observations","point estimates","Particels"),col=c("blue","red","green","black"),lty=c(1,1,1,NA),bg="white",pch=c(NA,NA,NA,1))


#Having the same weight equal to 1 in the particle filter (no correction), will result in a point estimator with a close trend as the state but not as accurate. 

