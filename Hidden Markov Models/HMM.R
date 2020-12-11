# Task 1 - Build a hidden Markov model (HMM)


# Defining the trasmissions and  the emissions probabilites matrices.

set.seed(12345)
library("HMM")
# creating the emissons probability matrix 
ep= diag(0.2,10)
ep[row(ep) == (col(ep) -1)] = 0.2
ep[row(ep) == (col(ep) -2)] = 0.2
ep[row(ep) == (col(ep) +1)] = 0.2
ep[row(ep) == (col(ep) +2)] = 0.2
ep[c(1,1,2,9,10,10),c(9,10,10,1,1,2)]=0.2
ep[9,2]=0;ep[2,9]=0
ep

# transProbs matrix
tp = diag(0.5,10)
tp[row(tp)==(col(tp)-1)]=0.5
tp[10,1]=0.5
tp

Sectors=paste0("S",1:10)
hmm = initHMM(Sectors,Sectors,
              startProbs = rep(0.1,10),transProbs = tp,
              emissionProbs = ep)

# Task 2 - Simulate the HMM for 100 time steps*

simulations <- simHMM(hmm, 100)


# Task 3 - Discard the hidden states from the sample obtained above.

# observations
simulations$observation

observation=simulations$observation
# filtering 
logfo = forward(hmm,observation)
loge = exp(logfo)

# smoothing
observation=simulations$observation
posterior = posterior(hmm,observation)

# the most probable path is 
mostpro = viterbi(hmm, simulations$observation)


# Task 4 - Compute the accuracy of the filtered and smoothed probability distributions, and of the most probable path.

#head(prop.table(exp(logfo)))
 forward = prop.table(loge,2)


# The accuracy for filtered probability distributions

#simulations$states
forpat = apply(forward,2,which.max)
forpath= rownames(forward)[forpat]
for_acc = length(which(forpath==simulations$states))/100
for_acc

# The accuracy for smoothed probability distributions

smoo = apply(posterior,MARGIN = 2,FUN=which.max)
smoopath= rownames(posterior)[smoo]
smoo_acc = length(which(smoopath==simulations$states))/100
smoo_acc

# The accuracy for most probable path:

mos_acc = length(which(mostpro==simulations$states))/100
mos_acc


# Smoothing gives us the highest accuracy and thats is due to the fact that it uses both previous and the future values to make a predictrion.




# Task 5  -Repeat the previous task with different simulated samples. In general, the smoothed distributions should be more accurate than the filtered distributions.

f100 = rep(0,100)
smoo100=rep(0,100)
mos100=rep(0,100)
for(i in 1:100){
simulations <- simHMM(hmm, 100)
observation=simulations$observation
logfo = forward(hmm,observation)
posterior = posterior(hmm,observation)
mostpro = viterbi(hmm, simulations$observation)
#forward acc
filtered = prop.table(exp(logfo),2)
forpat = apply(filtered, 2,which.max)
forwpath = rownames(logfo)[forpat]
#forwpath= table(forpat==simulations$states)
f100[i] = length(which(forwpath==simulations$states))/100
# post acc
smoo = apply(posterior,MARGIN = 2,FUN=which.max)
smoopath= rownames(posterior)[smoo]
smoo100[i] = length(which(smoopath==simulations$states))/100
# mo acc
mos100[i] = length(which(mostpro==simulations$states))/100
}

cat("The accuracy for filter is:",mean(f100))
cat("\nThe accuracy for smoothing is:",mean(smoo100))
cat("\nThe accuracy for path is:",mean(mos100))


# The reason for the smoothing to be the most accurate is that it uses all observations, including observations in the future, 
# which the filtering function doesn't. And viterbi has a logical constraint which could be useful but in this case will result in 
# less accuracy comparing to the smoothing.


# Task 6 - Is it true that the more observations you have the better you know where the robot is?

library(entropy)
ento=apply(prop.table(loge),2,entropy.empirical)
plot(ento,type = "l")

# As we can see in the plot the entropy does not converges to a stable value, which means having more data doesn't help us with accuracy.


# Task 7 - Consider any of the samples above of length 100. Compute the probabilities of the hidden states for the time step 101.

#filterProbs * trasmission to get 101
probablities = t(forward[,100])%*%tp
mostp = which.max(probablities)
mostp
probablities


# The probabilities of the hidden state for the time step 101 are: `r probablities` and the most probable is: `r mostp`


