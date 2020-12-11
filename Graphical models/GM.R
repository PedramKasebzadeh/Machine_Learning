
#1. Show that multiple runs of the hill-climbing algorithm can return non-equivalent Bayesian network (BN) structures. Explain why this happens.

set.seed(12345)
library(bnlearn)
data("asia")
res1 <- hc(asia, score="bde", iss = 30, restart = 4)
#res1$arcs
#vstructs(res1, arcs = FALSE, moral = FALSE, debug = FALSE)
cres1 <- cpdag(res1)
plot(cres1)
res2 <- hc(asia, score="bde", iss = 30, restart = 4)
#res2$arcs
#vstructs(res2, arcs = FALSE, moral = FALSE, debug = FALSE)
cres2 <- cpdag(res2)
plot(cres2)
all.equal(res1,res2)

# We can see two bayesian network above, both created by the same code that look differently.
# This is cause because the hill climbing algorithm does not gueratee us the best result. 
# The idea of hill climbing algorithm is in maximizing the probability given the data we have. 
# Unfortunately hc guarantees us just the local maximum and that is why we can get two different results with the same code. 
# On top of that adding imaginary sample influence the relationships between the nodes so that we get more edges. 
# Number of restarts however increases the possibility of getting the same result but with such a high number of imaginary sample we would need much higher number of restarts.

# 2. Learn both the structure and the parameters of a BN

# Below is the code I used to learn both the parameters and the structure of the BN. 
# Note that the reason for me using loglik as score is only due to the fact that when leaving it blank, 
# the BN get the same confusion matrix. In the begining this made me very confused since I did not expect the network to do as good as the "true" network.

library(bnlearn)
library(gRain)
set.seed(12345)
data("asia")
n = nrow(asia)*0.8
id = sample(nrow(asia),n)
train = asia[id,]
test = asia[-id,]
test <- data.frame(lapply(test, as.character), stringsAsFactors=FALSE)
exact  = hc(train)
plot(exact)
exact = bn.fit(exact,train)
exact = compile(as.grain(exact))
Prediction = rep(0,nrow(test))
for(i in 1:nrow(test)){
find=setEvidence(exact,nodes = colnames(test[,-2]) ,states = test[i,-2] )
q = querygrain(find)
Prediction[i] = names(which.max(q$S))
}
dag = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
plot(dag)
true = bn.fit(dag,train)
true_BN = compile(as.grain(true))
Prediction_true = rep(0,nrow(test))
for(i in 1:nrow(test)){
find_t=setEvidence(true_BN,nodes = colnames(test[,-2]) ,states = test[i,-2] )
q_t = querygrain(find_t)
Prediction_true[i] = names(which.max(q_t$S))
}
real = test[,2]
table(Prediction,real)
table(Prediction_true,real)



# We can see that even with hc which guarantees us just a local optimum we achieved the same result as with the true Asia BN. 
# This is because of the moralization and triangulation.
# The only problem is that with not optimal graph it might take more time but in case of a small network as this it does not matter. 



# 3. classify *S* given observations only for the so-called Markov blanket of *S*

###  the Markov blanket of S, which is a set
hcbn_train <- hc(x=train, score = "bde", iss = 3)
hcbn_train_fit <- bn.fit(x=hcbn_train, data=train)
test2 <- data.frame(lapply(test, as.character), stringsAsFactors=FALSE)
names(test2)[-2]
hcbn_train_grain <- as.grain(hcbn_train_fit)
hcbn_train_grain <- compile(hcbn_train_grain)
hcbn_train <- hc(x=train, score = "bde", iss = 3)
hcbn_train_fit <- bn.fit(x=hcbn_train, data=train)
S_markovblanket <- mb(hcbn_train_fit, "S")
predictS_mb <- rep(0, nrow(test2) )
for(i in 1:nrow(test2) ){
  finding_mb <- setFinding(hcbn_train_grain, nodes = S_markovblanket, 
                        states = test2[i, S_markovblanket])
  query_mb <- querygrain(finding_mb)
  predictS_mb[i] <- names(which.max(query_mb$S))
}
table(predictS_mb,test2$S)

# We still get the same confusion matrix with the "Markov blanket of S" as with the true directed acyclic graph model and with the hill climbing algorithm.
# This is because only the parents and the children in case there is no collider are required in obtaining the conditional probabilities as shown by the "Markov blanket of S".


# 4. Using a naive Bayes classifier

#2, Train on 80%
set.seed(12345)
n=dim(asia)[1]
id=sample(1:n, floor(n*0.8))
train=asia[id,]
y_test=asia[-id,]$S
test = subset(asia[-id,], select = -c(S))
# pred function from Jespers task 2 code
pred = function(model, state, target = "S", nodes = colnames(test)){
  #Give the observation/evidence
  e = setEvidence(object = model, states = state, nodes = nodes)
  #Get the probability for the different nodes given the evidence
  q = querygrain(e, target)
  #Return the state which has the highest probability
  return(names(which.max(q$S)))
}
#4, Naive Bayes
#The NB graph is S pointing to every other varible
#Using the code from the examples in bnlearn we get:
NB = empty.graph(colnames(asia))
arc.set = matrix(c("S", "A",
                   "S", "T",
                   "S", "L",
                   "S", "B",
                   "S", "E",
                   "S", "X",
                   "S", "D"),
                 ncol = 2, byrow = T,
                 dimnames = list(NULL, c("from", "to")))
arcs(NB) = arc.set
plot(NB)
NB = bn.fit(NB, data = train)
NB = compile(as.grain(NB))
predictions=c()
for (i in 1:nrow(test)) {
  predictions=c(predictions, pred(NB, state = as.character(unlist(test[i,]))))
}
table(predictions, y_test)

