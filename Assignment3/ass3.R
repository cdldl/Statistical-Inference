library(gbm)
library(ranger)
library(Metrics)



data(iris)
colnames(iris)[length(colnames(iris))] = 'y'
train_index <- sample(1:nrow(iris), nrow(iris)*0.75)
train = iris[train_index,]
test = iris[-train_index,]


error = function (w, misses) 
{
  sum(w * misses)
}

alpha = function (err, K) 
{
  if (err > 0.5) {
    1
  }
  else if (err <= 0) {
    20
  }
  else {
    log((1 - err)/err) + log(K - 1)
  }
}

weights = function (w, a, misses) 
{
  tmp_w <- w * exp(a * misses)
  tmp_w/(sum(tmp_w))
}
samplePrediction = function (sample, A) 
{
  votes.table <- votesTable(sample, A)
  prediction(votes.table)
}
weightedClassVote = function (k, sample, A) 
{
  sum(A * (sample == k))
}
votesTable = function (sample, A) 
{
  classes <- unique(sample)
  sapply(classes, function(k) weightedClassVote(k, sample, 
                                                A))
}
prediction = function (votes.table) 
{
  max.votes <- names(which(votes.table == max(votes.table)))
  ifelse(length(max.votes) == 1, max.votes, sample(max.votes,1))
}
finalPrediction  = function (C, A) 
{
  apply(C, 1, function(sample) samplePrediction(sample, A))
}



same = function (formula, data, test, m = 5) 
{
  outcome.label <- toString(formula[[2]])
  y <- data[, outcome.label]
  C <- matrix(nrow = nrow(test), ncol = m)
  n <- nrow(data)
  w <- rep(1/n, n)
  A <- numeric(m)
  K <- nlevels(y)
  for (i in 1:m) {
    t <- data[sample(n, n, replace = T, prob = w), ]
    t[, outcome.label] <- droplevels(t[, outcome.label])
    fit <- ranger(formula, data = t)
    C[, i] <- as.character(predict(fit, test)$predictions)
    h <- predict(fit, data)$predictions
    levels(h) <- levels(y)
    misses <- as.numeric(h != y)
    err <- error(w, misses)
    A[i] <- alpha(err, K)
    w <- weights(w, A[i], misses)
  }
  finalPrediction(C, A)
}

preds_same = same( y ~., train,test)
model_gbm = gbm( y ~.,data= train)
best.iter <- gbm.perf(model_gbm, method = "OOB")
preds_gbm = predict(model_gbm,newdata=test, n.trees = best.iter,type='link')
preds_gbm2 = matrix(preds_gbm,ncol=nlevels(train$y))
colnames(preds_gbm2) = colnames(preds_gbm[,,1])
preds_gbm2 = apply(preds_gbm2,1,function(x) colnames(preds_gbm2)[which.max(x)])
acc_gbm = accuracy(preds_gbm2,test$y)
acc_same = accuracy(preds_same,test$y)
