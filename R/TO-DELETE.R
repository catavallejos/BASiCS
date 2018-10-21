
library(mvtnorm)

library(mvnfast)

my.fun <- function(q) {
  x <- rep(0, q) 
  x[seq_len(q-1)] <- rmvn(n = 1, mu = rep(0, q-1), 
                          sigma = diag(q-1) - 1/q, ncores = 3)
  x[q] <- 0 * q - sum(x[seq_len(q-1)] ) 
  return(x)
}

nRep <- 500

q <- 5
set.seed(q)
A5 <- matrix(0, ncol = 1000, nrow = q)
for(i in seq_len(nRep)) {
  A5[,i] <- my.fun(q)  
}


q <- 10
set.seed(q)
A10 <- matrix(0, ncol = 1000, nrow = q)
for(i in seq_len(nRep)) {
  A10[,i] <- my.fun(q)  
}

q <- 100
set.seed(q)
A100 <- matrix(0, ncol = 1000, nrow = q)
for(i in seq_len(nRep)) {
  A100[,i] <- my.fun(q)  
}

q <- 1000
set.seed(q)
A1000 <- matrix(0, ncol = 1000, nrow = q)
for(i in seq_len(nRep)) {
  print(i)
  A1000[,i] <- my.fun(q)  
}

#q <- 10000
#set.seed(q)
#A10000 <- matrix(0, ncol = 100, nrow = q)
#for(i in seq_len(nRep)) {
#  A10000[,i] <- my.fun(q)  
#}

par(mfrow = c(2,2))
hist(A5[5,])
hist(A10[10,])
hist(A100[100,])
hist(A1000[1000,])

var(A5[5,])
var(A10[10,])
var(A100[100,])
var(A1000[1000,])

mean(A5[5,])
mean(A10[10,])
mean(A100[100,])
mean(A1000[1000,])

par(mfrow = c(2,2))
plot(apply(A5, 1, mean), apply(A5, 1, var))
points(mean(A5[5,]), var(A5[5,]), pch = 16, col = "red")

plot(apply(A10, 1, mean), apply(A10, 1, var))
points(mean(A10[10,]), var(A10[10,]), pch = 16, col = "red")

plot(apply(A100, 1, mean), apply(A100, 1, var))
points(mean(A100[100,]), var(A100[100,]), pch = 16, col = "red")

plot(apply(A1000, 1, mean), apply(A1000, 1, var))
points(mean(A1000[1000,]), var(A1000[1000,]), pch = 16, col = "red")


par(mfrow = c(1,2))
plot(c(mean(A5[5,]), mean(A10[10,]), mean(A100[100,]), mean(A1000[1000,])),
     type = "b", ylim = range(apply(A5, 1, mean)))
lines(c(mean(A5[1,]), mean(A10[1,]), mean(A100[1,]), mean(A1000[1,])),
      type = "b", col = "red")
abline(h = 0, col = "blue")

plot(c(var(A5[5,]), var(A10[10,]), var(A100[100,]), var(A1000[1000,])),
     type = "b")
lines(c(var(A5[1,]), var(A10[1,]), var(A100[1,]), var(A1000[1,])),
      type = "b", col = "red")
abline(h = 0, col = "blue")

par(mfrow = c(2,2))
hist(apply(A5, 1, mean))
abline(v = mean(A5[5,]), col = "red")
hist(apply(A10, 1, mean))
abline(v = mean(A10[10,]), col = "red")
hist(apply(A100, 1, mean))
abline(v = mean(A100[100,]), col = "red")
hist(apply(A1000, 1, mean))
abline(v = mean(A1000[1000,]), col = "red")

par(mfrow = c(2,2))
hist(apply(A5, 1, var))
abline(v = var(A5[5,]), col = "red")
hist(apply(A10, 1, var))
abline(v = var(A10[10,]), col = "red")
hist(apply(A100, 1, var))
abline(v = var(A100[100,]), col = "red")
hist(apply(A1000, 1, var))
abline(v = var(A1000[1000,]), col = "red")

plot(density(A1000[1,]), xlim = c(-0.1, 0.1))
lines(density(A1000[10,]))
lines(density(A1000[1000,]))

library(corrplot)
corrplot(var(t(A5)), is.corr = FALSE)
corrplot(cor(t(A10)))
corrplot(cor(t(A100)))

## I need to show that the generated random variables are iid,
## but that iid somehow breaks when we have too many variables
