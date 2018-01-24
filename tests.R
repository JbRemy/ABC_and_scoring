rm(list = ls())

setwd('/Users/marcetheve/Documents/MVA/Graphs in ML/ABC and Scoring')
list.files()
source('functions.R')

# Test
n <- 100
d <- 2
delta <- 20
mu.0 <- c(-1,1)
sigma2.0 <- diag(d)
mu.0.rank <- mu.0
sigma2.0.rank <- sigma2.0
X <- rmvnorm(n, c(1,3), diag(2))
sigma2.Y <- 1
Y <- -1*X[,1] + 1*X[,2] + rnorm(n,0,sigma2.Y)
theta.hidden <- c(-1,1)
df <- data.frame(Y = Y, x1 = X[,1], x2 = X[,2])
n.simul = 2000
n.accepted = 2000
n.obs = 100

#jpeg('sample_scoring.jpg')
par(mfrow=c(1,1))
r <- get.rank(X, theta.hidden)
plot(Y, r, pch=21, col = 'black', bg = 'purple', ylab = 'scoring')
abline(a = 0, b = 1, col = 'red', lty = 2, lwd = 2)
grid()
#dev.off()

a <- c()
ERR(Y, X, c(2,4))
for (i in seq(0, 8, 0.1)){
  a <- c(a, ERR(Y, X, c(2,i)))
}
plot(seq(0, 8, 0.1),a)
# Plus on met delta grand, plus le temps de simul augmente car le proposal dans l'accept reject est moins bon
#beta.simul <- rho.sampling(10, d, rho.pdf, Y, X, mu.0, sigma2.0, delta, prior.ranking.pdf)
#colMeans(beta.simul)

# hist(beta.simul[,1])
# hist(beta.simul[,2])


prior.theta.simul <- function(n, mu.0, sigma2.0){return(rmvnorm(n, mu.0, sigma2.0))}
f <- function(n, theta, sigma2.0){return(rnorm(n, X %*% theta, 1))}

stat.reg <- function(Y, X){return(solve(t(X)%*%X)%*%t(X)%*%Y)}

liste.stat.obs <- list('stat.reg(Y,X)', 'get.beta.pseudo.posterior(100, d, rho.pdf, Y, X, mu.0.rank, sigma2.0.rank, delta, prior.ranking.pdf)')
stat.obs <- compute.stat(Y,liste.stat.obs, 2)
liste.stat <- list('stat.reg(Y = x.simul,X = X)', 'get.beta.pseudo.posterior(1, d, rho.pdf, x.simul, X, mu.0.rank, sigma2.0.rank, delta, prior.ranking.pdf)')


theta.sim <- ABC(n.simul = n.simul, n.accepted = n.accepted, n.obs = n.obs, pi = prior.theta.simul, f = f,
                 mu.0 = mu.0, sigma2.0 = sigma2.0, stat.obs = stat.obs, liste.stat = liste.stat)

get.msse(theta.sim[[1]], theta.hidden)
get.msse(theta.sim[[2]], theta.hidden)

sigma.true <- solve(diag(d)/sigma2.0[1,1] + t(X) %*% X/sigma2.Y)
mu.true <- sigma.true %*% (t(X)%*%Y/sigma2.Y + mu.0/sigma2.0[1,1])
true.posterior <- mvrnorm(100000, mu.true, sigma.true)



# jpeg('posteriors_1000_100.jpg')
# par(mfrow = c(1,2))
hist(theta.sim[[1]][,1], freq = F, col = alpha('blue', 0.5), xlim = c(-4,2), ylim = c(0,4.5), 
     main = '', xlab = 'theta1')
hist(theta.sim[[2]][,1], freq = F, add = T, col = alpha('red', 0.5))
lines(density(true.posterior[,1]), col = 'green', lwd = 2)
legend(x = 'topright', legend = c('Regression', 'Scoring', 'True posterior'), col = c(alpha('blue', 0.5), alpha('red', 0.5), 'green'), lty = 0, lwd = c(1,1,2), pch = 15)

hist(theta.sim[[1]][,2], freq = F, col = alpha('blue', 0.5), xlim = c(-2,4), ylim = c(0,4.5), 
     main = '', xlab = 'theta2')
hist(theta.sim[[2]][,2], freq = F, add = T, col = alpha('red', 0.5))
lines(density(true.posterior[,2]), col = 'green', lwd = 2)
legend(x = 'topright', legend = c('Regression', 'Scoring', 'True posterior'), col = c(alpha('blue', 0.5), alpha('red', 0.5), 'green'), lty = 0, lwd = c(1,1,2), pch = 15)
# dev.off()


# write.csv(theta.sim[[1]], file = 'thetasim1_10000.csv', row.names = F)
# write.csv(theta.sim[[2]], file = 'thetasim2_10000.csv', row.names = F)
# 
# eps.max <- 1
# plot((100:(n.simul/eps.max))/n.simul, cumsum(theta.sim[[1]][,1])[100:(n.simul/eps.max)]/(100:(n.simul/eps.max)), type = 'l', ylim = c(-1.35, -0.92), 
#      col = 'blue', lwd = 2, ylab = 'Mean', xlab = 'Acceptance rate', main = 'N = 10000, delta = 10')
# lines((100:(n.simul/eps.max))/n.simul, cumsum(theta.sim[[2]][,1])[100:(n.simul/eps.max)]/(100:(n.simul/eps.max)), type = 'l', col = 'red', lwd = 2)
# abline(h = mu.true[1,1], col = 'green', lwd = 2)
# grid()

delta.seq <- seq(1,19, 2)
theta.liste <- list()
temps.calcul <- c()

pb <- txtProgressBar(0, length(delta.seq), style = 3)
for (i in 1:length(delta.seq)){
  setTxtProgressBar(pb, i)
  delta <- delta.seq[i]
  liste.stat.obs <- list('stat.reg(Y,X)', 'get.beta.pseudo.posterior(100, d, rho.pdf, Y, X, mu.0.rank, sigma2.0.rank, delta, prior.ranking.pdf)')
  stat.obs <- compute.stat(Y,liste.stat.obs, 2)
  liste.stat <- list('stat.reg(Y = x.simul,X = X)', 'get.beta.pseudo.posterior(1, d, rho.pdf, x.simul, X, mu.0.rank, sigma2.0.rank, delta, prior.ranking.pdf)')
  
  t1 <- Sys.time()
  
  theta.liste[[i]] <- ABC(n.simul = n.simul, n.accepted = n.accepted, n.obs = n.obs, pi = prior.theta.simul, f = f,
                   mu.0 = mu.0, sigma2.0 = sigma2.0, stat.obs = stat.obs, liste.stat = liste.stat)
  t2 <- Sys.time()
  temps.calcul <- c(temps.calcul, t2-t1)
  print(temps.calcul)
}
close(pb)
temps.calcul
t.calcul <- temps.calcul
t.calcul[1:2] <- t.calcul[1:2]/60
t.calcul[9] <- t.calcul[9] * 60

#jpeg('delta_time.jpg')
plot(delta.seq[1:9],t.calcul, type = 'o', pch = 21, col = 'black', bg = 'blue',
     ylab = 'time', xlab = 'delta', main = 'N = 2000' )
grid()
#dev.off()

write.csv(theta.liste[[9]][[1]], 'thetasim1_2000_delta17.csv')

theta1 <- matrix(NA, nr = 9, nc = 2)
theta2 <- matrix(NA, nr = 9, nc = 2)
for (i in 1:9){
  theta1[i,1] <- mean(theta.liste[[i]][[1]][,1])
  theta1[i,2] <- mean(theta.liste[[i]][[1]][,2])
  theta2[i,1] <- mean(theta.liste[[i]][[2]][,1])
  theta2[i,2] <- mean(theta.liste[[i]][[2]][,2])
  
}



s <- 1000:100
plot(cumsum(theta.liste[[i]][[1]][,1])[s] / s)
abline(h = mu.true)
plot(cumsum(theta.liste[[i]][[2]][,1])[s] / s)
s <- 1:200

jpeg('delta_posterior_2000_200.jpg')
par(mfrow = c(3,3))
for (i in 1:9){
  hist(theta.liste[[i]][[1]][,2][s], freq = F, col = alpha('blue', 0.5), xlim = c(-2,4), ylim = c(0,4.5), 
       xlab = 'theta2', main = paste('delta =',delta.seq[i]))
  hist(theta.liste[[i]][[2]][,2][s], freq = F, add = T, col = alpha('red', 0.5))
  lines(density(true.posterior[,2]), col = 'green', lwd = 2)
  
}
dev.off()

std <- c()
for (i in 1:9){
  std <- c(std, sd(theta.liste[[i]][[1]][,2][s]))
}
jpeg('delta_std.jpg')
par(mfrow = c(1,1))
plot(delta.seq[1:9], std, type = 'o', ylim = c(0, max(std)), xlab = 'delta')
abline(h = sqrt(sigma.true[2,2]))
grid()
dev.off()


n.simul.seq = seq(1000,8000,1000)
n.accepted = 200
delta = 5
liste.stat.obs <- list('stat.reg(Y,X)', 'get.beta.pseudo.posterior(100, d, rho.pdf, Y, X, mu.0.rank, sigma2.0.rank, delta, prior.ranking.pdf)')
stat.obs <- compute.stat(Y,liste.stat.obs, 2)
theta.liste.N <- list()
for (i in 1:length(n.simul.seq)){
  n.simul <- n.simul.seq[i]
  liste.stat <- list('stat.reg(Y = x.simul,X = X)', 'get.beta.pseudo.posterior(1, d, rho.pdf, x.simul, X, mu.0.rank, sigma2.0.rank, delta, prior.ranking.pdf)')
  theta.liste.N[[i]] <- ABC(n.simul = n.simul, n.accepted = n.accepted, n.obs = n.obs, pi = prior.theta.simul, f = f,
                          mu.0 = mu.0, sigma2.0 = sigma2.0, stat.obs = stat.obs, liste.stat = liste.stat)
  
}

var.N1 <- c()
var.N2 <- c()
for (i in 1:length(n.simul.seq)){
  var.N1 <- c(var.N1, var(theta.liste.N[[i]][[1]][,1]))
  var.N2 <- c(var.N2, var(theta.liste.N[[i]][[2]][,1]))
}

jpeg('var_N_delta5_accept200.jpg')
par(mfrow = c(1,2))
plot(n.simul.seq, var.N1, type = 'o', pch = 21, col = 'black', bg = 'blue', ylab = 'variance', xlab = 'N', main = 'Regression statistic')
grid()
plot(n.simul.seq, var.N2, type = 'o', pch = 21, col = 'black', bg = 'red', ylab = 'variance', xlab = 'N', main = 'Scoring statistic')
grid()
dev.off()

jpeg('posteriors_8000_200.jpg')
par(mfrow = c(1,2))
hist(theta.liste.N[[8]][[1]][,1], freq = F, col = alpha('blue', 0.5), xlim = c(-4,2), ylim = c(0,4.5), 
     main = '', xlab = 'theta1')
hist(theta.liste.N[[8]][[2]][,1], freq = F, add = T, col = alpha('red', 0.5))
lines(density(true.posterior[,1]), col = 'green', lwd = 2)
legend(x = 'topright', legend = c('Regression', 'Scoring', 'True posterior'), col = c(alpha('blue', 0.5), alpha('red', 0.5), 'green'), lty = 0, lwd = c(1,1,2), pch = 15)

hist(theta.liste.N[[8]][[1]][,2], freq = F, col = alpha('blue', 0.5), xlim = c(-2,4), ylim = c(0,4.5), 
     main = '', xlab = 'theta2')
hist(theta.liste.N[[8]][[2]][,2], freq = F, add = T, col = alpha('red', 0.5))
lines(density(true.posterior[,2]), col = 'green', lwd = 2)
legend(x = 'topright', legend = c('Regression', 'Scoring', 'True posterior'), col = c(alpha('blue', 0.5), alpha('red', 0.5), 'green'), lty = 0, lwd = c(1,1,2), pch = 15)
dev.off()


