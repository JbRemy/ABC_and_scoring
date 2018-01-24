
# 1) Implémenter une fonction ABC qui prend en argument une liste de fonctions pour les stats
# 2) Simuler selon la pseudo-posterior de la stat de ranking
# 3) Métriques d'évaluation

library(scales)


########################################################################################
# ABC Simulation
########################################################################################

# Compute the specified statistics on a sample
# INPUT : 
  # x : the sample
  # liste.stat : a list containing the stat functions
  # dim.stat : the dimension of each statistic. Each stat is a vector with equal dim
# OUTPUT:
  # matrix of size (dim.stat, nb.stats = length(liste.stat))
compute.stat <- function(x.simul, liste.stat, dim.stat){
  K <- length(liste.stat)
  res <- matrix(NA, nrow = dim.stat, ncol = K)
  for (k in 1:K){
    my.stat <- liste.stat[[k]]
    res[,k] <- eval(parse(text = my.stat))
    #res[,k] <- my.stat(x)
  }
  return (res)
}

# Run the ABC scheme
# INPUT : 
  # pi : function to simulate according to the prior on theta
  # f : function to simulate according to the conditional density of x wrt theta
  # mu0, sigma20 : prior parameters
  # stat.obs : matrix (dim.stat, nb.stats) containing the stat on the observed data
  # liste.stat : a list containing the stat functions
# OUTPUT:
  # matrix of size (n.accepted, nb.stats)
ABC <- function(n.simul, n.accepted, n.obs, pi, f, mu.0, sigma2.0, stat.obs, liste.stat){
  K <- length(liste.stat)
  dist <- matrix(NA, nrow = n.simul, ncol = K)
  
  pb <- txtProgressBar(0, n.simul, style = 3)
  theta.simul <- as.matrix(pi(n.simul, mu.0, sigma2.0))
  dim.stat <- nrow(stat.obs)
  
  for (i in 1:n.simul){
    setTxtProgressBar(pb, i)
    # Data simulation
    x.simul <- f(n.obs, theta.simul[i,], sigma2.0)
    # Statistics computations
    stat.simul <- compute.stat(x.simul, liste.stat = liste.stat, dim.stat = dim.stat)
    
    # Computing the distance between statistics
    dist[i,] <- colSums((stat.obs - stat.simul)^2)
  }
  theta.accepted <- list()
  for (k in 1:K){
    theta.accepted[[k]] <- theta.simul[order(dist[,k]),][1:n.accepted,]
  }
  close(pb)
  return (theta.accepted)  
}

# Example
if (FALSE) {
  # Experiment parameters, gaussian case
  # Model parameters
  mu.0 <- 1
  sigma2.0 <- 1
  sigma2 <- 1
  
  # pi(theta)
  pi <- function(n, mu.0, sigma2.0){
    return (rnorm(n, mu.0, sqrt(sigma2.0)))
  }
  
  # X | theta
  f <- function(n, theta, sigma2){
    return (rnorm(n, theta, sqrt(sigma2)))
  }
  
  # True posterior
  true_posterior <- function(n, x, sigma2, mu.0, sigma2.0){
    N <- length(x)
    post.sigma2 <- sigma2 * sigma2.0 / (sigma2 + N*sigma2.0)
    post.mu <- post.sigma2 * (mu.0/sigma2.0 + sum(x)/sigma2)
    return (rnorm(n, post.mu, sqrt(post.sigma2)))
  }
  
  # Dataset
  # Data parameters
  n.obs <- 1000
  n.accepted <- 100
  accept.rate <- 0.01
  n.simul <- n.accepted / accept.rate
  theta.hidden <- mu.0
  sigma2.hidden <- sigma2
  obs <- rnorm(n.obs, theta.hidden, sqrt(sigma2.hidden))
  
  stat.mean <- function(x){return(mean(x))}
  stat.var <- function(x){return(var(x))}
  
  liste.stat <- list('stat.mean(x)', 'stat.var(x)')
  stat.obs <- compute.stat(obs, liste.stat, 1)
  
  # Simul
  theta.accepted <- ABC(n.simul, n.accepted, n.obs, pi, f, mu.0, sigma2.0, stat.obs, liste.stat)
  true.posterior <- true_posterior(1000000, obs, sigma2, mu.0, sigma2.0)
  
  #Plots
  hist(theta.accepted[,1], freq = F, col = 'lightblue', border = 'blue')
  lines(density(true.posterior), col = 'red', lwd = 2)
  qqplot(true.posterior, theta.accepted[,1])
  abline(0,1, col = 'red')
  
}


########################################################################################
# Simulate according to the pseudo-posterior of the ranking stat
########################################################################################
library(MASS)
library(mvtnorm)

# Prior on the ranking parameters
prior.ranking.simul <- function(n, mu.0, sigma2.0){return(as.matrix(mvrnorm(n, mu.0, sigma2.0)))}
prior.ranking.pdf <- function(beta, mu.0, sigma2.0){return(as.matrix(dmvnorm(beta, mu.0, sigma2.0)))}


get.rank <- function(X, beta){
  return(X %*% beta)
}

# Empirical ranking risk
ERR <- function(Y, X, beta){
  n <- length(Y)
  res <- 0
  Sx <- get.rank(X, beta)
  for (i in 1:(n-1)){
    for (j in (i+1):(n-1)){
      res <- res + as.numeric(((Y[i]-Y[j])*(Sx[i]-Sx[j])) < 0)
    }
  }
  return(res * 2 / (n*(n-1)))
}

rho.pdf <- function(beta, Y, X, mu.0, sigma2.0, delta){
  return(exp(-delta*ERR(Y,X,beta))*prior.ranking.pdf(beta, mu.0, sigma2.0))
  }

# We whant to simulate from rho.pdf
# We know that rho.pdf < prior.ranking.pdf so we can use an accept-reject algo to sample from rho.pdf

rho.sampling <- function(n, d, rho.pdf, Y, X, mu.0, sigma2.0, delta, prior.ranking.pdf){
  beta <- matrix(NA, nrow = n, ncol = d)
  count <- 0
  while(count < n){
    n.to.draw <- n-count
    x <- prior.ranking.simul(n.to.draw, mu.0, sigma2.0)
    x <- matrix(x, ncol = 2)
    u <- runif(n.to.draw)
    ix <- c()
    for (i in 1:n.to.draw){
      if(u[i] < rho.pdf(x[i,], Y, X, mu.0, sigma2.0, delta) / prior.ranking.pdf(x[i,], mu.0, sigma2.0)){
        ix <- c(ix, i)
      }  
    }
    nb.accept <- length(ix)
    if (nb.accept > 0){
      beta[(count+1):(count+nb.accept),] <- x[ix,]
      count <- count + nb.accept
    }
    
  }
  return(beta)
}

get.beta.pseudo.posterior <- function(n, d, rho.pdf, Y, X, mu.0, sigma2.0, delta, prior.ranking.pdf){
  X <- apply(X, 2, function(x){return((x - mean(x))/sum(sd(x)))})
  return(colMeans(rho.sampling(n, d, rho.pdf, Y, X, mu.0, sigma2.0, delta, prior.ranking.pdf)))
}

if (FALSE){
  # Test
  n <- 100
  d <- 2
  delta <- 10
  mu.0 <- c(3,3)
  sigma2.0 <- diag(d)
  X <- rmvnorm(n, c(1,3), diag(2))
  Y <- 2*X[,1] + 4*X[,2] + rnorm(n,0,1)
  
  r <- get.rank(X, c(2,4))
  plot(Y, r)
  
  a <- c()
  ERR(Y, X, c(2,4))
  for (i in seq(0, 8, 0.1)){
    a <- c(a, ERR(Y, X, c(2,i)))
  }
  plot(seq(0, 8, 0.1),a)
  # Plus on met delta grand, plus le temps de simul augmente car le proposal dans l'accept reject est moins bon
  
  beta.simul <- rho.sampling(200, d, rho.pdf, Y, X, mu.0, sigma2.0, delta, prior.ranking.pdf)
  colMeans(beta.simul)
  hist(beta.simul[,1])
  
  hist(beta.simul[,2])
}



########################################################################################
# Evaluation Metrics
########################################################################################

get.msse <- function(theta.sim, theta.hidden){
  msse <- theta.sim - matrix(theta.hidden, nr = nrow(theta.sim), nc = length(theta.hidden), byrow = T)
  return(mean(diag(msse %*% t(msse))))
}









