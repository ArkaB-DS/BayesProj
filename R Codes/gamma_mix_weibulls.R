# Example 1: Gamma Mixture of Weibulls

library(mcmcse)
library(ellipse)
library(coda)

set.seed(100)

barker2coin <- function(y, x, k, beta, a1, b1, c.x, c.y)
{
  itr <- 0
  while(TRUE)
  {
    itr <- itr+1
    S <- rbinom(1, 1, beta)
    if(S == 0)
    {
      return(c(x, itr))
    } 
    else
    {
      C1 <- rbinom(1, 1, c.y/(c.x + c.y))
      if(C1 == 1)
      {
        lambda <- rgamma(1, a1, b1)
        p.y <- dweibull(y, k, lambda)/c.y
        
        C2 <- rbinom(1, 1, p.y)
        if(C2 == 1)return(c(y, itr))
      } 
      else 
      {
        lambda <- rgamma(1, a1, b1)
        p.x <- dweibull(x, k, lambda)/c.x
        
        C2 <- rbinom(1, 1, p.x)
        if(C2 == 1)
        {
          return(c(x, itr))
        }
      }
    }
  }
}


mcmc <- function(N = 1e3, beta = .95, start, h = 4, k = 5, a1 = 2, b1 = 5)
{
  samp <- numeric(length = N)
  loops <- numeric(length = N)
  samp[1] <- start
  accept <- 0
  for(i in 2:N)
  {
    prop <- rnorm(1, samp[i-1], sd = sqrt(h))

    if(prop < 0) 
    {
      samp[i] <- samp[i-1]
      next
    }
    c.y <- k/(exp(1)*prop)
    c.x <- k/(exp(1)*samp[i-1])
    
    re.iter <- barker2coin(x = samp[i-1], y = prop, k = k, beta = beta, a1 = a1, b1 = b1, c.x, c.y)
    samp[i] <-  re.iter[1]
    if(samp[i] != samp[i-1]) accept <- accept+1
    loops[i] <- re.iter[2]
  }
  return(list("mcmc" = samp, "loops" = loops, "accept" = accept/N))
}


a1 <- 10
b1 <- 100
k <- 10

true.mean <- (a1/b1 * gamma(1 + 1/k)) 
foo <- gamma(1 + 2/k) - (gamma(1 + 1/k))^2
true.var <-  ((gamma(1 + 1/k))^2* a1/b1^2 + foo*(a1/b1^2 + a1^2/b1^2)) 

bark <- mcmc(N = 1e5, start = true.mean, beta = 1,  k = k, a1 = a1, b1 = b1, h = true.var)
stable_99 <- mcmc(N = 1e5, start = true.mean, beta = .99,  k = k, a1 = a1, b1 = b1, h = true.var)
stable_90 <- mcmc(N = 1e5, start = true.mean, beta = .90,  k = k, a1 = a1, b1 = b1, h = true.var)
stable_75 <- mcmc(N = 1e5, start = true.mean, beta = .75,  k = k, a1 = a1, b1 = b1, h = true.var)

save(bark, stable_99, stable_90, stable_75, file = "./Data/weibull_one_run")


############### ESS and multiple reps ############
a1 <- 10
b1 <- 100
k <- 10

B <- 1000
time <- matrix(0, nrow = B, ncol = 4)
ess <- matrix(0, nrow = B, ncol = 4)
max.l <- matrix(0, nrow = B, ncol = 4)
mean.l <- matrix(0, nrow = B, ncol = 4)

for(b in 1:B)
{
  N <- 1e5
  
  time[b,1] <- system.time(bark <- mcmc(N = N, start = true.mean, beta = 1,  k = k, a1 = a1, b1 = b1, h = true.var))[3]
  time[b,2] <- system.time(stable_99 <- mcmc(N = N, start = true.mean, beta = .99,  k = k, a1 = a1, b1 = b1, h = true.var))[3]
  time[b,3] <- system.time(stable_90 <- mcmc(N = N, start = true.mean, beta = .90,  k = k, a1 = a1, b1 = b1, h = true.var))[3]
  time[b,4] <- system.time(stable_75 <- mcmc(N = N, start = true.mean, beta = .75,  k = k, a1 = a1, b1 = b1, h = true.var))[3]
  
  ess[b,1] <- effectiveSize(bark$mcmc)
  ess[b,2] <- effectiveSize(stable_99$mcmc)
  ess[b,3] <- effectiveSize(stable_90$mcmc)
  ess[b,4] <- effectiveSize(stable_75$mcmc)
  
  
  max.l[b,1] <- max(bark$loops)
  max.l[b,2] <- max(stable_99$loops)
  max.l[b,3] <- max(stable_90$loops)
  max.l[b,4] <- max(stable_75$loops)
  
  mean.l[b,1] <- mean(bark$loops)
  mean.l[b,2] <- mean(stable_99$loops)
  mean.l[b,3] <- mean(stable_90$loops)
  mean.l[b,4] <- mean(stable_75$loops)
  
  print(b)
  
  save(time, ess, max.l, mean.l, file = "./Data/wei_sim")
}

## Reproducing Tables and Plots

pdf("wei_trace.pdf", height = 5, width = 7)
par(mfrow = c(2,2))
par(mar=c(3, 2, 1.5, 2) + 0.1)
plot.ts(tail(bark$mcmc, 1e4), main = expression(paste(beta, " = 1")) , ylab = "X", xlab = "", ylim = c(0, .26))
text(5000, .25, bark$accept, col = "red")
plot.ts(tail(stable_99$mcmc, 1e4), main = expression(paste(beta, " = .99")) ,ylab = "", xlab = "",  ylim = c(0, .26))
text(5000, .25, stable_99$accept, col = "red")
par(mar=c(3, 2, 1.5, 2) + 0.1)
plot.ts(tail(stable_90$mcmc, 1e4), main = expression(paste(beta, " = .90")) , ylab = "",  ylim = c(0, .26))
text(5000, .25, stable_90$accept, col = "red")
plot.ts(tail(stable_75$mcmc, 1e4), main = expression(paste(beta, " = .75")) ,ylab = "",  ylim = c(0, .26))
text(5000, .25, stable_75$accept, col = "red")
dev.off()

acf1 <- acf(bark$mcmc, lag.max = 50, plot = FALSE)$acf
acf99 <- acf(stable_99$mcmc, lag.max = 50, plot = FALSE)$acf
acf90 <- acf(stable_90$mcmc, lag.max = 50, plot = FALSE)$acf
acf75 <- acf(stable_75$mcmc, lag.max = 50, plot = FALSE)$acf

pdf("wei_acf.pdf", height = 4, width = 6)
plot(0:50, acf1, type= 'l', xlab = "Lags", ylab = "Autocorrelation")
lines(0:50, acf99, lty = 2)
lines(0:50, acf90, lty = 3)
lines(0:50, acf75, lty = 4)
legend("topright", title = expression(beta), legend = c("1", ".99", ".90", ".75"), lty = 1:4, col = "black")
dev.off()

round(colMeans(ess)) ## ESS
# [1] 7484 6939 4320 2501
round(apply(ess, 2, sd)/sqrt(B), 2) ## ESS sd
# [1]  7.74 12.18 13.86  9.19
round(colMeans(ess/time), 2)
# [1] 251.40 616.60 717.79 658.05
round(apply(ess/time, 2, sd)/sqrt(B), 2)
# [1] 2.68 3.50 4.59 4.27
round(colMeans(mean.l), 2) ## Mean loops
# [1] 32.00  7.63  3.97  2.55
round(apply(mean.l, 2, sd)/sqrt(B), 2) ## Mean loops sd
# 1] 3.7 0.0 0.0 0.0

round(colMeans(max.l)) ## Max loops
# [1] 1315683     604      78      32
round(apply(max.l, 2, sd)/sqrt(B), 2) # Max loops sd
# [1] 366763.76      3.44      0.37      0.12
