set.seed(1234)
library(mvtnorm)
library(glmnet)

### Part 1 - Heterogeneity 

p <- 50
n <- 400
alpha <- (1:p)^{-2}
alpha_norm <- alpha/sqrt(sum(alpha^2))

gamma <- alpha_norm + sample(rnorm(p,sd = alpha_norm),p)
gamma_norm <- gamma/sqrt(sum(gamma^2))

X <-  rmvnorm(n,sigma = diag(rep(1,p)))
offset <- X%*%alpha_norm
mean_w <- X%*%gamma_norm
size_het <- 1
eff <- 1 + (X[,1] + X[,2])*size_het


sigma_w <- sqrt(0.1)
sigma_y <- sqrt(0.1)


noise_y <- rnorm(n,sd = sigma_w)
noise_w <- rnorm(n,sd = sigma_y)

W <- mean_w + noise_w
Y <- offset + eff*W + noise_y

index_1 <- sample(1:n,floor(n/2))
index_2 <- setdiff(1:n,index_1)

X_upd <- cbind(X,X*X[,1],X*X[,2])


## Problem 1: construct residuals \tilde Y_i, \tilde W_i, and estimate the heterogeneity
## in treatment effects using the method discussed in the class. Don't forget to use cross-fitting


### Part 2 - Optimal IV

tau <- 1
n <- 400
p <- floor(n/2)
snr <- 1
sigma_u <- 1
sigma_e <- 1
rho <- 0.5
beta_fs_unnorm <- (1:p)^(-2)
beta_fs <- snr*beta_fs_unnorm/sqrt(sum(beta_fs_unnorm))


Z <-  rmvnorm(n,sigma = diag(rep(1,p)))
mean_W <- as.numeric(Z%*%beta_fs)
noise_u <- rnorm(n,sd = sigma_u)
noise_e <- rho*noise_u + sqrt((1-rho^2))*rnorm(n,sd = sigma_e)
W <- mean_W + noise_u
Y <- tau*W + noise_e


### Problem 2: construct 5 different estimators for tau:
#1) ols (Y on W)
#2) TSLS using random 5 instruments
#3) TSLS using all the instruments
#4) Feasible optimal IV (don't forget about cross-fitting)
#5) Infeasible optimal IV (that uses mean_W)