#### Problem 1


rm(list = ls())
set.seed(1234)


n <- 20


sigma_0 <- 1
sigma_1 <- 2
rho <- 0.5

n_1 <- floor(n/2)
n_0 <- n-n_1
tau <- 5
B <- 10000


Y_0 <- rnorm(n,sd = sigma_0)
Y_1 <- tau + Y_0*rho + sqrt(1-rho^2)*rnorm(n,sd = sigma_1)

	
W <- rep(0,n)
W[sample(1:n,n_1)] <- 1


# Task: do three simulations (each B times). In the first one fix outcomes and treat W as random. In the second one -- do the opposite. In the final one -- treat both as random. 
# In each simulation compute the difference in means, and compute its variance (over simulations)



### Problem 2

rm(list = ls())
set.seed(1234)


n <- 50
k <- 10



sigma_mu <- 1
sigma_eps <- 2
tau <- 5
mu_k <- rnorm(k,sd = sigma_mu)
B <- 1000


eps <- rnorm(n*k,sd = sigma_eps)

Y_0 <- eps + rep(mu_k,each = n)
Y_1 <- Y_0 + tau

# Task: implement two different randomization schemes: half of cities or half of units in each city. 
# Simulate B times and compute the difference in means. Compute the variance over B simulations.



