#install.packages('CVXR')

library(mvtnorm)
library(CVXR)


rm(list = ls())
set.seed(1234)

my_density_function <- function(x,K,deg){

	x_min <- min(x)
	x_max <- max(x)
	range_x <- x_max - x_min
	low_x <- x_min - 0.2*range_x
	up_x <- x_max + 0.2*range_x
	range_full <-up_x- low_x
	splits <- seq(low_x,up_x,length.out = K+1)
	mesh_size <- splits[2] - splits[1]
	centers <- (splits[-1]+splits[-(K+1)])/2
	counts <- as.vector(table(cut(x,splits,include.lowest = TRUE)))
	scale <- sum(counts)*mesh_size
	
	data_matrix <- splines::ns(centers, df = deg)
	pois_reg_res <- glm(counts~data_matrix, family = 'poisson')
	freq_pois <- exp(pois_reg_res$linear.predictors)
	dens_pois <- freq_pois/scale
	

	return(cbind(centers,freq_pois,dens_pois))
}
	



p <- 10
n <- 100
beta_0 <- (1:p)^{-2}
beta_0_norm <- beta_0/sqrt(sum(beta_0^2))
gamma <- beta_0_norm 
tau_av <- 1
sigma_0 <- sqrt(0.2)
sigma_1 <- sqrt(0.2)


X <-  rmvnorm(n,sigma = diag(rep(1,p)))
logit <- X%*%gamma*2 
pi <- exp(logit)/(1+exp(logit))
pi_tau <- as.numeric(pi > 0.5)
tau_het <- 1 + (2*pi_tau-1)*1


noise_0 <- rnorm(n,sd = sigma_0)
noise_1 <- rnorm(n,sd = sigma_1)
W <- rbinom(n,1,pi)
Y_0 <- X%*%beta_0_norm + noise_0
Y_1 <- tau_het + X%*%beta_0_norm + noise_1
tau_cond <- mean(tau_het)
Y <- Y_0*(1-W) + Y_1*W


## Problem 1: Estimate tau using a linear regression of Y on W and X, and using CVXR to solve the balancing problem

tau_ols <- lm(Y~W + X)$coefficients[2]


w_unit <- Variable(n,name = 'weights')	
obj <- Minimize(sum(w_unit^2)/(n^2))
constr <- list(
	t(w_unit)%*%W/n ==1,
	t(w_unit)%*%X/n ==0,
	sum(w_unit)/n == 0
)

problem <- Problem(obj, constraints  = constr)
result <- psolve(problem)
weights_init <- result[[1]]

tau_bal <- t(weights_init)%*%Y/n


cbind(round(weights_conv,3),W)















## Problem 2: Estimate tau using CVXR imposing the non-negativeity constraints 


w_unit <- Variable(n,name = 'weights')	
obj <- Minimize(sum(w_unit^2)/(n^2))
constr <- list(
	t(w_unit)%*%W/n ==1,
	t(w_unit)%*%X/n ==0,
	sum(w_unit)/n == 0,
	W*w_unit >=0,
	(1-W)*w_unit <=0
	
)
problem <- Problem(obj, constraints  = constr)
result <- psolve(problem)
weights_conv <- result[[1]]

tau_bal_int <- t(weights_conv)%*%Y/n














## Problem 3: Repeat exericise 2 for B = 500 simulations (simulate only W and Y, fixing X) and plot the distribution of the estimators

B <- 500
results <- matrix(0, ncol = 2, nrow = B)
for (b in 1:B){


	W_b <- rbinom(n,1,pi)
	noise_0_b <- rnorm(n,sd = sigma_0)
	noise_1_b <- rnorm(n,sd = sigma_1)
	Y_0_b <- X%*%beta_0_norm + noise_0_b
	Y_1_b <- tau_het + X%*%beta_0_norm + noise_1_b

	Y_b <- Y_0_b*(1-W_b) + Y_1_b*W_b
	
	w_unit_b <- Variable(n,name = 'weights')	
	obj_b <- Minimize(sum(w_unit_b^2)/(n^2))
	constr_b <- list(
		t(w_unit_b)%*%W_b/n ==1,
		t(w_unit_b)%*%X/n ==0,
		sum(w_unit_b)/n == 0,
		W_b*w_unit_b >=0,
		(1-W_b)*w_unit_b <=0
	
	)
	problem <- Problem(obj_b, constraints  = constr_b)
	result_b <- psolve(problem)
	weights_b <- result_b[[1]]
	tau_bal_int_b <- as.vector(t(weights_b)%*%Y_b/n)
	
	tau_ols_b <- lm(Y_b~W_b+X)$coefficients[2]
	
	results[b,] <- c(tau_bal_int_b, tau_ols_b)
	print(b)
}

bal_est <- my_density_function(results[,1],100,deg = 3)[,c(1,3)]
ols_est <- my_density_function(results[,2],100,deg = 3)[,c(1,3)]


plot(bal_est,main = 'Distribution of the estimators', col = 'black',lty = 1,ylim = c(0,2.1),xlim = c(0,2),type = 'l',xlab = 'Estimates',ylab = 'Density')
lines(ols_est, col = 'red',lty = 2)
abline(v = mean(results[,1]), lty = 2,lwd = 0.5,col = 'black')
abline(v = mean(results[,2]), lty = 2,lwd = 0.5,col = 'red')
abline(v =0, lty = 2,lwd = 0.5,col = 'blue')
abline(v = 2, lty = 2,lwd = 0.5,col = 'blue')
abline(v =tau_cond, lty = 2,lwd = 0.5,col = 'blue')
legend('topleft',col <- c('black','red',), lty = c(1,2),legend = c('balancing','ols'))













