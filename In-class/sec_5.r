library(mvtnorm)
library(glmnet)

rm(list = ls())
set.seed(1234)

p <- 200
n <- 200
alpha <- (1:p)^{-2}
alpha_norm <- alpha/sqrt(sum(alpha^2))

gamma <- alpha_norm + sample(rnorm(p,sd = alpha_norm),p)
gamma_norm <- gamma/sqrt(sum(gamma^2))



X <-  rmvnorm(n,sigma = diag(rep(1,p)))
offset <- X%*%alpha_norm
mean_w <- X%*%gamma_norm
size_het <- 0
eff <- 1 + (X[,1] + X[,2])*size_het


sigma_w <- sqrt(0.1)
sigma_y <- sqrt(0.1)


B <- 100



results <- matrix(0, ncol = 4, nrow =B)

for (b in 1:B){


	noise_y <- rnorm(n,sd = sigma_w)
	noise_w <- rnorm(n,sd = sigma_y)

	W <- mean_w + noise_w
	Y <- offset + eff*W + noise_y
	
	index_1 <- sample(1:n,floor(n/2))
	index_2 <- setdiff(1:n,index_1)
	
	
	X_1 <- X[index_1,]
	W_1 <- W[index_1]
	Y_1 <- Y[index_1]
	
	X_2 <- X[index_2,]
	W_2 <- W[index_2]
	Y_2 <- Y[index_2]


	
	cv_glm_y <- cv.glmnet(X,Y, family = "gaussian")
	lambda_y <- cv_glm_y$lambda.min
	glm_opt_y <- glmnet(X, Y, family = "gaussian",lambda =lambda_y)
	sel_y <- (1:p)[as.numeric(glm_opt_y$beta) != 0]

	cv_glm_w <- cv.glmnet(X,W, family = "gaussian")
	lambda_w <- cv_glm_w$lambda.min
	glm_opt_w <- glmnet(X, W, family = "gaussian",lambda =lambda_w)
	sel_w <- (1:p)[as.numeric(glm_opt_w$beta) != 0]

	index_sel <- union(sel_w,sel_y)
	X_sel <- X[,index_sel]

	res_2 <- (lm(Y~W + X_sel)$coefficients)['W']
	

	cv_glm_y_1 <- cv.glmnet(X_1,Y_1, family = "gaussian")
	lambda_y_1 <- cv_glm_y_1$lambda.min
	glm_opt_y_1 <- glmnet(X_1, Y_1, family = "gaussian",lambda =lambda_y_1)
	er_y_2 <- Y_2 - predict(glm_opt_y_1,X_2)

	cv_glm_y_2 <- cv.glmnet(X_2,Y_2, family = "gaussian")
	lambda_y_2 <- cv_glm_y_2$lambda.min
	glm_opt_y_2 <- glmnet(X_2, Y_2, family = "gaussian",lambda =lambda_y_2)
	er_y_1 <- Y_1 - predict(glm_opt_y_2,X_1)
	
	er_y <- c(er_y_1,er_y_2)


	cv_glm_w_1 <- cv.glmnet(X_1,W_1, family = "gaussian")
	lambda_w_1 <- cv_glm_w_1$lambda.min
	glm_opt_w_1 <- glmnet(X_1, W_1, family = "gaussian",lambda =lambda_w_1)
	er_w_2 <- W_2 - predict(glm_opt_w_1,X_2)

	cv_glm_w_2 <- cv.glmnet(X_2,W_2, family = "gaussian")
	lambda_w_2 <- cv_glm_w_2$lambda.min
	glm_opt_w_2 <- glmnet(X_2, W_2, family = "gaussian",lambda =lambda_w_2)
	er_w_1 <- W_1 - predict(glm_opt_w_2,X_1)
	
	er_w <- c(er_w_1,er_w_2)


	res_1 <- lm(er_y~er_w)$coefficients[2]	
	
	

	res_3 <- (lm(Y~W+X[,sample(1:p,floor(p/2))])$coefficients)['W']
	res_4 <- (lm(Y~W+X[,1:(floor(p/2))])$coefficients)['W']

	
	results[b,] <- c(res_1-mean(eff), res_2-mean(eff), res_3-mean(eff),res_4 - mean(eff)) 
	print(b)
}


RMSE <- round(sqrt(colMeans(results^2)),3)
bias <- round(colMeans(results),3)



#####################





