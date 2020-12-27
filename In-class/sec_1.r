library(splines)


set.seed(1234)
rm(list = ls())


#### Loading data

data_africa <- read.table("http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/SAheart.data",
	sep=",",head=T,row.names=1)

Y <- data_africa[,'chd']
X_1 <- data_africa[,'sbp'] 
X_2 <- data_africa[,'alcohol']
X_3 <- data_africa[,'tobacco']


## Problem 1: Create a natural cubic spline with 4 degrees of freedom (5 knots) using ns function for X_1, X_2, and X_3

df_af <- 3

X_1s <- ns(X_1, df = df_af)
X_2s <- ns(X_2, df = df_af)
X_3s <- ns(X_3, df = df_af)



par(mfrow=c(2,3))

for(j in 1:df_af){
plot(cbind(X_1,X_1s[,j]),col = j+2,xlab = 'blood pressure', ylab = paste('base function',j,sep = ' '))
}
dev.off()




## Problem 2: run a logit regression (use glm function) with Y as a response and splines of X_1, X_2, and X_3 as regressors; compute the mean squared error;
## The predict.glm function (with type = 'response') will be handy

glm_res <- glm(Y~X_1s+X_2s+X_3s, family = 'binomial')
new_data <- cbind(X_1s,X_2s,X_3s)
glm_fit_base <- predict.glm(glm_res,newdata = as.data.frame(new_data),type = 'response')
MSE_base <- mean((Y -glm_fit_base)^2)




## Problem 3: do a bootstrap simulation (function sample will be useful) with B = 1000 replications, in each simulation compute
## run the logit exercise as in Problem 2, and compute the estimation error. Plot its distribution over B bootstrap samples.

B <- 1000
results_boot <- rep(0, B)
n <- length(Y)


for (b in 1:B){

	index_b <- sample(1:n,n,replace = TRUE)
	Y_b <- Y[index_b]
	X_1sb <- X_1s[index_b,]
	X_2sb <- X_2s[index_b,]
	X_3sb <- X_3s[index_b,]
	glm_res_b <- glm(Y_b~X_1sb+X_2sb+X_3sb, family = 'binomial')
	new_data <- cbind(X_1s,X_2s,X_3s)
	glm_fit_b <- predict.glm(glm_res_b,newdata = as.data.frame(new_data),type = 'response')
	MSE_b <-  mean((glm_fit_b -glm_fit_base)^2)
	results_boot[b] <- MSE_b
	print(b)
}

plot(density(results_boot),main = 'Distribution of the Estimation Error')



























































## Problem 4: Estimate the risk of the previous procedure by 2-fold cross-validation 

n <- length(Y)
random_perm <- sample(1:n, n)
index_1 <- random_perm[1:(floor(n/2))]
index_2 <- random_perm[(floor(n/2)+1):n]


glm_res_1 <- glm(Y[index_1]~X_1s[index_1]+X_2s[index_1]+X_3s[index_1], family = 'binomial')
glm_res_2 <- glm(Y[index_2]~X_1s[index_2]+X_2s[index_2]+X_3s[index_2], family = 'binomial')

new_data_1 <- cbind(X_1s[index_1],X_2s[index_1],X_3s[index_1])
new_data_2 <- cbind(X_1s[index_2],X_2s[index_2],X_3s[index_2])
glm_fit_1 <- predict.glm(glm_res_1,newdata = as.data.frame(new_data_2),type = 'response')
glm_fit_2 <- predict.glm(glm_res_2,newdata = as.data.frame(new_data_1),type = 'response')

error_glm <- rep(0,n)
error_glm[index_2] <- Y[index_2]- glm_fit_1
error_glm[index_1] <- Y[index_1]- glm_fit_2

MSE_cv <- mean(error_glm^2)






