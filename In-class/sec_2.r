library(glmnet)
library(splines)
library(pdp)
library(randomForest)

set.seed(1234)
rm(list = ls())


#### Loading data

data_africa <- read.table("http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/SAheart.data",
	sep=",",head=T,row.names=1)

Y <- data_africa[,'chd']
X_1 <- data_africa[,'sbp'] 
X_2 <- data_africa[,'alcohol']
X_3 <- data_africa[,'tobacco']
X_4 <- data_africa[,'obesity']



df_af <- 10

X_1s <- ns(X_1, df = df_af)
X_2s <- ns(X_2, df = df_af)
X_3s <- ns(X_3, df = df_af)
X_4s <- ns(X_4, df = df_af)



data_X <- cbind(X_1s ,X_2s ,X_3s, X_4s)
data_Y <- Y


## Problem 1

n <- length(Y)
n_train <- floor(n*0.8)
n_test <- n - n_train 
index_tr <- sample(1:n,n_train)
index_test <- setdiff(1:n,index_tr)

data_X_train <- data_X[index_tr,]
data_Y_train <- data_Y[index_tr]

data_X_test <- data_X[index_test,]
data_Y_test <- data_Y[index_test]


glm_res <- glmnet(data_X_train,data_Y_train, family = 'binomial',alpha = 1)
glm_cv <- cv.glmnet(data_X_train,data_Y_train, family = 'binomial',alpha = 1)
opt_lambda <- glm_cv$lambda.min
glm_res_opt <- glmnet(data_X_train,data_Y_train, family = 'binomial',alpha = 1,lambda = opt_lambda)
glm_res_nonr <- glmnet(data_X_train,data_Y_train, family = 'binomial',alpha = 1,lambda =0)

link_pred_reg <- predict(glm_res_opt,data_X_test,type = 'link')
link_pred_nonr <- predict(glm_res_nonr,data_X_test,type = 'link')
basic_pr <- log(mean(data_Y_train)/(1-mean(data_Y_train)))

emp_risk_test_reg <-mean(log(1+exp(link_pred_reg)) - data_Y_test*link_pred_reg)
emp_risk_test_nonr <-mean(log(1+exp(link_pred_nonr)) - data_Y_test*link_pred_nonr)
emp_risk_test_simp <-mean(log(1+exp(basic_pr)) - data_Y_test*basic_pr)



## Problem 2

B <- 400
results <- do.call(rbind,lapply(1:B, function(b){
	
	index_b <- sample(1:n,n,replace = TRUE)
	data_X_b <- data_X[index_b,]
	data_Y_b <- data_Y[index_b]

	
	index_tr_b <- sample(1:n,n_train)
	index_test_b <- setdiff(1:n,index_tr)

	data_X_train_b <- data_X_b[index_tr_b,]
	data_Y_train_b <- data_Y_b[index_tr_b]

	data_X_test_b <- data_X_b[index_test_b,]
	data_Y_test_b <- data_Y_b[index_test_b]

	glm_res_b <- glmnet(data_X_train_b,data_Y_train_b, family = 'binomial',alpha = 1)
	glm_cv_b <- cv.glmnet(data_X_train_b,data_Y_train_b, family = 'binomial',alpha = 1)
	opt_lambda_b <- glm_cv_b$lambda.min
	glm_res_opt_b <- glmnet(data_X_train_b,data_Y_train_b, family = 'binomial',alpha = 1,lambda = opt_lambda)
	glm_res_nonr_b <- glmnet(data_X_train_b,data_Y_train_b, family = 'binomial',alpha = 1,lambda =0)

	link_pred_reg_b <- predict(glm_res_opt_b,data_X_test_b,type = 'link')
	link_pred_nonr_b <- predict(glm_res_nonr_b,data_X_test_b,type = 'link')
	basic_pr_b <- log(mean(data_Y_train_b)/(1-mean(data_Y_train_b)))

	emp_risk_test_reg_b <-mean(log(1+exp(link_pred_reg_b)) - data_Y_test_b*link_pred_reg_b)
	emp_risk_test_nonr_b <-mean(log(1+exp(link_pred_nonr_b)) - data_Y_test_b*link_pred_nonr_b)
	emp_risk_test_simp_b <-mean(log(1+exp(basic_pr_b)) - data_Y_test_b*basic_pr_b)
	return(c(emp_risk_test_reg_b,emp_risk_test_nonr_b,emp_risk_test_simp_b))
}))



plot(density(results[,1]),main = 'Distribution of the Empirical Risk', col = 'black',xlim = c(0.4,0.75), ylim =c(0,14))
lines(density(results[,2]), col = 'red')
lines(density(results[,3]), col = 'green')

## Problem 3


data_rf <- cbind(Y,X_1,X_2,X_3,X_4)
rf_res <- randomForest(Y~X_1+X_2+X_3+X_4, data = data_rf,importance = TRUE)
pd <- partial(rf_res, pred.var = c("X_1", "X_4"))

pdp1 <- plotPartial(pd)
pdp1
