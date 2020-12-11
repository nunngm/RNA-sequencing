# Loading the library
library(glmnet)
library(caret)
set.seed(31138)

# Method to evaluate final results
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}

#Find the AUGs at 0.25 hpi
data = compare.group(hpi = "0", altH = "greater")

# Generate a list of AUGs at this timepoint using Tstat
genesOI = (-data[[1]]$stat+data[[2]]$stat-data[[3]]$stat+data[[4]]$stat)/sqrt(4)
y_resp = as.numeric(genesOI)

# Collect the appropriate count data
countData = counts(allData, normalized = T)[, sub_0] #only 0.25 hpi time point

# Randomly sample the genes
sample = sample(1:nrow(countData), round(0.8*nrow(countData)))
trn = countData[sample, ]
tst = countData[-sample, ]

# Linear model for the genes of interest, Recommended folds = 10, all samples are in the same units so standardize = F, covariate regression therefore alpha must be 1
lambdas <- 10^seq(2, -3, by = -.1)
lasso_reg <- cv.glmnet(trn , y_resp[sample], alpha = 1, lambda = lambdas, standardize = F, nfolds = 10)

# Find the best lambda to use moving forward
lambda_best <- lasso_reg$lambda.min 
lambda_best

#Which samples are important for the model?
x = as.matrix(coef(lasso_reg, lasso_reg$lambda.min))
x = x[x[, 1]!=0, ]
x

# Generate the optimal model
lasso_model <- glmnet(trn, y_resp[sample], alpha = 1, lambda = lambda_best, standardize = F)

# Test the performance of the model
predictions_test = predict(lasso_model, s = lambda_best, newx = tst)
eval_results(y_resp[-sample], predictions_test, tst)
#Large RMSE and low R-squared - not very accurate model

# let's check it anyways
predictions = predict(lasso_model, s = lambda_best, newx = countData)
eval_results(y_resp, predictions, countData)


##Find genes which are AUGs at a given timepoint -> to be fed into LASSO
# Perform a conservative pearson correlation on 
data = compare.group(altH = "greater")

# Now incorporates directionality ... what genes 
genesOI = rowMeans(cbind((-data[[1]]$stat+data[[2]]$stat-data[[3]]$stat+data[[4]]$stat)/sqrt(4), (data[[8]]$stat-data[[7]]$stat+data[[6]]$stat-data[[5]]$stat)/sqrt(4), (data[[12]]$stat-data[[11]]$stat+data[[10]]$stat-data[[9]]$stat)/sqrt(4)))

## 
genesOI = rowMeans(cbind((-data[[1]]$stat+data[[2]]$stat-data[[3]]$stat+data[[4]]$stat)/sqrt(4), (data[[8]]$stat-data[[7]]$stat+data[[6]]$stat-data[[5]]$stat)/sqrt(4), (data[[12]]$stat-data[[11]]$stat+data[[10]]$stat-data[[9]]$stat)/sqrt(4)))


genesOI = as.data.frame(cbind(genesOI, pt(q = genesOI, df = 3, lower.tail = F), objectSymbol[rownames(data[[1]])]))
rownames(genesOI) = rownames(data[[1]])

genesOI = genesOI[order(genesOI[, 1], decreasing = T), ]

temp = genesOI[genesOI$genesOI > 1.6, ]




data = counts(allData, normalized = T)[,]
y_resp = as.numeric(genesOI$genesOI)
# y_resp = factor(rep(c(rep("no", 9), rep("yes", 3)), 3))

lambdas <- 10^seq(2, -3, by = -.1)



# Setting alpha = 1 implements lasso regression
sample = sample(1:nrow(data), round(0.8*nrow(data)))
trn = data[sample, ]
tst = data[-sample, ]
lasso_reg <- cv.glmnet(trn, y_resp[sample], alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)

# Best 
lambda_best <- lasso_reg$lambda.min 
x = as.matrix(coef(lasso_reg, lasso_reg$lambda.min))
x = x[x[, 1]!=0, ]
lambda_best
x
lasso_model <- glmnet(trn, y_resp[sample], alpha = 1, lambda = lambda_best, standardize = TRUE)

predictions_train <- predict(lasso_model, s = lambda_best, newx = trn)
eval_results(y_resp[sample], predictions_train, trn)
x = as.matrix(coef(lasso_model, lasso_model$lambda.min))
x = x[x[, 1]!=0, ]

predictions_test = predict(lasso_model, s = lambda_best, newx = tst)
eval_results(y_resp[-sample], predictions_test, tst)



