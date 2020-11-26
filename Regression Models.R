# Loading the library
library(glmnet)
library(caret)

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
# Loading the data
data(swiss)

x_vars <- model.matrix(Fertility~. , swiss)[,-1]
y_var <- swiss$Fertility
lambda_seq <- 10^seq(2, -2, by = -.1)

# Splitting the data into test and train
set.seed(86)
train = sample(1:nrow(x_vars), nrow(x_vars)/2)
x_test = (-train)
y_test = y_var[x_test]

cv_output <- cv.glmnet(x_vars[train,], y_var[train],
                       alpha = 1, lambda = lambda_seq, 
                       nfolds = 5)
cv_output <- cv.glmnet(x_vars, y_var,
                       alpha = 1, lambda = lambda_seq, 
                       nfolds = 5)

# identifying best lamda
best_lam <- cv_output$lambda.min
best_lam


##Find genes which are AUGs at a given timepoint -> to be fed into LASSO
# Perform a conservative pearson correlation on 
temp = compare.group()
data = temp[[4]]
data$padj =  -2*(log(1-temp[[4]]$padj)+log(1-temp[[2]]$padj)+log(temp[[1]]$padj))
data$padj =  pchisq(data$padj,2*3)


genesOI = find.volcano("12")

genesOI$padj = (1-genesOI$padj)*(sign(genesOI$log2FoldChange))
genesOI = genesOI[order(genesOI$padj, decreasing = T), ]

data = counts(allData, normalized = T)[,]
y_resp = factor(rep(c(rep("no", 9), rep("yes", 3)), 3))

lambdas <- 10^seq(2, -3, by = -.1)

# Setting alpha = 1 implements lasso regression
lasso_reg <- cv.glmnet(t(data), y_resp, alpha = 1, family = "binomial", lambda = lambdas, standardize = TRUE, nfolds = 5)

# Best 
lambda_best <- lasso_reg$lambda.min 
x = as.matrix(coef(lasso_reg, lasso_reg$lambda.min))
x = x[x[, 1]!=0, ]
lambda_best

lasso_model <- glmnet(t(data), y_resp, alpha = 1, lambda = lambda_best, standardize = TRUE, family = "binomial")

predictions_train <- predict(lasso_model, s = lambda_best, newx = t(data))
eval_results(as.integer(y_resp), predictions_train, t(data))
x = as.matrix(coef(lasso_model, lasso_model$lambda.min))
x = x[x[, 1]!=0, ]

predictions_test <- predict(lasso_model, s = lambda_best, newx = x_test)
eval_results(y_test, predictions_test, test)


cvfit = cv.glmnet(t(data), y_resp, family = "binomial")
warnings(cvfit = cv.glmnet(t(data), y_resp, family = "binomial"))
