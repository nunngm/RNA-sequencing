library(glmnet)
lib

x = matrix(rnorm(100 * 20), 100, 20)
g4 = sample(1:4, 100, replace = TRUE)
fit3 = glmnet(x, g4, family = "multinomial")
fit3a = glmnet(x, g4, family = "multinomial", type.multinomial = "grouped")
# poisson
N = 500
p = 20
nzc = 5
x = matrix(rnorm(N * p), N, p)
beta = rnorm(nzc)
f = x[, seq(nzc)] %*% beta
mu = exp(f)
y = rpois(N, mu)
fit = glmnet(x, y, family = "poisson")
plot(fit)
pfit = predict(fit, x, s = 0.001, type = "response")
plot(pfit, y)

allgenes = counts(allData, normalized = T)
model_full = model.matrix(~age+infection:hpi+age:infection+age:infection:hpi)
dds = slpr(model_full , allgenes)


slpr = function(statistics, 
                var.groups,       
                alpha=1, 
                family="gaussian",
                lambda.via.CV=T, 
                num.cv.iter=1,
                cv.criteria="lambda.min",
                thresh=1e-7, # glmnet default
                maxit=1e5, # glmnet default
                parallel=F,
                ncores=4,
                num.predictors=NA,
                gauss.lasso=T,
                model.assessment=T) { 
  
  if (missing(statistics)) {
    stop("Univariate test statistics must be specified!")
  }
  if (missing(var.groups)) {
    stop("var.groups must be specified!")
  }  
  if (model.assessment & !gauss.lasso) {
    stop("model.assessment was requested but gauss.lasso was false!")
  }  
  
  num.groups = nrow(var.groups)
  group.names = rownames(var.groups)
  x=t(var.groups)  
  y=statistics
  intercept=T      
  cv.glmnet.results = NA
  
  if (lambda.via.CV) {
    # If requested, execute CV in parallel
    if (parallel) {
      require(doMC)
      registerDoMC(cores=ncores)
    }    
    lambda.sum = 0
    # If requested, perform CV multiple times to reduce variance associated with
    # random splitting of the data
    for (i in 1:num.cv.iter) {
      message("Executing cv.glmnet for iteration ", i)
      cv.glmnet.results = cv.glmnet(x=x, y=y, standardize=F, 
                                    alpha=alpha, family=family, intercept=intercept,
                                    thresh=thresh, maxit=maxit, parallel=parallel)
      message("...finished executing cv.glmnet for iteration ", i)
      if (cv.criteria == "lambda.min") {
        lambda = cv.glmnet.results$lambda.min
      } else {
        lambda = cv.glmnet.results$lambda.1se
      } 
      lambda.sum = lambda.sum + lambda      
    }
    mean.lambda = lambda.sum/num.cv.iter
    message("Mean CV lambda: ", mean.lambda)
    
    # Extract glmnet result on full data
    glmnet.results = cv.glmnet.results$glmnet.fit
    
    # Get index of first lambda in sequence greater than or equal to the avg. lambda
    larger.lambdas = which(glmnet.results$lambda >= mean.lambda)
    if (length(larger.lambdas) == 0) {
      lambda.index = 1
    } else {
      lambda.index = larger.lambdas[length(larger.lambdas)]
    }
    lambda = glmnet.results$lambda[lambda.index]
    
  } else {
    message("Executing glmnet...")
    glmnet.results = glmnet(x=x, y=y, standardize=F, alpha=alpha, family=family, intercept=intercept,
                            thresh=thresh, maxit=maxit)
    message("...finished executing glmnet")
    lambda = getLambdaForNumPredictors(glmnet.results, num.predictors)
  }
  
  results = list()
  results$glmnet.results=glmnet.results
  results$cv.glmnet.results=cv.glmnet.results
  results$coef.lasso = coef(glmnet.results, s=lambda)[2:(num.groups+1)] # eliminate the intercept    
  
  # if gauss.lasso is true and there was at least one
  # non-zero predictor in the lasso fit, perform an OLS regression using just the predictors
  # with non-zero coefficients in the lasso fit
  nonzero.predictors = which(results$coef.lasso != 0)
  num.nonzero.predictors = length(nonzero.predictors)
  message("Number of non-zero predictors at optimal lambda: ", num.nonzero.predictors, ", ",
          paste(nonzero.predictors, collapse=", "))
  results$slpr.aic=NA
  results$binary.aic=NA   
  results$coef.ols = rep(0, length(results$coef.lasso))
  
  if (gauss.lasso & num.nonzero.predictors > 0) {   
    
    message("Performing two-stage Gauss-Lasso estimation...")
    
    ols.formula = "y ~ 1 "
    
    # limit x to just non-zero predictors
    x = as.matrix(x[,nonzero.predictors])
    
    for (j in 1:num.nonzero.predictors) { 
      predictor.name = group.names[nonzero.predictors[j]]
      colnames(x)[j] = predictor.name
      ols.formula = paste(ols.formula, "+", predictor.name, sep="")
    }
    #message("ols.formula: ", ols.formula)
    design.mat = cbind(y, x)
    colnames(design.mat)[1] = "y"
    
    # Fit unpenalized regression model
    results$lm.results = lm(ols.formula, as.data.frame(design.mat), model=F, x=F, y=F)
    
    # Update the coef.ols to the OLS values
    
    results$coef.ols[nonzero.predictors] = results$lm.results$coefficients[2:(num.nonzero.predictors+1)]
    
    # if model.assessment is true, compute AIC values for SLPR model and 
    # a model representing multiset methods like MGSA
    
    # Get AIC for SLPR model
    results$slpr.aic = AIC(results$lm.results)
    
    if (model.assessment) {
      
      # Compute a single predictor that is 1 if the var belongs to any sets with non-zero lasso estimates
      col.sum = apply(x,1,sum)
      binary.predictor = sapply(col.sum, function(x) {
        if (x > 0) {
          return (1)
        } else {
          return (0)
        }
      })
      
      #message("Total number of vars: ", length(col.sum))
      #message("Number of vars belonging to more than one non-zero set: ", length(which(col.sum > 0)))            
      
      # Fit a model to the single binary predictor
      binary.fit = lm("y ~ x",data.frame(y=design.mat[,1], x=binary.predictor))
      
      results$binary.aic = AIC(binary.fit)
      message("SLPR OLS AIC: ", results$slpr.aic, ", binary predictor AIC: ", results$binary.aic)
    }    
  } else {
    warning("Gauss-Lasso requested but no non-zero predictors!")
  }
  
  return (results)
}

getLambdaForNumPredictors = function(glmnet.results, num.predictors) {  
  potential.lambda = which(glmnet.results$df >= num.predictors)
  if (length(potential.lambda) == 0) {
    message("Warning: no lambda had num non-zero above ", num.predictors, ", max non-zero: ", max(glmnet.results$df))
    lambda.index = which(glmnet.results$df == max(glmnet.results$df))[1]
    message("Selected lambda: ", lambda.index, ", df: ", glmnet.results$df[lambda.index])
  } else {
    lambda.index = min(potential.lambda)
  }          
  lambda = glmnet.results$lambda[lambda.index]
  return (lambda)
}
