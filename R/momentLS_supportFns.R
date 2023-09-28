varTruthCalc <- function(F_weights, F_support){
  
  varTruth = sum(F_weights * (1 + F_support)
                 /(1 - F_support))
  
  return(varTruth)
  
}

avar_momentLS = function(x,delta){
  mLS = SR1(r = autocov(x),delta = delta)
  return(asympVariance(weights = mLS$weights,support = mLS$support))
}

tuningAVE = function(B, L, c, chainParams, nBest){
  
  #' tuning function for momentLS asymptotic variance estimator
  #' runs all combinations of provided L and c
  #' provides requested number of best estimates
  #' ...calculated via lowest n asymptotic variance MSE
  #' 
  #' @param B number of chains to generate via MH alg
  #' @param L vector of 'number of splits' to try
  #' @param c vector of c values to try
  #' @param chainParams list of requested chain parameters
  #' @param nBest number of best estimator combinations to return
  
  # creating combinations of L and c
  momentLSestimators <- data.table(expand.grid(L,c))
  names(momentLSestimators) <- c("L","c")
  
  # initializing vectors and matrices
  resultMat <- matrix(nrow=B,ncol = nrow(momentLSestimators)+1)
  delta_tilde <- vector(length=nrow(momentLSestimators))
  momentLS_avar <- vector(length = nrow(momentLSestimators))
  
  for(b in 1:B){
    if(b/10==round(b/10)){cat(b,"/",B,"\n",sep = "")}
    
    set.seed(497+b)
    
    # generates chain with given parameters
    ch = generateMHChain(M = 10000, nStates = 100, 
                         discreteMC = chainParams$discreteMC, g = chainParams$g, 
                         d = chainParams$d)
    # true parameter value
    varTruth = varTruthCalc(ch$F_weights, ch$F_support)
    
    # tuning delta parameter
    for (i in 1:nrow(momentLSestimators)){
      delta_tilde[i] = as.numeric(tune_delta(ch$x,nSplits = 
                                               as.numeric(
                                                 momentLSestimators[i,1]))$delta*
                                    momentLSestimators[i,2])
    }
    
    # computing momentLS AVE with delta_tilde
    for (i in 1:nrow(momentLSestimators)){
      momentLS_avar[i] = avar_momentLS(x = ch$x,delta = delta_tilde[i])
    }
    
    resultMat[b,] =
      c("Truth"=varTruth, momentLS_avar)
  }
  
  aVarMSE <- vector(length = nrow(momentLSestimators))
  
  # MSE of AVE
  for (i in 1:nrow(momentLSestimators)){
    aVarMSE[i] = mean((resultMat[,i+1] - resultMat[,1])^2)
  }
  
  # returns n best combinations and respective aVarMSEs
  best_estimators = momentLSestimators[order(aVarMSE)[1:nBest],]
  
  best_aVarMSE = aVarMSE[order(aVarMSE)[1:nBest]]
  
  results_df = data.frame(best_estimators, best_aVarMSE)
  
  return(results_df)
  
}
