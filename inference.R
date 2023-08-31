# ---- check/install/load packages ----
# ------------------------------------------------------------------------.
#install.packages("remotes")
#remotes::install_github("tye27/RobinCar")
library(RobinCar)     # Semi-parametric approach
library(sandwich)     # robust sandwich 
library(future.apply) # parallel computing
expit <- function(x) exp(x)/(1+exp(x))
# ------------------------------------------------------------------------.

########################################################################## .
# ---- Simulate Data ----
#'@param N numeric, sample size
#'@param vRandRatio vector, randomization ratio (treat:ctrl)
#'@param vCoef vector, model coefficients
#'@param bStrata logical, indicating if randomization is stratified
#'@returns A data frame containing the simulated dataset.
########################################################################## .
SimData <- function(N,
                    vRandRatio,
                    vCoef = c(-3.4,1,3,-3,0,1.5,-0.5),
                    bStrata = TRUE){
  
  if (length(N) !=1 || !is.numeric(N) || N <= 0 || N!=round(N) )
    stop("Sample size not valid, provide a positive integer number.\n")
  if (length(vRandRatio) !=2 || !is.numeric(vRandRatio) || any(vRandRatio <= 0))
    stop("Randomization ratio not valid, provide a 2x1 numeric vector of positive numbers.\n")
  if (length(vCoef) !=7 || !is.numeric(vCoef))
    stop("vCoef not valid, provide a 7x1 numeric vector.\n")
  
  # Define covariates
  vCont <- rnorm(N)
  vCat  <- rbinom(N,1,0.5)
  vCont2 <- vCont^2
  
  # Design matrix
  dfX <- data.frame(vCont,vCat, vCont2)
  
  if(bStrata){
    # Stratified Randomization
    lStrata <- split(dfX, dfX$vCat)
    lDat <- lapply(lStrata, function(x) {
      x$vTrt <- rbinom(nrow(x),size=1,prob=vRandRatio[1]/sum(vRandRatio))
      x
    })
    dfDat <-  do.call(rbind, lDat)
  }else{
    dfX$vTrt <- rbinom(N,size=1,prob=vRandRatio[1]/sum(vRandRatio))
    dfDat <- dfX
  }
  
  dfDat$vContZ <- vCont * dfDat$vTrt
  dfDat$vCont2Z <- vCont2 * dfDat$vTrt
  
  dfDat$vProb <- expit(cbind(1, as.matrix(dfDat[,c("vTrt","vCont","vCat", "vCont2", "vContZ", "vCont2Z")]))%*%vCoef)
  dfDat$vY <- rbinom(N,1,dfDat$vProb)
  dfDat$vCat <- as.factor(dfDat$vCat)
  dfDat$vTrt <- as.factor(dfDat$vTrt)
  return(dfDat)
}
########################################################################## .

########################################################################## .
# ---- Model based or sandwich covariance ----
#'@param lGlmFit 	a fitted object of class inheriting from "glm"
#'@param strCovType type of estimated covariance in the form of `method_type`. The default is model-based. 
#'@returns A matrix containing the covariance matrix estimate.
########################################################################## .
GetCovariance <- function( lGlmFit, strCovType = NULL){
  if(is.null(strCovType)) {
    strCovType <- "model"
    message("Covariance type not provided.
                Model based covariance is used by default.\n")
  }
  if(strCovType %in% c("model","RobinCar_model","Ge_model"))
  {
    mCov <- vcov(lGlmFit)
  }else if(grepl("^RobinCar_{1}.*$", strCovType)){
    mCov <- NULL
  }else if(grepl("^Ge_{1}.*$", strCovType)){
    mCov <- sandwich::vcovHC(lGlmFit,
                             gsub( "^Ge_{1}", "", strCovType))
  }else if(grepl("^EIF_{1}.*$", strCovType)){
    mCov <- vcov(lGlmFit)
  }else{
    mCov <- sandwich::vcovHC(lGlmFit,strCovType)
  }
  return( mCov )
}
########################################################################## .

########################################################################## .
# ---- A wrapper of function of glm ----
# if not convergent or probability predicted to be all 0 or 1, reduce covariate
#'@param cFormula an object of class "formula" in the form `Y~trt+covariates`
#'@param dfDat a data frame containing the variables in the model. 
#'@returns an object of class inheriting from "glm". 
########################################################################## .
GlmFit <- function(cFormula, dfDat ){
  
  nP <- length( attr(terms(cFormula),"term.labels") )
  bWarn <- FALSE
  while( nP >= 1){
    withCallingHandlers(lGlmFit <- glm(cFormula,family=binomial(link='logit'),
                                       data=dfDat),
                        warning = function(w) {
                          message(conditionMessage(w))
                          message("\n")
                          bWarn <<- TRUE})
    if(bWarn){
      cFormula <- drop.terms(terms(cFormula),
                             nP,
                             keep.response=TRUE)
      nP <- nP - 1
      bWarn <- FALSE
    }else{
      break
    }
  }
  return(lGlmFit)
}

if(FALSE){
  set.seed(3)
  dfDat <- SimData(N          = 30                  ,
                   vRandRatio = c(2,1)              ,
                   vCoef      = c(-1.3,0,2,-2,0,0,0),
                   bStrata    = T
  )
  lFit <- GlmFit(vY ~ vTrt + vCont + vCat, dfDat)
  summary( lFit )
}
########################################################################## .

########################################################################## .
# ---- tidy results up ----
#'@param dEst numeric, estimate
#'@param dStdErr numeric, estimated standard error 
#'@param dConfLevel numeric, confidence level
#'@returns a data frame with estimate, estimated standard error, lower and upper bound. 
########################################################################## .
tidyRes <- function( dEst, dStdErr, dConfLevel = 0.95){
  dfEst <- data.frame(estimate = dEst,            # estimate
                      stderr = dStdErr,           # standard deviation
                      lower  = dEst + dStdErr*qnorm(1/2-dConfLevel/2), # lower CI
                      upper  = dEst + dStdErr*qnorm(1/2+dConfLevel/2), # upper CI
                      measure = "RD"
  )
  rownames(dfEst) <- NULL
  return(dfEst)
} 
########################################################################## .

########################################################################## .
# ---- Efficient influence function (EIF) ----
#'@param dP1Est numeric, estimated proportion for Group 1
#'@param dP0Est numeric, estimated proportion for Group 0
#'@param vEIF1 vector, estimated efficient influence function for estimated proportion 1
#'@param vEIF0 vector, estimated efficient influence function for estimated proportion 0
#'@param dConfLevel numeric, confidence level
#'@param dAdjust numeric, `HC0` = 1, `HC1` = n/(n-p), `HC3` = (n/(n-p))^2, same formulation as in `RobinCar`
#'@returns a data frame with estimate, estimated standard error, lower and upper bound. 
########################################################################## .
EIFInf <- function( dP1Est, dP0Est, vEIF1, vEIF0,
                    dConfLevel = 0.95,
                    dAdjust = 1){
  dfRes <-  NULL
  dRD <- dP1Est - dP0Est
  dStdErrRD <- sqrt(cov(vEIF1-vEIF0)/length(vEIF1)*dAdjust)
  dfRes <- tidyRes(dRD, dStdErrRD, dConfLevel)
  return( dfRes )
}
########################################################################## .

########################################################################## .
# ----  A wrapper of the function robincar_glm ----
#'@param cFormula an object of class "formula" in the form `Y~trt+covariates`
#'@param dfDat a data frame containing the variables in the model. 
#'@param dConfLevel numeric, confidence level
#'@param dAdjust numeric, `HC0` = 1, `HC1` = n/(n-p), `HC3` = (n/(n-p))^2, same formulation as in `RobinCar`.
#'@returns a data frame with estimate, estimated standard error, lower and upper bound. 
########################################################################## .
RobinCarInf <- function( cFormula, dfDat, 
                         dConfLevel=0.95, dAdjust=1){
  
  vAllVars <- all.vars( cFormula )
  
  lRes <- RobinCar::robincar_glm(
    dfDat,
    treat_col = vAllVars[2],
    response_col = vAllVars[1],
    car_scheme = "simple",
    adj_method = "homogeneous",
    vcovHC = "HC0",
    g_family = stats::binomial,
    formula = cFormula ,
    contrast_h = "diff"
  )
  
  dRD <- lRes$contrast$result$estimate
  dSd <- lRes$contrast$result$se * sqrt(dAdjust)
  dfRes <- tidyRes(dRD, dSd, dConfLevel)
  
  return(dfRes)
}
########################################################################## .

########################################################################## .
# ----  Proposed/Delta method based Variance Estimation ----
#'@param dP1Est numeric, estimated proportion for Group 1
#'@param dP0Est numeric, estimated proportion for Group 0
#'@param vD1 vector, estimated partial derivative for estimated proportion 1
#'@param vD0 vector, estimated partial derivative for estimated proportion 0
#'@param dConfLevel numeric, confidence level
#'@param dCorrect numeric, 0 for the Delta method. 
#'@returns a data frame with estimate, estimated standard error, lower and upper bound. 
########################################################################## .
GcompInf <- function( dP1Est, dP0Est, vD1, vD0, mCov,
                      dConfLevel = 0.95,
                      dCorrect = 0){
  dfRDEst <- NULL
  dRD <- dP1Est - dP0Est
  dSd <- sqrt((vD1 - vD0)%*% mCov %*%t(vD1-vD0) + dCorrect)
  dfRes <- tidyRes(dRD, dSd, dConfLevel)
  return( dfRes )
}
########################################################################## .

########################################################################## .
# ---- Description: G Computation ----
# variance computed using Delta
#'@param dfDat a data frame containing the variables in the model. 
#'@param cFormula an object of class "formula" in the form `Y~trt+covariates`
#'@param vVcov vector, methods of variance estimation. Select from `c("HC3", 
#'"const", "HC", "HC0", "HC1", "HC2", "HC4", "HC4m", "HC5","Ge_HC3", "Ge_const", 
#'"Ge_HC", "Ge_HC0", "Ge_HC1", "Ge_HC2", "Ge_HC4", "Ge_HC4m", "Ge_HC5", 
#'"RobinCar_HC0","RobinCar_HC1", "RobinCar_HC3", "EIF_HC0", "EIF_HC1", "EIF_HC3")`. 
#'@param dConfLevel numeric, confidence level
#'@returns a data frame with estimate, estimated standard error, lower and upper bound. 
########################################################################## .
GcompFit <- function( dfDat,
                      cFormula,
                      vVcov = c("HC3") ,
                      dConfLevel = 0.95
){
  
  mX <- model.matrix( cFormula, data = dfDat)
  if( qr(mX)$rank < ncol(mX) ){
    # if not full rank, reduce covariates first
    vRemove <- caret::findLinearCombos(mX)$remove - 1
    cFormula <- drop.terms(terms(cFormula), vRemove, keep.response=TRUE)
    mX <- model.matrix( cFormula, data = dfDat)
  }
  
  # Fit a model using complete cases
  lGlmFit <- GlmFit(cFormula, dfDat)
  cFormula <- lGlmFit$formula
  mX <- model.matrix( cFormula, data = dfDat)
  
  vBeta <- lGlmFit$coefficients
  strModel <- class(lGlmFit)
  
  nTrtInd <- 2
  nN <- nrow(mX)
  
  # ---- G computation ----
  dfDat0 <- dfDat1 <- dfDat
  dfDat1[, all.vars(cFormula)[2]] <- factor(1, levels = c("0","1"))
  dfDat0[, all.vars(cFormula)[2]] <- factor(0, levels = c("0","1"))

  mX1 <- model.matrix( cFormula, data = dfDat1)
  mX0 <- model.matrix( cFormula, data = dfDat0)
  
  vY1Pred <- predict(lGlmFit, newdata = dfDat1,type = "response")
  vY0Pred <- predict(lGlmFit, newdata = dfDat0,type = "response")
  vYPred  <- predict(lGlmFit, type = "response")
  
  mLamBetaRHS <- nrow(mX) * t(mX * c(lGlmFit$y - vYPred))
  
  dP1Est <- mean(vY1Pred)
  dP0Est <- mean(vY0Pred)
  
  vA1 <- vY1Pred*(1-vY1Pred)
  vA0 <- vY0Pred*(1-vY0Pred)
  nN <- nrow(mX)
  nDim <- ncol(mX)
  
  vD1 <- t(vA1) %*% mX1 / nN
  vD0 <- t(vA0) %*% mX0 / nN
  vTrtInd <- mX[,nTrtInd] == 1
  
  dCorrect <- var(vY1Pred - vY0Pred)/nN
  
  dfRes <- NULL
  
  for(j in 1:length(vVcov)){
    
    mCov <- GetCovariance(lGlmFit, vVcov[j])
    
    if(grepl("^EIF_{1}.*$", vVcov[j])){
      mLamBeta <- mCov %*% mLamBetaRHS
      
      vEIF1 <- vY1Pred - dP1Est + t( vD1 %*% mLamBeta )
      vEIF0 <- vY0Pred - dP0Est + t( vD0 %*% mLamBeta )
      p <- lGlmFit$rank
      
      dAdjust <- switch (gsub("^EIF_{1}","", vVcov[j]),
                         "HC3" = (nN/(nN-p))^2,
                         "HC1" = (nN/(nN-p)),
                         "HC0" = 1
      )
      dfResTmp <- EIFInf( dP1Est, dP0Est, vEIF1, vEIF0,
                          dConfLevel = dConfLevel,
                          dAdjust = dAdjust)
      dfResTmp$`SEMethod` <- paste0("Gcomp_",vVcov[j])
      
    }else if(vVcov[j] %in% c("model","HC3", "const", "HC", "HC0",
                             "HC1", "HC2", "HC4", "HC4m", "HC5")){
      
      dfResTmp <- GcompInf( dP1Est, dP0Est, vD1, vD0, mCov,
                            dConfLevel = dConfLevel,
                            dCorrect = dCorrect)
      dfResTmp$`SEMethod` <- paste0("Gcomp_Imp_",vVcov[j])
      
    }else if(grepl("^Ge_{1}.*$", vVcov[j])){
      dfResTmp <- GcompInf( dP1Est, dP0Est, vD1, vD0, mCov,
                            dConfLevel = dConfLevel,
                            dCorrect = 0)
      dfResTmp$`SEMethod` <- paste0("Gcomp_",vVcov[j])
    }else if(grepl("^RobinCar_{1}.*$", vVcov[j])){
      p <- lGlmFit$rank
      dAdjust <- switch(gsub("^RobinCar_{1}","", vVcov[j]),
                        "HC3" = (nN/(nN-p))^2,
                        "HC1" = (nN/(nN-p)),
                        "HC0" = 1
      )
      dfResTmp <- RobinCarInf( cFormula, dfDat, 
                               dConfLevel = dConfLevel,
                               dAdjust = dAdjust)
      dfResTmp$`SEMethod` <- paste0("Gcomp_",vVcov[j])
    }else{
      stop("vVcov not supported.")
    }
    dfRes <- rbind(dfRes,dfResTmp)
  }
  
  return(dfRes)
  
}
########################################################################## .
