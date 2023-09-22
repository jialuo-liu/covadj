## G-computation for RR and OR
```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r message=FALSE, warning=FALSE}
rm(list = ls(all=T))
source("inference.R")
```

### Simulate Data
```{r message=FALSE, warning=FALSE}
set.seed(21)
dfDat <- SimData(N = 60                        ,
                 vRandRatio = c(1,1)           ,
                 vCoef =  c( -2,  0, 2, -2,   0,   0,    0 )   ,
                 bStrata = F)
```

### Generate Results
```{r message=FALSE, warning=FALSE}
dfRes <- GcompFit(dfDat = dfDat,
                  cFormula = vY~vTrt+vCont+vCat,
                  vVcov = c("Ge_model",     # M1: Delta (model) 
                            "Ge_HC2",       # M2: Delta (HC2) 
                            "Ge_HC3",       # M3: Delta (HC3) 
                            "EIF_HC0",      # M4: EIF
                            "RobinCar_HC3", # M5: Semi-parametric
                            "HC2",          # M6: Proposed (HC2)
                            "HC3"           # M7: Proposed (HC3)
                            ),
                  vMeasure = c("RR","OR"),
                  dConfLevel = 0.95)
```

```{r, echo=FALSE, results='asis'}
knitr::kable(dfRes)
```

```{r}
lGlmFit <- GlmFit(vY ~ vTrt + vCont + vCat, dfDat)
```
### use `stdReg::stdGlm` to reproduce M4

```{r}
# Rate Ratio (RR)
library(stdReg)
lstdLogit <- stdReg::stdGlm(fit= lGlmFit, data=dfDat, X="vTrt")
lSummary <- summary(lstdLogit, transform = "log",
                    contrast = "difference",
                    reference = 0, CI.level = 0.95)
tidyRes(exp(lSummary$est.table[2,1]),lSummary$est.table[2,2], strMeasure = "RR")
#>    estimate    stderr     lower    upper measure
#> 1 0.6010343 0.4529889 0.2473514 1.460441      RR
subset( dfRes, measure == "RR" & SEMethod == "Gcomp_EIF_HC0" )
#>    estimate    stderr     lower    upper measure      SEMethod
#> 7 0.6010343 0.4529889 0.2473514 1.460441      RR Gcomp_EIF_HC0


# Odds Ratio (OR)
lstdLogit <- stdReg::stdGlm(fit= lGlmFit, data=dfDat, X="vTrt")
lSummary <- summary(lstdLogit, transform = "logit",
                    contrast = "difference",
                    reference = 0, CI.level = 0.95)
tidyRes(exp(lSummary$est.table[2,1]),lSummary$est.table[2,2], strMeasure = "OR")
#>    estimate    stderr     lower    upper measure
#> 1 0.5527547 0.5217658 0.1987946 1.536952      OR
subset( dfRes, measure == "OR" & SEMethod == "Gcomp_EIF_HC0" )
#>    estimate    stderr     lower    upper measure      SEMethod
#> 8 0.5527547 0.5217658 0.1987946 1.536952      OR Gcomp_EIF_HC0
```

### use `RobinCar::robincar_glm` to reproduce M5
```{r}
# Rate Ratio (RR)
# use robincar_glm directly to get semi-parametric
logRR <- function(est){
    base <- est[1]
    cont <- log( est[-1] / base )
    return(cont)
}

lRobin <- RobinCar::robincar_glm(
    dfDat,
    treat_col = "vTrt",
    response_col = "vY",
    car_scheme = "simple",
    adj_method = "homogeneous",
    vcovHC = "HC3",
    g_family = stats::binomial,
    formula = lGlmFit$formula,
    contrast_h = logRR
)
tidyRes(exp(lRobin$contrast$result$estimate), lRobin$contrast$result$se, strMeasure = "RR")
#>    estimate    stderr     lower    upper measure
#> 1 0.6010343 0.5391235 0.2089279 1.729028      RR
subset( dfRes, measure == "RR" & SEMethod == "Gcomp_RobinCar_HC3" )
#>    estimate    stderr     lower    upper measure           SEMethod
#> 9 0.6010343 0.5391235 0.2089279 1.729028      RR Gcomp_RobinCar_HC3


# Odds Ratio (OR)
# use robincar_glm directly to get semi-parametric
logOR <- function(est){
    base <- est[1]
    cont <- log( est[-1]*(1-est[1]) / est[1]/(1-est[-1]) )
    return(cont)
}

lRobin <- RobinCar::robincar_glm(
    dfDat,
    treat_col = "vTrt",
    response_col = "vY",
    car_scheme = "simple",
    adj_method = "homogeneous",
    vcovHC = "HC3",
    g_family = stats::binomial,
    formula = lGlmFit$formula,
    contrast_h = logOR
)
tidyRes(exp(lRobin$contrast$result$estimate), lRobin$contrast$result$se, strMeasure = "OR")
#>    estimate    stderr     lower    upper measure
#> 1 0.5527547 0.6176909 0.1647225 1.854863      OR
subset( dfRes, measure == "OR" & SEMethod == "Gcomp_RobinCar_HC3" )
#>     estimate    stderr     lower    upper measure           SEMethod
#> 10 0.5527547 0.6176909 0.1647225 1.854863      OR Gcomp_RobinCar_HC3
```
