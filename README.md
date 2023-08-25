## Comparison to existing packages

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
                  dConfLevel = 0.95)
lGlmFit <- GlmFit(vY ~ vTrt + vCont + vCat, dfDat)
```

```{r, echo=FALSE, results='asis'}
knitr::kable(dfRes)
```


### Use `margins::margins` to reproduce M1-M3
```{r}
# M1: Delta (model) 
lGe_model <-
  margins::margins(
    model = lGlmFit,
    variables = "vTrt",
    vcov = vcov(lGlmFit)
  )
dfGe_model <- summary(object = lGe_model)
tidyRes(dfGe_model$AME, dfGe_model$SE)
dfRes[1,]

# M2: Delta (HC2) 
lGe_HC2 <-
  margins::margins(
    model = lGlmFit,
    variables = "vTrt",
    vcov = vcovHC(lGlmFit,type = "HC2")
  )
dfGe_HC2 <- summary(object = lGe_HC2)
tidyRes(dfGe_HC2$AME, dfGe_HC2$SE)
dfRes[2,]

# M3: Delta (HC3) 
lGe_HC3 <-
  margins::margins(
    model = lGlmFit,
    variables = "vTrt",
    vcov = vcovHC(lGlmFit,type = "HC3")
  )
dfGe_HC3 <- summary(object = lGe_HC3)
tidyRes(dfGe_HC3$AME, dfGe_HC3$SE)
dfRes[3,]
```

### Use `margins::margins` to reproduce M6-M7
```{r}
#M6: Proposed (HC2)
library(margins)
lGe_HC2 <-
  margins::margins(
    model = lGlmFit,
    variables = "vTrt",
    vcov = vcovHC(lGlmFit,type = "HC2")
  )
dfGe_HC2 <- summary(object = lGe_HC2)
se_HC2 <- sqrt(dfGe_HC2$SE^2 + var(lGe_HC2$dydx_vTrt1)/nrow(dfDat))
tidyRes(dfGe_HC2$AME, se_HC2)
dfRes[6,]

# M7: Proposed (HC3) 
lGe_HC3 <-
  margins::margins(
    model = lGlmFit,
    variables = "vTrt",
    vcov = vcovHC(lGlmFit,type = "HC3")
  )
dfGe_HC3 <- summary(object = lGe_HC3)
se_HC3 <- sqrt(dfGe_HC3$SE^2 + var(lGe_HC3$dydx_vTrt1)/nrow(dfDat))
tidyRes(dfGe_HC2$AME, se_HC3)
dfRes[7,]
```

### use `stdReg::stdGlm` to for to reproduce M4

```{r}
library(stdReg)
lstdLogit <- stdReg::stdGlm(fit= lGlmFit, data=dfDat, X="vTrt")
lSummary <- summary(lstdLogit, contrast = "difference",
                    reference = "0", CI.level = 0.95)
tidyRes(lSummary$est.table[2,1],lSummary$est.table[2,2])
dfRes[4,]
```

### use `RobinCar::robincar_glm` to reproduce M5
```{r}
# use robincar_glm directly to get semi-parametric
lRobin <- RobinCar::robincar_glm(
  dfDat,
  treat_col = "vTrt",
  response_col = "vY",
  car_scheme = "simple",
  adj_method = "homogeneous",
  vcovHC = "HC3",
  g_family = stats::binomial,
  formula = vY ~ vTrt + vCont + vCat,
  contrast_h = "diff"
)
tidyRes(lRobin$contrast$result$estimate, lRobin$contrast$result$se)
dfRes[5,]
```



