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
#>      estimate     stderr     lower      upper measure
#> 1 -0.07165622 0.06655213 -0.202096 0.05878357      RD

dfRes[1,]
#>      estimate     stderr     lower      upper measure       SEMethod
#> 1 -0.07165622 0.06655213 -0.202096 0.05878356      RD Gcomp_Ge_model

# M2: Delta (HC2) 
lGe_HC2 <-
  margins::margins(
    model = lGlmFit,
    variables = "vTrt",
    vcov = vcovHC(lGlmFit,type = "HC2")
  )
dfGe_HC2 <- summary(object = lGe_HC2)
tidyRes(dfGe_HC2$AME, dfGe_HC2$SE)
#>      estimate     stderr      lower      upper measure
#> 1 -0.07165622 0.06793246 -0.2048014 0.06148896      RD
dfRes[2,]
#>      estimate     stderr      lower      upper measure     SEMethod
#> 2 -0.07165622 0.06793246 -0.2048014 0.06148895      RD Gcomp_Ge_HC2

# M3: Delta (HC3) 
lGe_HC3 <-
  margins::margins(
    model = lGlmFit,
    variables = "vTrt",
    vcov = vcovHC(lGlmFit,type = "HC3")
  )
dfGe_HC3 <- summary(object = lGe_HC3)
tidyRes(dfGe_HC3$AME, dfGe_HC3$SE)
#>      estimate     stderr      lower      upper measure
#> 1 -0.07165622 0.07495167 -0.2185588 0.07524635      RD
dfRes[3,]
#>      estimate     stderr      lower      upper measure     SEMethod
#> 3 -0.07165622 0.07495167 -0.2185588 0.07524635      RD Gcomp_Ge_HC3
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
#>      estimate     stderr      lower      upper measure
#> 1 -0.07165622 0.06889201 -0.2066821 0.06336964      RD
dfRes[6,]
#>      estimate     stderr      lower      upper measure      SEMethod
#> 6 -0.07165622 0.06889201 -0.2066821 0.06336964      RD Gcomp_Imp_HC2

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
#>      estimate     stderr      lower      upper measure
#> 1 -0.07165622 0.07582244 -0.2202655 0.07695303      RD
dfRes[7,]
#>      estimate     stderr      lower      upper measure      SEMethod
#> 7 -0.07165622 0.07582244 -0.2202655 0.07695303      RD Gcomp_Imp_HC3
```

### use `stdReg::stdGlm` to for to reproduce M4

```{r}
library(stdReg)
lstdLogit <- stdReg::stdGlm(fit= lGlmFit, data=dfDat, X="vTrt")
lSummary <- summary(lstdLogit, contrast = "difference",
                    reference = "0", CI.level = 0.95)
tidyRes(lSummary$est.table[2,1],lSummary$est.table[2,2])
#>      estimate     stderr      lower     upper measure
#> 1 -0.07165622 0.06254651 -0.1942451 0.0509327      RD
dfRes[4,]
#>      estimate     stderr      lower     upper measure      SEMethod
#> 4 -0.07165622 0.06254651 -0.1942451 0.0509327      RD Gcomp_EIF_HC0
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
#>      estimate     stderr      lower      upper measure
#> 1 -0.07165622 0.07143424 -0.2116648 0.06835233      RD
dfRes[5,]
#>      estimate     stderr      lower      upper measure           SEMethod
#> 5 -0.07165622 0.07143424 -0.2116648 0.06835233      RD Gcomp_RobinCar_HC3
```
