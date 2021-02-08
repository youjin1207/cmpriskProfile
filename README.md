# Facility Profiling under Competing Risks using Multivariate Prognostic Scroes Weighting

Youjin Lee and Douglas E. Schaubel

## Data

We use the Scientific Registry or Transplant Recipients (SRTR) data for data application study. Data access can be requested [here](https://www.srtr.org/requesting-srtr-data/data-requests).


## Code

- `Code/auxfunctions.R` : provides auxiliary functions to construct risk classes as a form of the matrix with a fixed dimension, to derive a prognostic score from a stratified additive hazards model, and to estimate the variance of a CIF. 

- `Code/sim1.R` : provides the data generating process under scenario (i) presented in the paper.

- `Code/sim2.R` : provides the data generating process under scenario (ii) presented in the paper.

- `Code/read_SRTR` : reads the data and implements some exploratory analysis. 
 
- `Code/analysis_SRTR` : estimates the center effects using the exact macthing, prognostic scores, and the naive method with the SRTR data. 

- `Code/result_SRTR` : summarizes ther results from the analysis. 

- `Code/timevarying_SRTR`: explores the time-varying center effects. 

## Instructions for the use of sample data

In `Data/sample.csv`, we provide a hypothetical data with the same data structure as the SRTR data. The data is provided for illustrative purpose only, not for reproducing the results. You may use this sample data to run the codes for SRTR data analysis.

```{r}
library(MASS)
library(survival)
source("Code/auxfunctions.R")

matchingdat = read.csv("Data/sampledat.csv", header = TRUE, sep = ",")
print(names(matchingdat)) # variable names

# exact matching with three key covariates: age, diabetes, bloodtype
age.cate = ifelse(matchingdat$age + 50 < 40, 1, 
                  ifelse(matchingdat$age + 50 < 50, 2, 
                         ifelse(matchingdat$age + 50 < 60, 3, 4)))
diabetes.cate = ifelse(matchingdat$diabetes == 0, 1, 2)
bloodtype.cate = ifelse(matchingdat$blood_a == 1, 1, 
                        ifelse(matchingdat$blood_b == 1, 2, 
                               ifelse(matchingdat$blood_ab == 1, 3, 4)))
matchingdat$age.cate = age.cate
matchingdat$diabetes.cate = diabetes.cate
matchingdat$bloodtype.cate = bloodtype.cate
rownames(matchingdat) = as.character(c(1:nrow(matchingdat)))

# population level risk class probability constructed by the exact matching
table.mat = ftable(matchingdat$age.cate, matchingdat$diabetes.cate, matchingdat$bloodtype.cate)
template.mat = as.matrix(table.mat) 
prop.mat = template.mat / nrow(matchingdat)


scoredat = read.csv("Data/sampledat.csv", header = TRUE, sep = ",")
rownames(scoredat) = as.character(c(1:nrow(scoredat)))
# fit a cause-specific hazards
fit1 = coxph(Surv(timeto, cause1) ~  age + bmi + female +     
               race_black + race_hispanic + race_asian +  
               race_other + blood_a + blood_ab +     
               blood_b + diag_poly + diag_diab +     
               diag_hyper + diag_other + copd +        
               hyperten + malig + pvd +       
               diabetes + strata(center), 
              data = scoredat)
fit2 = coxph(Surv(timeto, cause2) ~  age + bmi + female +     
               race_black + race_hispanic + race_asian +  
               race_other + blood_a + blood_ab +     
               blood_b + diag_poly + diag_diab +     
               diag_hyper + diag_other + copd +        
               hyperten + malig + pvd +       
               diabetes + strata(center), 
             data = scoredat)

first.scores.whole = predict(fit1, reference = "sample")
second.scores.whole = predict(fit2, reference = "sample")
scoredat$first.scores.whole = first.scores.whole
scoredat$second.scores.whole = second.scores.whole
first.scores.raw = c(quantile(scoredat$first.scores.whole, seq(0.2, 1, 0.2))[1:4], max(scoredat$first.scores.whole) + 0.1)
second.scores.raw = c(quantile(scoredat$second.scores.whole, seq(0.2, 1, 0.2))[1:4], max(scoredat$second.scores.whole) + 0.1)  

# population-level risk class probability constructed via multivariate prognostic scores
dist.mat = matrix(0, nrow(scoredat), 2)
score1 = scoredat$first.scores.whole
score2 = scoredat$second.scores.whole
dist.mat[,1] = apply(as.matrix(score1), 1, function(x) sum(x > first.scores.raw))
dist.mat[,2] = apply(as.matrix(score2), 1, function(x) sum(x > second.scores.raw))
table.raw = table(dist.mat[,1], dist.mat[,2])
table.raw = as.matrix(table.raw) / nrow(scoredat)
```

