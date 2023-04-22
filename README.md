# survMISA : Mediation Identification by Splitting and Aggregation with survival data

*survMISA* is an R package for identifying the most promising mediators within a broad array of candidates for survival data.
It constructs statistics with symmetry and deploys the symmetry to obtain a data-driven threshold. And it makes inference without p-values and is able to control the false discovery rate.

# Install the package

```R
library(devtools)
install_github("wangchengLi6/survMISA")
```

Dependencies: The `devtools` package is required for `survMISA` installation. Additionally, the following packages are required: `foreach`, `glmnet`, and `survival` for package's function. If parallel computation is needed, we recommend installing the `parallel` and `doParallel` packages.

Platform Compatibility: We have successfully tested our software on Windows 11 (R version >= 4.2.2) and MacOS (R version >= 4.0.2). If you plan to use it on Linux, please ensure that you have installed the latest version of R (for assistance, please refer to https://cran.r-project.org/).

If you encounter any other issues during installation, please contact me at 2290726821@qq.com .

# Citation

Still working.

# Example

Following code will generate a toy data for illustration:

```R
n = 500;p = 1000;
A = matrix(c(rep(1,20),rep(-1,20),rep(0,p-40)),nrow = 1,ncol = p)/2
eta = matrix(c(0.5,-0.5),nrow = 2,ncol = p);gam = 0.5;xi = matrix(c(0.5,0.5),2,1)
B = matrix(c(rep(1,20),rep(0,20),rep(-1,20),rep(0,940)),nrow = p,ncol = 1)/2

set.seed(1)
# Exposures and covariates
X = matrix(rnorm(n),nrow = n)
Z = matrix(rnorm(n*2),n,2)
# Candidate mediators
eM = matrix(rnorm(n*p),n,p)
M = X %*% A + Z %*% eta + eM
# Survival outcomes
expp = exp(X %*% gam + Z %*% xi + M %*% B)
DT = sqrt(log(1-matrix(runif(n),ncol = 1))/expp*(-1) * 2)# Death time
CT = runif(n,max = 10) # Censoring time
Y = pmin(DT,CT)
Delta = (Y==DT)*1

# check the model
require(survival)
coxph(Surv(Y,Delta)~X+Z+M[,c(1:20,41:60)])
```

You can implement *survMISA* and *survMISA-DR* procedure all in one function:

```R
library(survMISA)
atest = AllinOne(
  in.X = X,
  in.Z = Z,
  in.M = M,
  in.Y = Surv(Y,Delta),
  para.kappa = 1,
  para.nu = 2,
  para.seed = 1:5,
  para.rs = "half",
  para.q = 0.1,
  para.parallel = T,
  para.freq = 0.5)
atest$deran$id.S3.deran ## The real mediators is 1,2,...,20
```

Remark: The length of `para.seed` determine the times of multiple splitting. If length is 1, then derandomization method is not applied (i.e. *survMISA* ). Otherwise the derandomization method is applied (i.e. *survMISA-DR* ).

Or you can implement *survMISA* procedure step by step:

```R
## Step1: prescreening
id.S1 = PreScr_step(in.X = X,in.Z = Z,in.M = M,para.kappa = 1,para.parallel = T);id.S1
## Step2: sample splitting
indv = SampleSplit_step(in.n = nrow(M),para.nu = 2,para.seed = 1)
indv1 = indv$indv1;indv1
indv2 = indv$indv2;indv2
## Step3: construct T1 and T2 on D1 and D2, separately
t1list = T1onD1(in.X = X,in.Z = Z,in.M = M,in.Y = Surv(Y,Delta),in.id.S1 = id.S1,in.indv1 = indv1,para.rs = "half",para.parallel = T)
t2list = T2onD2(in.X = X,in.Z = Z,in.M = M,in.Y = Surv(Y,Delta),in.id.S2 = t1list$id.S2,indv2)
## Step4: aggregation
Aggr_step(t1list$T1_part,t2list$T2_part,t1list$id.S2,0.1) # The real mediators is 1,2,...,20
```
