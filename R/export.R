#' @importFrom survival coxph Surv
#' @importFrom glmnet glmnet cv.glmnet
#' @import foreach
#' @importFrom iterators iter

#' @title Step1: Pre-screening by SIS procedure
#' @description The first step of survMISA involves pre-screening candidate mediators to a moderate level using the SIS procedure. While we recommend utilizing the X-M side information, you may also specify your preferred indicator. In cases where the number of candidate variables is extremely high, the pre-screening step can significantly impact the final result, so careful parameter selection is essential.
#' @param in.X vector, with `n` elements or matirx with `1` column, the input of one-dimensional exposures.
#' @param in.Z matrix, with `n` rows and `s` colnumns, the input of s-dimensional covariates.
#' @param in.M matrix, with `n` rows and `p` colnumns, the input of p-dimensional candidate mediators.
#' @param para.d integer, the number of mediators will be retained by SIS program. If `para.d` is not specified, we will use `para.kappa` and the sample size `n` to calculate it. Refer to `para.kappa`.
#' @param para.kappa numeric, suggested to be `1` or `2`. When `para.d` is not specified, it will be used to calculate `para.d` as \eqn{d = n/log(n) * \kappa}. We recommended to determine `para.d` through `para.kappa` rather than manually. Refer to `para.d`.
#' @param para.parallel logical, default as `FALSE`. Determines whether to use parallel computation for acceleration.
#' @param para.coresnum integer, the number of cores to use when `para.parallel` is set to `TRUE`. Default as all the cores detected by `detectCores()`.
#' @param in.stat_prescr vector with length p, the user-specified filter indicator. The larger absolute value, the more likely to be retained. When it is provided, the marginal effect of X-M will no longer be calculated.
#' @return vector with length `para.d`, a set of indices containing the mediators retained by the SIS program.
#' @export
PreScr_step = function(in.X,in.Z, # 输入X和Z可以是矩阵
                       in.M, # 中介变量输入必须是矩阵，每一列代表一个中介变量，同一行为一个样本。
                       para.kappa = 1, # 如果未指定 para.d ，则使用 para.kappa 计算 para.d
                       para.parallel = F, para.coresnum, # 是否使用并行运算
                       para.d, # 预计在预筛选步骤中保留的个数，建议小于中介变量的个数，否则预筛选不会发生。
                       in.stat_prescr ## 用户指定的筛选指标序列，绝对值越大将越有可能被保留.
){
  # in.X = dt1$X
  # in.Z = dt1$Z
  # in.M = dt1$M
  # para.kappa = 1
  # para.parallel = F
  # para.d = 100
  # in.stat_prescr
  if(para.parallel){check.parallel()}
  if(!is.matrix(in.M)){
    stop("Mediators should has a form of matrix, the row is an individual and the column is a mediator.")
  }

  n = nrow(in.M)
  p = ncol(in.M)
  mat.X = matrix(in.X,nrow = n)
  if(!missing(in.Z)){mat.Z = matrix(in.Z,nrow = n)}

  mat.M = in.M


  ## 计算需要保留的个数，并判断是否会全部保留
  if(missing(para.d)){para.d = round(n/log(n)*para.kappa)}
  if(para.d >= p){
    warning("Prescreening doesn't work as the dimensionality of candidates mediators
            is smaller than the number allowed to be retained.")
    return(1:p)
  }


  ## 直接调用用户指定的筛选指标，或者计算筛选指标
  if(!missing(in.stat_prescr)){
    ind.prescr = abs(in.stat_prescr)
  }else{
    if(missing(in.Z)){
      mat.xz = mat.X
    }else{
      mat.xz = cbind(mat.X,mat.Z)
    }

    est_a = est.linear(mat.xz,mat.M,
                       para.parallel = para.parallel,para.coresnum)
    ind.prescr = abs(est_a[1,]/est_a[2,])
  }

  ## 选择前d大的标号
  id.S1 = sort(order(ind.prescr,decreasing = T)[1:para.d]) # 预筛选步骤得到的中介变量子集的标号 S1.
  return(id.S1)
}

#' @title Step2: sample splitting
#' @description A function designed to randomly divide the sample of size n into two complementary subsets, with sizes \eqn{n_1} and \eqn{n_2}, respectively.
#' @param in.n integer, the size of whole sample.
#' @param para.nu numeric, to calculate the size of subsample as:
#' \deqn{n_1 : n_2 = \nu.}
#' @param para.seed numeric, random seed for random splitting. Use to fix the result.
#' @return a list containing two sets, `indv1` and `indv2`, which represent the indices of the two subsamples, respectively.
#' @examples
#' \code{SampleSplit_step(100,1,1)}
#' \code{SampleSplit_step(100,2,1)}
#' SampleSplit_step(100,2,2)
#' @export
SampleSplit_step = function(
    in.n # 样本量大小
    ,para.nu  = 1 # 分割比例，nu:1 = n1:n2
    ,para.seed = 1 # 随机数种子
){
  n1 = floor(para.nu/(para.nu+1) * in.n)
  n2 = in.n - n1
  print(paste0("n = ",in.n,";n1 = ",n1,";n2 = ",n2,"."))
  set.seed(para.seed)
  indv1 = sort(sample(1:in.n,n1,replace = F))
  indv2 = (1:in.n)[-indv1]

  return(list(
    indv1 = indv1,
    indv2 = indv2
  ))
}

#' @title Step3.1: Calculate \eqn{T_1^{(k)}} on \eqn{\mathcal{D}_1}
#' @description In step3.1, we will calculate \eqn{T_1^{(k)}} on \eqn{\mathcal{D}_1}. We first fit partially penalized Cox model on \eqn{\mathcal{D}_1} and \eqn{\mathcal{S}_1} for further screening and denote the retained mediators as \eqn{\mathcal{S}_2}. Then calculate \eqn{\tilde{\beta}_k} and \eqn{\hat{\alpha}_k^{(1)}} on \eqn{\mathcal{S}_2}. Finally calculate
#' \deqn{T_1^{(k)} = \begin{cases} \frac{ \hat{\alpha}_k^{(1)} \tilde{\beta}_k}{\hat{se}_{\alpha,k}^{(1)}}, & k \in \mathcal{S}_2 \\ 0, & \mbox{else}\end{cases}.}
#' @param in.X vector, with `n` elements or matirx with `1` column, the input of one-dimensional exposures.
#' @param in.Z matrix, with `n` rows and `s` colnumns, the input of s-dimensional covariates.
#' @param in.M matrix, with `n` rows and `p` colnumns, the input of p-dimensional candidate mediators.
#' @param in.Y `Surv` object, with `n` elements, the input of survival information.
#' @param in.id.S1 vector, indices of the mediators retained by the SIS program.
#' @param in.indv1 vector, indices of the individuals in the first subsample \eqn{\mathcal{D}_1}
#' @param para.rs one of `half`, `BIC` and `1se`. The strategy for selecting the penalty parameter \eqn{\lambda} in the penalized Cox model. The `half` option retains half of the mediators in \eqn{\mathcal{S}_1}. The `BIC` and `1se` options represent `lambda.min` and `lambda.1se` from the `cv.glmnet` function, respectively.
#' @param para.parallel logical, default as `FALSE`. Determines whether to use parallel computation for acceleration.
#' @param para.coresnum integer, the number of cores to use when `para.parallel` is set to `TRUE`. Default as all the cores detected by `detectCores()`.
#' @return a list with elements `T1` , `T1_part` , `id.S1` , `id.S2` , `indv` :
#'
#' `T1` is \eqn{T_1^{(k)},\forall k \in \mathcal{S}_1}.
#'
#' `T1_part` is \eqn{T_1^{(k)},\forall k \in \mathcal{S}_2}.
#'
#' `id.S1` and `id.S2` contains the indices of sets \eqn{\mathcal{S}_1} and \eqn{\mathcal{S}_2}, respectively.
#'
#' `indv` contains the indices in subsample \eqn{\mathcal{D}_1}.
#' @export
T1onD1 = function(
    in.X,in.Z # 输入X和Z可以是矩阵，也可以是变量
    ,in.M # 中介变量输入必须是矩阵，且每一列代表一个中介变量，同一行为一个样本。
    ,in.Y # 必须是 survival::Surv 对象
    ,in.id.S1 # 在prescreening步骤中被保留的中介变量的序号
    ,in.indv1 # 样本分割得到的第一部分样本个体的标号
    ,para.rs = c("half","BIC","1se") # 使用惩罚估计筛选变量的策略。
    # half策略将选择恰好保留一半变量的lambda，BIC和1se是指使用cv方法后，依据BIC选择lambda或BIC-1se。具体可以参考cv.glmnet中对两种结果的解释。
    ,para.parallel = F, para.coresnum # 是否使用并行计算以提高计算速度，可以使用并行计算的步骤是惩罚估计的cv方法。
){
  if(para.parallel){check.parallel()}
  if(!is.matrix(in.M)){
    stop("Mediators should has a form of matrix, the row is an individual and the column is a mediator.")
  }

  searchLambda = function(m,xz,y,rs,parallel,coresnum){
    mxz = cbind(m,xz)
    if(parallel){
      if(missing(coresnum)){
        cl = makeCluster(detectCores())
      }else{
        cl = makeCluster(coresnum)
      }
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
      glm1 = cv.glmnet(mxz,y,family = "cox",parallel = parallel,
                       penalty.factor = c(rep(1,ncol(m)),rep(0,ncol(xz))))
    }else{
      glm1 = cv.glmnet(mxz,y,family = "cox",parallel = parallel,
                       penalty.factor = c(rep(1,ncol(m)),rep(0,ncol(xz))))
    }
    lam = switch(
      rs,
      "BIC" = glm1$lambda.min,
      "1se" = glm1$lambda.1se,
      "half" =
        glm1$lambda[which.min(abs(glm1$nzero-ncol(xz) - ncol(m)/2))]
    )

    return(lam)
  }
  # in.X = dt1$X
  # in.Z = dt1$Z
  # in.M = dt1$M
  # in.Y = Surv(dt1$Y,dt1$Delta)
  # in.id.S1 = id1
  # in.indv1 = indv1
  # para.rs = "1se"
  # para.parallel = T

  n = nrow(in.M)
  n1 = length(in.indv1)
  p = ncol(in.M)
  size.s1 = length(in.id.S1)

  mat.X = matrix(in.X,nrow = n)[in.indv1,,drop = F]
  if(!missing(in.Z)){mat.Z = matrix(in.Z,nrow = n)[in.indv1,,drop = F]}
  if(!is.matrix(in.M)){
    stop("Mediators should has a form of matrix, the row is an individual and the column is a mediator.")
  }
  mat.sM = in.M[in.indv1,in.id.S1,drop = F]
  surv.Y = in.Y[in.indv1,,drop = F]
  if(missing(in.Z)){
    mat.xz = mat.X
  }else{
    mat.xz = cbind(mat.X,mat.Z)
  }

  # set.seed(1)
  lam = searchLambda(mat.sM,mat.xz,surv.Y,para.rs,para.parallel,para.coresnum)
  mat.mxz = cbind(mat.sM,mat.xz)
  glm2 = glmnet(mat.mxz,surv.Y,family = "cox",lambda = lam,
                penalty.factor = c(rep(1,size.s1),rep(0,ncol(mat.xz))))
  tildebeta = as.vector(coef(glm2)[1:size.s1])
  id.nzero = which(tildebeta != 0)
  id.S2 = in.id.S1[id.nzero]

  est_a = est.linear(mat.xz,mat.sM,para.parallel,para.coresnum)
  est.alpha1 = est_a[1,]/est_a[2,]

  T1 = tildebeta * est.alpha1
  T1.part = T1[id.nzero]
  return(list(
    T1 = T1,
    T1_part = T1.part,
    id.S1 = in.id.S1,
    id.S2 = id.S2,
    indv = in.indv1
  ))
}

#' @title Step3.2: Calculate \eqn{T_2^{(k)}} on \eqn{\mathcal{D}_2}
#' @description In step3.2, we will calculate \eqn{T_2^{(k)}} on \eqn{\mathcal{D}_2}.
#' On \eqn{\mathcal{D}_2}, we obtain \eqn{\hat{\alpha}_k^{(2)}}, \eqn{\hat{\beta}_k^{(2)}}, \eqn{\hat{se}_{\alpha,k}^{(2)}} and \eqn{\hat{se}_{\beta,k}^{(2)}}, with selected mediators in \eqn{\mathcal{S}_2}. Here \eqn{\hat{\alpha}_k^{(2)}} and \eqn{\hat{\beta}_k^{(2)}} are the OLS and the maximum partial likelihood estimators of \eqn{\alpha_k} and \eqn{\beta_k}, respectively. Then we construct \eqn{T_2} on \eqn{\mathcal{D}_2} as follow:
#' \deqn{ T_2^{(k)}= \begin{cases} \frac{ \hat{\alpha}_k^{(2)} \hat{\beta}_k^{(2)} }{ \hat{se}_{\alpha,k}^{(2)} \hat{se}_{\beta,k}^{(2)} }, & k \in \mathcal{S}_2 \\ 0, & \mbox{else} \end{cases}.}
#' @param in.X vector, with `n` elements or matirx with `1` column, the input of one-dimensional exposures.
#' @param in.Z matrix, with `n` rows and `s` colnumns, the input of s-dimensional covariates.
#' @param in.M matrix, with `n` rows and `p` colnumns, the input of p-dimensional candidate mediators.
#' @param in.Y `Surv` object, with `n` elements, the input of survival information.
#' @param in.id.S2 vector, indices of the mediators retained by the penalized Cox model.
#' @param in.indv2 vector, indices of the individuals in the first subsample \eqn{\mathcal{D}_2}
#' @param para.parallel logical, default as `FALSE`. Determines whether to use parallel computation for acceleration.
#' @param para.coresnum integer, the number of cores to use when `para.parallel` is set to `TRUE`. Default as all the cores detected by `detectCores()`.
#' @return a list with elements `T2_part` , `id.S2` , `indv2` :
#'
#' `T2_part` is \eqn{T_2^{(k)},\forall k \in \mathcal{S}_2}.
#'
#' `id.S2` contains the indices of the set \eqn{\mathcal{S}_2}.
#'
#' `indv2` contains the indices in subsample \eqn{\mathcal{D}_2}.
#' @export
T2onD2 = function(
    in.X,in.Z # 输入X和Z可以是矩阵，也可以是变量
    ,in.M # 中介变量输入必须是矩阵，且每一列代表一个中介变量，同一行为一个样本。
    ,in.Y # 必须是 survival::Surv 对象
    ,in.id.S2 # 在惩罚估计后，仍被保留的中介变量的序号
    ,in.indv2 # 样本分割得到的第二部分样本个体的标号
    ,para.parallel = F, para.coresnum
){
  #
  # in.X = dt1$X
  # in.Z = dt1$Z
  # in.M = dt1$M
  # in.Y = Surv(dt1$Y,dt1$Delta)
  # in.id.S2 = t1l$id.S2
  # in.indv2 = indv2
  # para.parallel = F
  if(para.parallel){check.parallel()}
  if(!is.matrix(in.M)){
    stop("Mediators should has a form of matrix, the row is an individual and the column is a mediator.")
  }
  n = nrow(in.M)
  # n2 = length(in.indv2)
  # p = ncol(in.M)
  size.s2 = length(in.id.S2)

  mat.X = matrix(in.X,nrow = n)[in.indv2,,drop = F]
  if(!missing(in.Z)){mat.Z = matrix(in.Z,nrow = n)[in.indv2,,drop = F]}
  if(!is.matrix(in.M)){
    stop("Mediators should has a form of matrix, the row is an individual and the column is a mediator.")
  }
  mat.sM = in.M[in.indv2,in.id.S2,drop = F]
  surv.Y = in.Y[in.indv2,,drop = F]
  if(missing(in.Z)){
    mat.xz = mat.X
  }else{
    mat.xz = cbind(mat.X,mat.Z)
  }

  est_a = est.linear(mat.xz,mat.sM,para.parallel,para.coresnum)
  est.alpha = est_a[1,]/est_a[2,]

  mat.mxz = cbind(mat.sM,mat.xz)
  est_b = coxph(surv.Y~mat.mxz)
  est.beta = (coef(est_b)/sqrt(diag(est_b$var)))[1:size.s2]

  T2.part = est.alpha * est.beta

  return(list(
    T2_part = T2.part,
    id.S2 = in.id.S2,
    indv2 = in.indv2
  ))
}

#' @title Threshold Calculation
#' @description A function to calculate the threshold \eqn{L} for Mirror statistics.
#' @param Wk vector, the value of mirror statistics \eqn{W_k}.
#' @param g.q numeric, from 0 to 1. The target FDR level.
#' @return A threshold value \eqn{L} according to specified target FDR level \eqn{q}.
#' @section Details:
#' The threshold \eqn{L} is calculated as follow:
#' \deqn{L = \inf \left\{ t: \frac{ |\{k: W_k \le -t\}| }{ \max(|\{k: W_k \ge t\}|,1) } \le q , k \in \mathcal{S}_2 \right\}}
#' @export
get.MISAthreshold = function(Wk,g.q){ # Wk should be vector or 1-col matrix
  Wk_abs = abs(Wk)
  n = length(Wk)

  ord = order(Wk_abs,decreasing = T)
  Wk = Wk[ord]
  Wk_abs = Wk_abs[ord]

  sign_w = Wk>=0
  t_ = cumsum(sign_w)
  # fdp = 1-(t_/1:n)  ### %v 似乎之前算错了
  fdp = ((1:n)/t_)-1
  fdp = fdp[ Wk_abs > 0]

  if(length(which(fdp <= g.q))==0){
    thre = max(Wk_abs) + 1
    message("No t satisfy the FDP condition!")
  }else{
    thre = Wk_abs[max(which(fdp<=g.q))] ## %v 离散修正 - --
  }

  return(thre)
}

#' @title Step4: Aggregation
#' @description Aggregate the \eqn{T_1} and \eqn{T_2} derived in step 3 for inference. First construct \eqn{W_k} and calculate the threshold \eqn{L}. Then select
#' \deqn{\hat{\mathcal{H}}_1 = \{k : W_k \ge L\}.}
#'
#' @param in.T1D1_S2 numeric vector, contains \eqn{T_1^{(k)}} for all \eqn{k} in \eqn{\mathcal{S}_2}.
#' @param in.T2D2_S2 numeric vector, contains \eqn{T_2^{(k)}} for all \eqn{k} in \eqn{\mathcal{S}_2}.
#' @param in.id.S2 vector, indices of the mediators retained by the penalized Cox model, that is \eqn{\mathcal{S}_2}.
#' @param para.q numeric value from 0 to 1. The target FDR level.
#' @return `id.S3`, the finally selected mediators, that is \eqn{\hat{\mathcal{H}}_1}.
#' @export
Aggr_step = function(
    in.T1D1_S2
    ,in.T2D2_S2
    ,in.id.S2
    ,para.q
){
  Wk = in.T1D1_S2 * in.T2D2_S2
  thre = get.MISAthreshold(Wk,para.q)
  id.rej = which(Wk >= thre)
  id.S3 = in.id.S2[id.rej]
  return(id.S3)
}


#' @title Derandomization Operation
#' @description The function to aggregate the results from multiple sample splitting.
#' First we generate the voting set and choose the set which has max overlap with voting set.
#' Voting set is defined as
#' \deqn{ \hat{\mathcal{H}}_1^{(vote)} = \left\{ k:\sum_{b = 1}^{B}{ I(k \in \hat{\mathcal{H}}_1^{(b)}) } \ge B \times freq \right\}.}
#' @param in.S3list list, with `B` elements. `b`-th element is the set of selected mediators derived from the `b`-th sample splitting.
#' @param freq numeric value from `0` to `1`, default as `0.5`. The frequency threshold which is used to determine whether an index should be included in the voting set.
#' @return
#' A list with elements `id.S3.deran` , `B` , `freq` , `S_vote` , `bstar` :
#'
#' `id.S3.deran`  the indices of mediators selected by derandomization operation.
#'
#' `B`  the total times of sample splitting.
#'
#' `freq`  a table which is used to illustrate the frequency of each mediator present in results of multiple sample splitting.
#'
#' `S_vote`  the voting set.
#'
#' `bstar`  the indice of which set has the max overlap with the voting set.
#' @export
Deran_step = function(in.S3list,para.freq = 0.5){
  # in.S3list = s3list
  B = length(in.S3list)
  tab.freq = table(unlist(in.S3list))/B
  S.vote = as.numeric(names(tab.freq))[which(tab.freq >= para.freq)]
  r.overlap = foreach(id.S3 = in.S3list,.combine = "c")%do%{
    # c(sum(id.S3 %in% S.vote),length(id.S3))
    c(sum(id.S3 %in% S.vote)/(length(id.S3)+1),length(id.S3))
  }
  r.overlap = matrix(r.overlap,ncol = 2,byrow = T)
  ind.maxol = which(r.overlap[,1] == max(r.overlap[,1]))
  bstar = ind.maxol[which.min(r.overlap[ind.maxol,2])]
  # ind.maxol = which(r.overlap[,1] == max(r.overlap[,1]))
  # bstar = ind.maxol[which.min(abs(r.overlap[ind.maxol,2] - mean(r.overlap[ind.maxol,2])))]
  # in.S3list[[2]]

  toret = list(
    "id.S3.deran" = in.S3list[[bstar]],
    "B" = B,
    "freq" = tab.freq,
    "S_vote" = S.vote,
    "bstar" = bstar
  )
  return(toret)
}

#' @title All in One function
#' @description Complete the survMISA procedure all in one function. It includes all the steps of survMISA.
#' @param in.X vector, with `n` elements or matirx with `1` column, the input of one-dimensional exposures.
#' @param in.Z matrix, with `n` rows and `s` colnumns, the input of s-dimensional covariates.
#' @param in.M matrix, with `n` rows and `p` colnumns, the input of p-dimensional candidate mediators.
#' @param in.Y `Surv` object, with `n` elements, the input of survival information.
#' @param para.kappa numeric, suggested to be `1` or `2`. When `para.d` is not specified, it will be used to calculate `para.d` as \eqn{d = n/log(n) * \kappa}. Refer to `PreScr_step()`
#' @param para.nu numeric, to calculate the size of subsamples. Refer to `SampleSplit_step`.
#' @param para.seed numeric or vector, random seed for random splitting. A vector of seeds means multiple sample splitting with each seed. Refer to `SampleSplit_step()`.
#' @param para.rs one of `half`, `BIC` and `1se`. The strategy for selecting the penalty parameter \eqn{\lambda} in the penalized Cox model. Refer to `T1onD1()`.
#' @param para.q numeric value from 0 to 1. The target FDR level.
#' @param para.parallel logical, default as `FALSE`. Determines whether to use parallel computation for acceleration.
#' @param para.coresnum integer, the number of cores to use when `para.parallel` is set to `TRUE`. Default as all the cores detected by `detectCores()`.
#' @param para.d integer, the number of mediators will be retained by SIS program. If `para.d` is not specified, we will use `para.kappa` and the sample size `n` to calculate it. Refer to `PreScr_step()`.
#' @param in.stat_prescr vector with length p, the user-specified filter indicator. Refer to `PreScr_step()`.
#' @return A list with elements `deran` , `S3.mulsplit` , `all.mulsplit` :
#'
#' `deran`  the result of derandomization operation. If the length of `para.seed` is 1, than it return the result from single sample splitting.
#'
#' `S3.mulsplit` the list of `B` sets. Each of them includes the indices of mediators selected by the `b`-th sample splitting.
#'
#' `all.mulsplit`  the list of `B` list. Each of them includes detailed information derived in the `b`-th sample splitting.
#'
#' @examples
#' \code{n = 500;p = 1000;
#' A = matrix(c(rep(1,20),rep(-1,20),rep(0,p-40)),nrow = 1,ncol = p)/2
#' eta = matrix(c(0.5,-0.5),nrow = 2,ncol = p);gam = 0.5;xi = matrix(c(0.5,0.5),2,1)
#' B = matrix(c(rep(1,20),rep(0,20),rep(-1,20),rep(0,940)),nrow = p,ncol = 1)/2

#' set.seed(1)
#' # Exposures and covariates
#' X = matrix(rnorm(n),nrow = n)
#' Z = matrix(rnorm(n*2),n,2)
#' # Candidate mediators
#' eM = matrix(rnorm(n*p),n,p)
#' M = X %*% A + Z %*% eta + eM
#' # Survival outcomes
#' expp = exp(X %*% gam + Z %*% xi + M %*% B)
#' DT = sqrt(log(1-matrix(runif(n),ncol = 1))/expp*(-1) * 2)# Death time
#' CT = runif(n,max = 10) # Censoring time
#' Y = pmin(DT,CT)
#' Delta = (Y==DT)*1

#' # check the model
#' require(survival)
#' coxph(Surv(Y,Delta)~X+Z+M[,c(1:20,41:60)])
#' library(survMISA)
#' atest = AllinOne(
#'   in.X = X,
#'   in.Z = Z,
#'   in.M = M,
#'   in.Y = Surv(Y,Delta),
#'   para.kappa = 1,
#'   para.nu = 2,
#'   para.seed = 1:2,
#'   para.rs = "half",
#'   para.q = 0.1,
#'   para.parallel = F,
#'   para.freq = 0.5)
#' atest$deran$id.S3.deran ## The real mediators is 1,2,...,20}
#' @export
AllinOne = function(
    in.X
    ,in.Z
    ,in.M
    ,in.Y
    ,para.kappa = 1 # 如果未指定 para.d ，则使用 para.kappa 计算 para.d
    ,para.nu  = 1 # 分割比例，nu:1 = n1:n2
    ,para.seed = 1 # 样本分割时的随机数种子；如果长度为1，则使用survMISA；若长度大于1，则使用survMISA-DR
    ,para.rs = c("half","BIC","1se") # 使用惩罚估计筛选变量的策略。
    ,para.q ## 预期 FDR 水平
    ,para.parallel = F, para.coresnum # 是否使用并行运算
    ,para.freq = 0.5 # 去随机化步骤中的超参数
    ,para.d # 预计在预筛选步骤中保留的个数，建议小于中介变量的个数，否则预筛选不会发生。
    ,in.stat_prescr ## 用户指定的预筛选指标序列，绝对值越大将越有可能被保留
){
  if(para.parallel){check.parallel()}
  if(!is.matrix(in.M)){
    stop("Mediators should has a form of matrix, the row is an individual and the column is a mediator.")
  }
  n = nrow(in.M)
  p = ncol(in.M)

  mat.X = matrix(in.X,nrow = n)
  # if(!missing(in.Z)){mat.Z = matrix(in.Z,nrow = n)}

  mat.M = in.M
  surv.Y = in.Y

  message("Procedure start: Step 1 precreening")
  id.prescr = PreScr_step(
    in.X = mat.X,
    in.Z = in.Z,
    in.M = mat.M,
    in.stat_prescr = in.stat_prescr,
    para.kappa = para.kappa,
    para.d = para.d,
    para.parallel = para.parallel,
    para.coresnum = para.coresnum
  )

  message("Step1 finish, step2 sample splitting or multiple sample splitting start.")
  multiple_split = list()
  id.S3.list = list()
  ct = 1
  for(para.iseed in para.seed){
    tomess = paste0("Split :: ",ct,";Random seed :: ",para.iseed)
    message(tomess)
    indv.ss = SampleSplit_step( # individual index in sample splitting
      in.n = n,
      para.nu = para.nu,
      para.seed = para.iseed)

    indv1 = indv.ss$indv1
    indv2 = indv.ss$indv2

    if(length(para.seed) == 1){
      tomes = "Sample splitting finish, step 3.1 start"
    }else{
      tomes = paste0("Sample splitting finish, step 3.1 start for ",ct," / ",length(para.seed))
    }
    message(tomes)
    t1d1.list = T1onD1(
      in.X = mat.X,
      in.Z = in.Z,
      in.M = mat.M,
      in.Y = surv.Y,
      in.id.S1 = id.prescr,
      in.indv1 = indv1,
      para.rs = para.rs,
      para.parallel = para.parallel,
      para.coresnum = para.coresnum
    )

    if(length(para.seed) == 1){
      tomes = "Step3.1 finish, step 3.2 start"
    }else{
      tomes = paste0("Step 3.1 finish, step 3.2 start for ",ct," / ",length(para.seed))
    }
    message(tomes)
    t2d2.list = T2onD2(
      in.X = mat.X,
      in.Z = in.Z,
      in.M = mat.M,
      in.Y = surv.Y,
      in.id.S2 = t1d1.list$id.S2,
      in.indv2 = indv2,
      para.parallel = para.parallel,
      para.coresnum = para.coresnum
    )

    if(length(para.seed) == 1){
      tomes = "Step3.2 finish, step 4 start"
    }else{
      tomes = paste0("Step 3.2 finish, step 4 start for ",ct," / ",length(para.seed))
    }
    message(tomes)
    id.S3 = Aggr_step(
      in.T1D1_S2 = t1d1.list$T1_part,
      in.T2D2_S2 = t2d2.list$T2_part,
      in.id.S2 = t1d1.list$id.S2,
      para.q = para.q
    )
    single_split = list(
      id.S1 = id.prescr,
      id.S2 = t1d1.list$id.S2,
      id.S3 = id.S3,
      seed = para.iseed,
      indv1 = indv1,
      indv2 = indv2,
      Wk.S2 = t1d1.list$T1_part * t2d2.list$T2_part
    )
    multiple_split[[ct]] = single_split
    id.S3.list[[ct]] = id.S3
    ct = ct+1
  }

  list.deran = Deran_step(id.S3.list,para.freq = para.freq)


  toret = list(
    deran = list.deran,
    S3.mulsplit = id.S3.list,
    all.mulsplit = multiple_split
  )
  return(toret)
}
