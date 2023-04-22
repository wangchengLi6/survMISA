est.linear = function(x,y,para.parallel = F,coresnum){
  get.linear = function(x,y){
    mls = summary(lm(y~x))
    mc = mls$coefficients
    return(c(mc[2,1],mc[2,2]))
  }

  if(para.parallel){check.parallel()}

  if(para.parallel){
    if(missing(coresnum)){
      cl = makeCluster(detectCores())
    }else{
      cl = makeCluster(coresnum)
    }

    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    regOutput =
      foreach(yi = iter(y,by = "column")
              # ,.export = c("get.linear")
              ,.combine = "cbind") %dopar% {
                get.linear(x,yi)
              }
  }else{
    regOutput =
      foreach(yi = iter(y,by = "column")
              # ,.export = c("get.linear")
              ,.combine = "cbind") %do% {
                get.linear(x,yi)
              }
  }
  return(regOutput)
}

check.parallel = function(){
  ## 当需要进行并行运算时，我们需要载入 parallel 和 doParallel 软件包。
  ## 这项检查将在para.parallel参数为 TRUE 时进行
  if(!requireNamespace("parallel", quietly = TRUE)){
    stop(
      "Package `parallel` must be installed to use parallel computation.",
      call. = FALSE
    )
  }
  if(!requireNamespace("doParallel", quietly = TRUE)){
    stop(
      "Package `doParallel` must be installed to use parallel computation.",
      call. = FALSE
    )
  }
  require(parallel)
  require(doParallel)
}

#' @export
gene_CoxTime = function(expp,a=1,b=0,inverse_int = NULL){
  # generate standard uniform variable
  n = length(expp)
  temp_uniform = matrix(runif(n),ncol = 1)

  # culculate inverse function to get survival time
  h0 = log(1-temp_uniform)/expp*(-1)

  # base hazard is au^b
  time_cox = (h0 * (b+1) / a)**(1/(b+1))

  return(time_cox)
}
