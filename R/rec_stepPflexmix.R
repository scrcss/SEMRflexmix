
rec_stepPFlexmix <- function(x, y, z, min.comp=10, max.coef = 10, try.times = 5, Wupdate = 1)
{
  N = length(y)
  W = rep(1, N)
  k = 1
  res = list()
  Model <- FLXMRParticalglm(as.formula(paste0('~', paste(paste0("x", 1:dim(x)[2]), collapse = '+'))), Wupdate = Wupdate)
  if(is.null(z)){
    conModel <- NULL 
    dataFrame <-data.frame(y=y,x)
    colnames(dataFrame) <- c('y', paste0("x", 1:dim(x)[2]))
  }else{
    conModel <- FLXPmultinom(as.formula(paste0('~', paste(paste0("z", 1:dim(z)[2]), collapse = '+'))))
    dataFrame <-data.frame(y=y,x,z)
    colnames(dataFrame) <- c('y', paste0("x", 1:dim(x)[2]), paste0("z", 1:dim(z)[2]))
  }
  mycont <- list(minprior = 1)
  as(mycont, "FLXcontrol")
  while(sum(W) > min.comp){
    singular = 1 
    try.init = 0
    invalid.est = T
    flag  = F
    flag.err = F
    while((singular!=0)|invalid.est){
      try.init = try.init + 1
      res[[k]] = try(stepPFlexmix(y~1,k = 2,model=Model,concomitant = conModel,data=dataFrame, weights = W, control = mycont, nrep = 3, unique = T))   # 若遇到错误呢
      if("try-error" %in% class(res[[k]])){
        flag.err = T
        res[[k]] = NULL
        break
      }else{
        singular = res[[k]]@model[[1]]@singular
        invalid.est = any(abs(res[[k]]@concomitant@coef) > max.coef)
      }
      if(try.init > try.times) {
        flag = T
        break
      }
    }
    if(((flag == T) & (singular)) | (flag.err)) break
    W=W*res[[k]]@posterior$scaled[, 1]
    k = k+1
  }
  em.beta = sapply(res,function(x) parameters(x)$Comp.2)
  p = dim(em.beta)[1]
  beta = cbind(em.beta[1:(p-1),])
  sigma = em.beta[p,]
  U = cbind(sapply(res,function(x) x@posterior$scaled[,2]))
  final.U = t(rbind(1,apply(1-U,1,cumprod)))
  final.U = final.U*cbind(U,1)
  cluster=apply(final.U,1,which.max)
  singular = sapply(res, function(x) x@model[[1]]@singular)
  if (is.null(z)){
    temp_v = sapply(res, function(x) x@concomitant@coef[, 2])
    v = cumprod(c(1, 1-temp_v)) * c(temp_v,1)
  }else{
    v = sapply(res, function(x) x@concomitant@coef[, 2])
  }
  
  return(list(beta=beta,sigma=sigma,v = v, U=final.U,cluster=cluster,k=k,err = flag.err, singular = singular))
}
